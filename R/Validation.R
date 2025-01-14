#' Run `reverseGWAS` on a discovery and validation split of the data.
#'
#' Carry out the discovery-validation pipeline, whose arguments parallel those of the mainDriver and similar functions.
#'
#' @inheritParams mainDriver
#' @param Klist,Llist Integer vectors of values for `K` (number of clauses) and `L` (size of clauses); they must have the same length!
#' @param shuffle Perform an RCT? If `TRUE` (the default), simulate an RCT; otherwise, the first half of the patients (input rows) go into discovery and the second half, into validation.
#' @param SEED Integer seed for the random number generator; keeping it consistent between runs generates the same RCT patient split, thus ensuring reproducible results; default: `987654321`.
#'
#' @examples
#' # Y <- rgwas::validationDriver(inputFile = "TestInputN500000P10_3.csv.gz", extremeValue = log(5e-8), shuffle = TRUE)
#'
#' @export
validationDriver = function(inputFile, type = "CNF", objective = "agreement", complement = FALSE, extremeValue = log(MAX_P),
                            Klist = KLIST, Llist = LLIST, index = 0, outputAssociations = FALSE, shuffle = TRUE, SEED = 987654321L) {
  stopifnot(extremeValue <= 0)
  ext <- stringr::str_sub(inputFile, start = -4)
  if (stringr::str_sub(inputFile, -3) == ".gz") {
    ext <- stringr::str_sub(inputFile, -7, -4)
    miniFile <- stringr::str_sub(inputFile, 1, -4)
  } else {
    miniFile <- inputFile
  }
  readFunction <- ifelse(ext == '.tsv', readr::read_tsv, readr::read_csv)
  inputTab <- readFunction(inputFile, col_types = readr::cols(ID = "c", .default = "l"))
  numSNPs <- stringr::str_extract(miniFile, paste0("([0-9]+)", ext)) %>%
    stringr::str_remove(ext) %>%
    readr::parse_integer()
  IDs <- inputTab %>%
    dplyr::select(ID)
  myMatrices <- prepareMatrices(inputTab, numSNPs, complement = complement)
  phenotypes <- myMatrices$phenotypes %>%
    tibble::as_tibble()
  genotypes <- myMatrices$genotypes %>%
    tibble::as_tibble()
  SNPs <- colnames(genotypes)
  results <- vector("list", length(SNPs)) %>%
    set_names(SNPs)
  nPheno <- ncol(phenotypes)
  for (ind in 1:length(SNPs)) {
    curSNP <- SNPs[ind]
    targetOutFile <- stringr::str_replace(miniFile, ext, paste0(curSNP, "_", 1, "_", "FullResults", ext))
    if (file.exists(targetOutFile)) {
      print(paste("Skipping SNP", curSNP))
    } else {
      redTab <- dplyr::bind_cols(IDs, genotypes %>% dplyr::select(all_of(curSNP)), phenotypes)
      redFilename <- stringr::str_replace(miniFile, ext, paste0(curSNP, "_", 1, ext))
      readr::write_csv(redTab, redFilename)
      results[[curSNP]] <- fullValidation(redFilename, type = type, objective = objective, extremeValue = extremeValue,
                           Klist = Klist, Llist = Llist, index = index * THOUSAND + ind, nPheno = nPheno, shuffle = shuffle, seed = SEED)
    }
  }
  if (outputAssociations) {
    singleRes <- computeBestPValues(myMatrices$genotypes, myMatrices$phenotypes, fisher = TRUE)
    readr::write_csv(singleRes$allValues %>%
                tibble::as_tibble(rownames = "association"), stringr::str_replace(miniFile, ext, '_Associations.csv'))
    N <- nrow(myMatrices$genotypes)
    if (shuffle) {
      set.seed(SEED)
      first <- which(runif(N) <= 0.5)
    } else {
      first <- 1:(round(N/2))
    }
    second <- setdiff(1:N, first)
    singleTrain <- computeBestPValues(myMatrices$genotypes[first, , drop = FALSE], myMatrices$phenotypes[first, , drop = FALSE],
                                      fisher = TRUE)
    readr::write_csv(singleTrain$allValues %>%
                tibble::as_tibble(rownames = "association"), stringr::str_replace(miniFile, ext, '_Train_Associations.csv'))
    singleTest <- computeBestPValues(myMatrices$genotypes[second, , drop = FALSE], myMatrices$phenotypes[second, , drop = FALSE],
                                      fisher = TRUE)
    readr::write_csv(singleTest$allValues  %>%
                tibble::as_tibble(rownames = "association"), stringr::str_replace(miniFile, ext,  '_Test_Associations.csv'))
  }
  results
}

fullValidation = function(inputFile, type, objective, extremeValue, Klist, Llist, index = 0, removeInput = TRUE, nPheno = 0,
                          shuffle = FALSE, seed) {
  stopifnot(extremeValue <= 0)
  splitNames <- splitFile(inputFile, type = type, shuffle = shuffle, seed = seed)
  stopifnot(length(Klist) == length(Llist))
  initRes <- tibble::tibble()
  valRes  <- tibble::tibble()
  print(paste("There are", length(Klist), "options to process"))
  for (ind in 1:length(Klist)) {
    print(paste("Processing option", ind, "of", length(Klist)))
    curK <- Klist[ind]
    curL <- Llist[ind]
    curResult <- optimizingDriver(splitNames[1], startPValue = extremeValue, factor = 2, type = type, complement = 0,
                                  K = curK, L = curL, startIndex = index * THOUSAND + ind)
    curSummary <- curResult$summary
    if (!is.null(curSummary)) {
      initRes <- initRes %>%
        dplyr::bind_rows(dplyr::bind_cols(tibble::tibble(K = curK, L = curL), curSummary))
      if (curSummary$status[1] %in% c("feasible", "optimal")) {
        ext <- stringr::str_sub(splitNames[1], -4)
        specString <- paste(c("", type, "K", curK, "L", curL), collapse = "_")
        readr::write_csv(curResult$summary,                  stringr::str_replace(splitNames[1], ext, stringr::str_c(specString,  '_Summary.csv')))
        readr::write_csv(as.data.frame(curResult$phenotype), stringr::str_replace(splitNames[1], ext, stringr::str_c(specString,  '_Phenotype.csv')))
        curFormula <- curSummary$formula
        BD <- ifelse("LogCMHExactP" %in% names(curSummary), curSummary$bestSingle, NA)
        valResults <- generateValidationPhenotype(splitNames[2], curFormula, type = type, K = curK, L = curL, bestDiscovery = BD,
                                                  outputSummary = FALSE, outputPhenotype = FALSE, outputAssociations = FALSE)
        valRes <- valRes %>%
          dplyr::bind_rows(dplyr::bind_cols(tibble::tibble(K = curK, L = curL), valResults$summary))
      }
    }
  }
  if (nrow(valRes) > 0) {
    fullRes <- left_join(initRes, valRes, by = c("K", "L", "SNP"), suffix = c("_Discovery", "_Validation"))
  } else {
    fullRes <- initRes
  }
  readr::write_csv(fullRes, stringr::str_replace(inputFile, paste0("\\", stringr::str_sub(inputFile, -4)), '_FullResults.csv'))
  if (removeInput) { file.remove(c(splitNames, inputFile)) }
  fullRes
}

splitFile = function(inputFile, ext = paste0("\\", stringr::str_sub(inputFile, -4)), type = "CNF", shuffle = FALSE, seed) {
  readFunction <- ifelse(ext == '.tsv', readr::read_tsv, readr::read_csv)
  inputTab <- readFunction(inputFile, col_types = readr::cols(ID = "c", .default = "l"))
  N <- nrow(inputTab)
  if (shuffle) {
    set.seed(seed)
    firstHalf <- which(runif(N) <= 0.5)
  } else {
    firstHalf <- 1:(round(N/2))
  }
  secondHalf <- setdiff(1:N, firstHalf)
  Tab1 <- inputTab %>%
    dplyr::slice(firstHalf)
  Tab2 <- inputTab %>%
    dplyr::slice(secondHalf)
  inputFile1 <- stringr::str_replace(inputFile, ext, '_Train_1.csv')
  inputFile2 <- stringr::str_replace(inputFile, ext,  '_Test_1.csv')
  stopifnot(!file.exists(inputFile1) && !file.exists(inputFile2))
  readr::write_csv(Tab1, inputFile1)
  readr::write_csv(Tab2, inputFile2)
  output <- c(inputFile1, inputFile2)
  output
}

generateValidationPhenotype = function(inputFile, Formula, type, objective = "agreement", K = 3, L = 3, bestDiscovery = NA,
                                       outputSummary = FALSE, outputPhenotype = FALSE, outputAssociations = FALSE, SNP = "") {
  ext <- stringr::str_sub(inputFile, -4)
  if (stringr::str_sub(inputFile, -3) == ".gz") {
    ext <- stringr::str_sub(inputFile, -7, -4)
    miniFile <- stringr::str_sub(inputFile, 1, -4)
  } else {
    miniFile <- inputFile
  }
  numSNPs <- stringr::str_extract(miniFile, paste0("([0-9]+)", ext)) %>%
    stringr::str_remove(ext) %>%
    readr::parse_integer()
  readFunction <- ifelse(ext == '.tsv', readr::read_tsv, readr::read_csv)
  inputTab <- readFunction(inputFile, col_types = readr::cols(ID = "c", .default = "l"))
  myMatrices <- prepareMatrices(inputTab, numSNPs, complement = -1)
  phenotypes <- myMatrices$phenotypes
  parsedFormula <- parseFormula(Formula, type = type)
  outputPheno <- applyFormula(parsedFormula, phenotypes, type = type)
  genotypes <- myMatrices$genotypes
  singleRes <- computeBestPValues(genotypes, phenotypes, fisher = TRUE)
  curValues <- rep(NA, numSNPs)
  if (objective == "agreement") {
    curValues <- as.vector(computeAgreements(genotypes, outputPheno))
  } else {
    curValues <- as.vector(computeCovariances(genotypes, outputPheno, scaleUp = FALSE))
  }
  outputS <- vector("list", numSNPs) %>%
    set_names(colnames(genotypes))
  MH <- rep(NA, numSNPs)
  for (ind in 1:numSNPs) {
    curName <- colnames(genotypes)[ind]
    curValue <- curValues[ind]
    curOutput <- tibble::tibble(time = NA, gap = NA, opt = curValue, status = "optimal", formula = Formula)
    curGeno <- genotypes[, ind]
    myT <- matrix(0, 2, 2, dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE")))
    extraT <- table(as.logical(curGeno), as.logical(outputPheno))
    myT[rownames(extraT), colnames(extraT)] <- myT[rownames(extraT), colnames(extraT)] %>%
      magrittr::add(extraT)
    checkR <- tibble::tibble(TP = myT["TRUE", "TRUE"], FP = myT["FALSE","TRUE"], FN = myT["TRUE","FALSE"], TN = myT["FALSE","FALSE"])
    curOutput <- curOutput %>%
      dplyr::bind_cols(TP = checkR$TP, FP = checkR$FP, FN = checkR$FN, TN = checkR$TN)
    if (!is.na(bestDiscovery)) {
      bestSingle = as.logical(phenotypes[, bestDiscovery, drop = TRUE])
      if (sum(bestSingle) > 1 && sum(!bestSingle) > 1) { # if the sum is wrong, a strata is too small and the test fails!
        MH[ind] <- computeStratifiedPVal(as.logical(curGeno), as.logical(outputPheno), bestSingle)
      }
    }
    outputS[[curName]] <- curOutput
  }
  outputSum <- matrix(unlist(outputS), nrow = length(outputS), dimnames = list(NULL, names(outputS[[1]])), byrow = TRUE) %>%
    tibble::as_tibble() %>%
    dplyr::mutate_at(c("time", "gap", "opt"), as.double) %>%
    dplyr::mutate(SNP = colnames(genotypes)) %>%
    dplyr::select(SNP, everything()) %>%
    dplyr::mutate_at(c("TP", "FP", "TN", "FN"), as.integer) %>%
    computeStats()
  if (!all(is.na(MH))) {
    outputSum <- outputSum %>%
      dplyr::mutate(LogCMHExactP = MH)
  }
  outputSum <- outputSum %>%
    dplyr::mutate(bestSingle = singleRes$phenotype, bestSingleStat = singleRes$value)
  complexPhenotypes = tibble::tibble(phenotype = outputPheno) %>%
    magrittr::set_colnames(colnames(genotypes)) %>%
    dplyr::mutate_all(as.numeric) %>%
    dplyr::mutate(ID = inputTab$ID) %>%
    dplyr::select(ID, everything())
  output <- list(phenotype = complexPhenotypes, summary = outputSum, associations = singleRes$allValues %>%
                   tibble::as_tibble(rownames = "association"))
  specString <- paste0(paste(c("", type, "K", K, "L", L), collapse = "_"), ifelse(SNP != "", paste0("_", SNP), ""))
  if (outputSummary)      {readr::write_csv(output$summary,      stringr::str_replace(miniFile, ext, stringr::str_c(specString,      '_Summary.csv')))}
  if (outputPhenotype)    {readr::write_csv(output$phenotype,    stringr::str_replace(miniFile, ext, stringr::str_c(specString,    '_Phenotype.csv')))}
  if (outputAssociations) {readr::write_csv(output$associations, stringr::str_replace(miniFile, ext, stringr::str_c(specString, '_Associations.csv')))}
  output
}
