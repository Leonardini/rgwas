KLIST = c(rep(1:3, 3),        1, 5)
LLIST = c(rep(1:3, each = 3), 5, 1)

### Driver function for the discovery-validation pipeline, whose arguments parallel those of the mainDriver and similar functions.
### If shuffle = FALSE, the first half of the patients go into discovery and the second into validation; if TRUE, simulate an RCT.
validationDriver = function(inputFile, type = "CNF", objective = "agreement", complement = 2, extremeValue = log(MAX_P),
                            Klist = KLIST, Llist = LLIST, index = 0, outputAssociations = TRUE, shuffle = TRUE) {
  stopifnot(extremeValue <= 0)
  ext <- str_sub(inputFile, start = -4)
  if (str_sub(inputFile, -3) == ".gz") {
    ext <- str_sub(inputFile, -7, -4)
    miniFile <- str_sub(inputFile, 1, -4)
  } else {
    miniFile <- inputFile
  }
  readFunction <- ifelse(ext == '.tsv', read_tsv, read_csv)
  inputTab <- readFunction(inputFile, col_types = cols(ID = "c", id = "c", .default = ifelse(type == "LCP", "d", "l")))
  if ("id" %in% colnames(inputTab)) { inputTab %<>% rename(ID = id) }
  numSNPs <- str_extract(miniFile, paste0("([0-9]+)", ext)) %>% str_remove(ext) %>% parse_integer
  IDs <- inputTab %>% 
    dplyr::select(ID)
  myMatrices <- prepareMatrices(inputTab, numSNPs, complement = complement)
  phenotypes <- myMatrices$phenotypes %>%
    as_tibble()
  genotypes <- myMatrices$genotypes %>%
    as_tibble()
  SNPs <- colnames(genotypes)
  results <- vector("list", length(SNPs)) %>%
    set_names(SNPs)
  nPheno <- ncol(phenotypes)
  for (ind in 1:length(SNPs)) {
    curSNP <- SNPs[ind]
    targetOutFile <- str_replace(miniFile, ext, paste0(curSNP, "_", 1, "_", "FullResults", ext))
    if (file.exists(targetOutFile)) {       
      print(paste("Skipping SNP", curSNP))
    } else {
      redTab <- bind_cols(IDs, genotypes %>% dplyr::select(all_of(curSNP)), phenotypes)
      redFilename <- str_replace(miniFile, ext, paste0(curSNP, "_", 1, ext))
      write_csv(redTab, redFilename)
      results[[curSNP]] <- fullValidation(redFilename, type = type, objective = objective, extremeValue = extremeValue, 
                           Klist = Klist, Llist = Llist, index = index * THOUSAND + ind, nPheno = nPheno, shuffle = shuffle)
    }
  }
  if (outputAssociations) {
    singleRes <- computeBestPValues(myMatrices$genotypes, myMatrices$phenotypes, cor = (type == "LCP"), fisher = (type != "LCP"))
    write_csv(singleRes$allValues %>% as_tibble(rownames = "association"), str_replace(miniFile, ext, '_Associations.csv'))
    N <- nrow(myMatrices$genotypes)
    if (shuffle) {
      set.seed(MY_SEED)
      first <- which(runif(N) <= 0.5)
    } else {
      first <- 1:(round(N/2))
    }
    second <- setdiff(1:N, first)
    singleTrain <- computeBestPValues(myMatrices$genotypes[first, , drop = FALSE], myMatrices$phenotypes[first, , drop = FALSE], 
                                      cor = (type == "LCP"), fisher = (type != "LCP"))
    write_csv(singleTrain$allValues %>% as_tibble(rownames = "association"), str_replace(miniFile, ext, '_Train_Associations.csv'))
    singleTest <- computeBestPValues(myMatrices$genotypes[second, , drop = FALSE], myMatrices$phenotypes[second, , drop = FALSE], 
                                     cor = (type == "LCP"), fisher = (type != "LCP"))
    write_csv(singleTest$allValues  %>% as_tibble(rownames = "association"), str_replace(miniFile, ext,  '_Test_Associations.csv'))
  }
  results
}

fullValidation = function(inputFile, type, objective, extremeValue, Klist, Llist, index = 0, removeInput = TRUE, nPheno = 0, shuffle = FALSE) {
  stopifnot(extremeValue <= 0)
  splitNames <- splitFile(inputFile, type = type, shuffle = shuffle)
  stopifnot(length(Klist) == length(Llist))
  initRes <- tibble()
  valRes  <- tibble()
  for (ind in 1:length(Klist)) {
    curK <- Klist[ind]
    curL <- Llist[ind]
    curResult <- optimizingDriver(splitNames[1], startPValue = extremeValue, factor = 2, type = type, complement = 0, 
                                  K = curK, L = curL, startIndex = index * THOUSAND + ind)
    curSummary <- curResult$summary
    if (!is.null(curSummary)) {
      initRes %<>% bind_rows(bind_cols(tibble(K = curK, L = curL), curSummary))
      if (curSummary$status[1] %in% c("feasible", "optimal")) {
        ext <- str_sub(splitNames[1], -4)
        specString <- paste(c("", type, "K", curK, "L", curL), collapse = "_")
        write_csv(curResult$summary,                  str_replace(splitNames[1], ext, str_c(specString,      '_Summary.csv')))
        write_csv(as.data.frame(curResult$phenotype), str_replace(splitNames[1], ext, str_c(specString,    '_Phenotype.csv')))
        curFormula <- curSummary$formula
        BD <- ifelse("LogCMHExactP" %in% names(curSummary), curSummary$bestSingle, NA)
        valResults <- generateValidationPhenotype(splitNames[2], curFormula, type = type, K = curK, L = curL, bestDiscovery = BD,
                                                  outputSummary = FALSE, outputPhenotype = FALSE, outputAssociations = FALSE)
        valRes %<>%
          bind_rows(bind_cols(tibble(K = curK, L = curL), valResults$summary))
      }
    }
  }
  if (nrow(valRes) > 0) {
    fullRes <- left_join(initRes, valRes, by = c("K", "L", "SNP"), suffix = c("_Discovery", "_Validation"))
  } else {
    fullRes <- initRes
  }
  write_csv(fullRes, str_replace(inputFile, paste0("\\", str_sub(inputFile, -4)), '_FullResults.csv'))
  if (removeInput) { file.remove(c(splitNames, inputFile)) }
  fullRes
}

splitFile = function(inputFile, ext = paste0("\\", str_sub(inputFile, -4)), type = "LCP", shuffle = FALSE) {
  readFunction <- ifelse(ext == '.tsv', read_tsv, read_csv)
  inputTab <- readFunction(inputFile, col_types = cols(ID = "c", .default = ifelse(type == "LCP", "d", "l")))
  N <- nrow(inputTab)
  if (shuffle) {
    set.seed(MY_SEED)
    firstHalf <- which(runif(N) <= 0.5)
  } else {
    firstHalf <- 1:(round(N/2))
  }
  secondHalf <- setdiff(1:N, firstHalf) 
  Tab1 <- inputTab %>% slice(firstHalf)
  Tab2 <- inputTab %>% slice(secondHalf)
  inputFile1 <- str_replace(inputFile, ext, '_Train_1.csv')
  inputFile2 <- str_replace(inputFile, ext,  '_Test_1.csv')
  stopifnot(!file.exists(inputFile1) && !file.exists(inputFile2))
  write_csv(Tab1, inputFile1)
  write_csv(Tab2, inputFile2)
  output <- c(inputFile1, inputFile2)
  output
}

parseFormula = function(Formula, type) {
  outerString   <- ifelse(type == "LCP", " \\+ ", ifelse(type == "CNF", " AND ", " OR "))
  innerString   <- ifelse(type == "LCP", " \\* ", ifelse(type == "DNF", " AND ", " OR "))
  splitFormula <- str_split(Formula, outerString)[[1]]
  fullySplitFormula <- map(splitFormula, function(x) { 
    x %>% str_trim %>% str_remove("^\\(") %>% str_remove("\\)$") %>% str_split(innerString) %>% extract2(1) })
  extractFunction <- ifelse(type == "LCP", function(x) {x[[2]]},                                         function(x) {x})
  extractCoeffs   <- ifelse(type == "LCP", function(x) {x[[1]] %>% str_remove("\\)$") %>% parse_double}, function(x) {1})  
  usedNames  <- map(fullySplitFormula, extractFunction)
  usedValues <- map(fullySplitFormula, extractCoeffs) 
  allUsedNames <- sort(unique(unlist(usedNames)))
  M <- length(allUsedNames)
  N <- ifelse(type == "LCP", 1, length(splitFormula))
  parsedFormula <- matrix(0, M, N, dimnames = list(allUsedNames, 1:N))
  if (type == "LCP") {
    parsedFormula[unlist(usedNames), 1] <- unlist(usedValues)
  } else {
    for (ind in 1:N) {
      parsedFormula[usedNames[[ind]], ind] <- usedValues[[ind]]
    }
  }
  parsedFormula
}

applyFormula = function(parsedFormula, phenotypes, type) {
  stopifnot(all(rownames(parsedFormula) %in% colnames(phenotypes)))
  if (type == "LCP") {
    finalOutput <- phenotypes[, rownames(parsedFormula), drop = FALSE] %*% parsedFormula
  } else {
    phenotypes %<>%
      as.data.frame %>%
      dplyr::select(rownames(parsedFormula))
    outerFunction <- ifelse(type == "CNF", and, or)
    innerFunction <- ifelse(type == "DNF", and, or)
    allTerms <- map(1:ncol(parsedFormula), ~{reduce(phenotypes[, (parsedFormula[, .] == 1), drop = FALSE], innerFunction)})
    finalOutput   <- reduce(allTerms, outerFunction)
  }
  if (class(finalOutput) != "matrix") { finalOutput <- matrix(finalOutput, ncol = 1) }
  finalOutput
}

generateValidationPhenotype = function(inputFile, Formula, type, objective = ifelse(type == "LCP", "correlation", "agreement"),
                      K = 3, L = 3, bestDiscovery = NA, outputSummary = FALSE, outputPhenotype = FALSE, outputAssociations = FALSE, 
                      SNP = "") {
  ext <- str_sub(inputFile, -4)
  if (str_sub(inputFile, -3) == ".gz") {
    ext <- str_sub(inputFile, -7, -4)
    miniFile <- str_sub(inputFile, 1, -4)
  } else {
    miniFile <- inputFile
  }
  numSNPs <- str_extract(miniFile, paste0("([0-9]+)", ext)) %>% str_remove(ext) %>% parse_integer
  readFunction <- ifelse(ext == '.tsv', read_tsv, read_csv)
  inputTab <- readFunction(inputFile, col_types = cols(ID = "c", .default = ifelse(type == "LCP", "d", "l")))
  myMatrices <- prepareMatrices(inputTab, numSNPs, complement = -1) ### previously was 0
  phenotypes <- myMatrices$phenotypes
  parsedFormula <- parseFormula(Formula, type = type)
  outputPheno <- applyFormula(parsedFormula, phenotypes, type = type)
  genotypes <- myMatrices$genotypes
  singleRes <- computeBestPValues(genotypes, phenotypes, cor = (type == "LCP"), fisher = (type != "LCP"))
  curValues <- rep(NA, numSNPs)
  if (type == "LCP") {
    stopifnot(objective == "correlation")
    curValues <- as.vector(computeCorr(genotypes, outputPheno, scaleUp = FALSE))
  } else {
    if (objective == "agreement") {
      curValues <- as.vector(computeAgreements(genotypes, outputPheno))
    } else {
      curValues <- as.vector(computeCovariances(genotypes, outputPheno, scaleUp = FALSE))
    }
  }
  outputS <- vector("list", numSNPs) %>%
    set_names(colnames(genotypes))
  MH <- rep(NA, numSNPs)
  for (ind in 1:numSNPs) {
    curName <- colnames(genotypes)[ind]
    curValue <- curValues[ind]
    curOutput <- tibble(time = NA, gap = NA, opt = curValue, status = "optimal", formula = Formula)
    curGeno <- genotypes[, ind]
    if (type != "LCP") {
      myT <- matrix(0, 2, 2, dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE")))
      extraT <- table(as.logical(curGeno), as.logical(outputPheno))
      myT[rownames(extraT), colnames(extraT)] %<>%
        add(extraT)
      checkR <- tibble(TP = myT["TRUE", "TRUE"], FP = myT["FALSE","TRUE"], FN = myT["TRUE","FALSE"], TN = myT["FALSE","FALSE"])
      curOutput %<>% 
        bind_cols(TP = checkR$TP, FP = checkR$FP, FN = checkR$FN, TN = checkR$TN)
      if (!is.na(bestDiscovery)) {
        bestSingle = as.logical(phenotypes[, bestDiscovery, drop = TRUE])
        if (sum(bestSingle) > 1 && sum(!bestSingle) > 1) { ### if the sum is wrong, a strata is too small and the test fails! 
          MH[ind] <- computeStratifiedPVal(as.logical(curGeno), as.logical(outputPheno), bestSingle)
        }
      }
    } else {
      curPValue <- cor.test(curGeno, outputPheno, method = "pearson", alternative = "g")$p.value
      curOutput %<>%
        bind_cols(pValue = curPValue)
    }
    outputS[[curName]] <- curOutput
  }
  outputSum <- matrix(unlist(outputS), nrow = length(outputS), dimnames = list(NULL, names(outputS[[1]])), byrow = TRUE) %>%
    as_tibble %>%
    mutate_at(c("time", "gap", "opt"), as.double) %>%
    mutate(SNP = colnames(genotypes)) %>%
    dplyr::select(SNP, everything())
  if (type != "LCP") {
    outputSum %<>%
      mutate_at(c("TP", "FP", "TN", "FN"), as.integer) %>%
      computeStats()
    if (!all(is.na(MH))) {
      outputSum %<>%
        mutate(LogCMHExactP = MH)
    }
  }
  outputSum %<>%
    mutate(bestSingle = singleRes$phenotype, bestSingleStat = singleRes$value)
  complexPhenotypes = tibble(phenotype = outputPheno) %>%
    set_colnames(colnames(genotypes)) %>%
    mutate_all(as.numeric) %>%
    mutate(ID = inputTab$ID) %>%
    select(ID, everything())
  output <- list(phenotype = complexPhenotypes, summary = outputSum, associations = singleRes$allValues %>% 
                   as_tibble(rownames = "association"))
  specString <- paste0(paste(c("", type, "K", K, "L", L), collapse = "_"), ifelse(SNP != "", paste0("_", SNP), ""))
  if (outputSummary)      {write_csv(output$summary,      str_replace(miniFile, ext, str_c(specString,      '_Summary.csv')))}
  if (outputPhenotype)    {write_csv(output$phenotype,    str_replace(miniFile, ext, str_c(specString,    '_Phenotype.csv')))}
  if (outputAssociations) {write_csv(output$associations, str_replace(miniFile, ext, str_c(specString, '_Associations.csv')))}
  output
}
