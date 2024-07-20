source("Settings.R")
source("Preprocess.R")
source("Utilities.R")
source("Validation.R")

## This function generates the specified number of randomizations of the phenotype matrix (KAPPA = 10)
randomizePhenotypes = function(fname = "CPLEXInputUKBB_all_unrelated_Autoimmune_11_HapMap300_chr3_subset11b_281591_149.csv.gz", 
                               numRandomizations = 10) {
  set.seed(MY_SEED)
  fullTab = read_csv(fname, guess_max = Inf, show_col_types = FALSE)
  if (!("ID" %in% colnames(fullTab))) {
    fullTab %<>% rename(ID = id)
  }
  numSNPs = str_extract_all(fname, "[0-9]+") %>%
    map_chr(dplyr::last) %>%
    as.integer()
  prepMat = prepareMatrices(fullTab, numSNPs = numSNPs, complement = 0, keepDuplicates = TRUE)
  allPheno = prepMat$phenotypes
  allPheno = allPheno[, setdiff(colnames(allPheno), "Autoimmune.pheno"), drop = FALSE]
  birewire.sampler.bipartite(allPheno, K = numRandomizations, path = "../", exact = TRUE, write.sparse = FALSE)
  setwd("../1/")
  allRandomPheno = vector("list", numRandomizations)
  for (index in 1:numRandomizations) {
    ## allRandomPheno[[index]] = randomizeBinaryMatrix(allPheno, keepRownames = TRUE, keepColnames = TRUE)
    allRandomPheno[[index]] = read.delim(paste0("network_", index), sep = " ", header = TRUE)
  }
  setwd("../")
  save(allRandomPheno, file = "RandomPhenotypes.RData")
  allRandomPheno
}

## This function generates all the triage runs; if validation = TRUE, they are also going to undergo a validation on half the data 
prepAllTriage = function(Kvec = 1:3, Lvec = 1:3, epsVec = c(0, 0.1, 0.2), patVec = as.integer(10^(3:5)), type = "DNF", 
                         binWidth = 0.05, extra = 49, validation = FALSE) {
  set.seed(MY_SEED)
  L1 = length(Kvec)
  L2 = length(Lvec)
  L3 = length(epsVec)
  L4 = length(patVec)
  e = new.env()
  fname = paste0("RandomSets", ifelse(validation, "Validation", ""), ".RData")
  if (!file.exists(fname)) {
    allSets = vector("list", L4)
    load("RandomPhenotypes.RData", e)
    randomPheno = get("allRandomPheno", envir = e)
    randomPheno = randomPheno[[1]]
    rSums = rowSums(randomPheno)
    for (ind in 1:L4) {
      curSize = patVec[ind]
      curSet = sample(names(rSums), size = curSize * ifelse(validation, 2L, 1L), replace = FALSE, prob = rSums + 1L) %>%
        as.integer()
      allSets[[ind]] = curSet
    }
    save(allSets, file = fname)
  } else {
    load(fname, e)
    allSets = get("allSets", envir = e)
  }
  for (ind in 1:L4) {
    if (!file.exists(paste0("BinnedGenotypes1Size", length(allSets[[ind]]), ".RData"))) {
      print(paste("Binning the genotypes for size", patVec[ind]))
      Z = binGenotypesByFrequency(directory = "CPLEXInput_HapMap300k_AutoImmune_Traits", pattern = "chr1_", Subset = allSets[[ind]])
    }
  }
  allOpts = hcube(c(L1, L2, L3, L4))
  allPvals = vector("list", nrow(allOpts))
  print("Processing the option settings")
  for (index in 1:nrow(allOpts)) {
    curRow = allOpts[index, ]
    curK = Kvec[curRow[1]]
    curL = Lvec[curRow[2]]
    curEps = epsVec[curRow[3]]
    curPat = patVec[curRow[4]]
    pvals = prepSingleTriage(K = curK, L = curL, epsilon = curEps, numPatients = curPat, type = type, binWidth = binWidth, 
                             extra = extra, validation = validation)
    allPvals[[index]] = pvals
    print(index)
  }
  save(allPvals, file = "PValues.RData")
  allPvals
}

## This function generates a single triage run; if validation = TRUE, it is also going to undergo a validation on half the data 
prepSingleTriage = function(K = 3, L = 3, epsilon = 0, numPatients = 1000, type = "DNF", binWidth = 0.05, 
                            extra = 1, validation = FALSE) {
  e = new.env()
  load("RandomPhenotypes.RData", e)
  randomPheno = get("allRandomPheno", envir = e)
  load(paste0("RandomSets", ifelse(validation, "Validation", ""), ".RData"), e)
  randomSets = get("allSets", envir = e)
  curSet = randomSets[[which(sapply(randomSets, length) == numPatients * ifelse(validation, 2L, 1L))]] 
  innerF = ifelse(type == "DNF", and, or)
  outerF = ifelse(type == "DNF", or, and)
  Q = length(randomPheno)
  pvals = rep(NA, Q)
  for (index in 1:Q) {
    specString = paste0("K", K, "_L", L, "_M", numPatients * ifelse(validation, 2L, 1L), "_eps", epsilon, "_", type, "_I", index - 1, "_N", extra + 1)
    if (!file.exists(paste0("TriageInput_", specString, ".csv"))) {
      curPheno = randomPheno[[index]]
      curPheno = curPheno[as.character(curSet), , drop = FALSE]
      P = ncol(curPheno)
      allTerms = matrix(NA, K, L)
      fullGeno = rep(1 - as.integer(type == "DNF"), numPatients)
      allTerms[1:K, 1:L] = sample(1:P, K * L, replace = TRUE)
      for (ind in 1:K) {
        curTerm = allTerms[ind, ]
        curGeno = rep(as.integer(type == "DNF"), numPatients)
        for (ind2 in 1:L) {
          curGeno = innerF(curGeno, curPheno[, curTerm[ind2]])
        }
        fullGeno = outerF(fullGeno, curGeno)
      }
      names(fullGeno) = NULL
      dimnames(allTerms) = list(paste0("K", 1:K), paste0("L", 1:L))
      numBins = as.integer(1/binWidth) + 1L
      binNumber = floor(sum(fullGeno) / (numPatients * ifelse(validation, 2L, 1L) * binWidth)) + 1
      load(paste0("BinnedGenotypes", as.integer(binNumber), "Size", as.integer(numPatients * ifelse(validation, 2L, 1L)), ".RData"))
      uSNPs = unique(allSNPs[,2])
      useSNPs = sample(uSNPs, extra, replace = FALSE)
      goodSNPs = allSNPs[allSNPs[,2] %in% useSNPs, , drop = FALSE]
      complement = as.integer(binNumber >= numBins/2)
      goodMat = matrix(complement, nrow = numPatients * ifelse(validation, 2L, 1L), ncol = extra, dimnames = list(curSet, useSNPs))
      mode(goodSNPs) = "character"
      goodMat[goodSNPs] = 1 - complement
      ## goodMat = randomizeBinaryMatrix(goodMat, keepRownames = FALSE, keepColnames = FALSE)
      goodMat = birewire.rewire.bipartite(goodMat, exact = TRUE)
      colnames(goodMat) = paste0("S", 1:extra)
      goodMat = as_tibble(goodMat)
      noisyGeno = fullGeno
      perturbPositions = c()
      extraPerturbation = c()
      if (epsilon > 0) {
        perturbPositions = sample(numPatients, epsilon * numPatients, replace = FALSE)
        noisyGeno[perturbPositions] = !(noisyGeno[perturbPositions])
        if (validation) {
          extraPerturbation = sample(numPatients, epsilon * numPatients, replace = FALSE)
          noisyGeno[numPatients + extraPerturbation] = !(noisyGeno[numPatients + extraPerturbation])
        }
      }
      miniTab = matrix(0, 2, 2, dimnames = list(c(FALSE, TRUE), c(FALSE, TRUE)))
      crossTab = table(noisyGeno, fullGeno, useNA = "no")
      miniTab[rownames(crossTab), colnames(crossTab)] = crossTab
      pval = phyper(miniTab[1,1] - 1, sum(miniTab[1,]), sum(miniTab[2,]), sum(miniTab[,1]), lower.tail = FALSE, log.p = TRUE)
      finalMat = bind_cols(ID = paste0("P", 1:(numPatients * ifelse(validation, 2L, 1L))), truth = noisyGeno, goodMat, curPheno)
      write_csv(finalMat, paste0("TriageInput_", specString, ".csv"))
      save(allTerms, useSNPs, perturbPositions, extraPerturbation, miniTab, pval, file = paste0("TriageDraws_", specString, ".RData"))
      pvals[index] = pval
    }
  }
  pvals
}
     
## This function bins the genotypes in a specified directory according to their frequency
## To save space, we are only recording the SNP-ID pairs present in the data (if freq < 0.5)
## and those absent in the data (if freq >= 0.5); this may be converted to a SLAM format later
binGenotypesByFrequency = function(directory = "CPLEXInputs", ext = "csv.gz", binWidth = 0.05, pattern = NULL, Subset = NULL) {
  initDir = getwd()
  setwd(directory)
  List = list.files() %>%
    str_subset(ext)
  if (!is.null(pattern)) {
    List %<>% 
      str_subset(pattern)
  }
  numBins = as.integer(1/binWidth) + 1L
  allBins = vector("list", numBins)
  numSNPs = str_extract_all(List, "[0-9]+") %>%
    map_chr(dplyr::last) %>%
    as.integer()
  L = length(List)
  print(paste("There are", L, "files to process"))
  for (index in 1:L) {
    curFname = List[[index]]
    print(paste("Processing", curFname))
    curSize  = numSNPs[index]
    print(paste("There are", curSize, "SNPs to process"))
    curTab = read_csv(curFname, guess_max = Inf, show_col_types = FALSE)
    if (!("ID" %in% colnames(curTab))) {
      curTab %<>% rename(ID = id)
    }
    if (!is.null(Subset)) {
      curTab %<>% filter(ID %in% Subset)
    }
    curGeno = prepareMatrices(curTab, numSNPs = curSize, complement = 0, keepDuplicates = TRUE)$genotypes
    curIDs  = as.integer(rownames(curGeno))
    curSums = colSums(curGeno, na.rm = TRUE)
    curOrder = order(curSums)
    curGeno = curGeno[, curOrder]
    allSNPs = colnames(curGeno) %>%
      str_sub(3, -3) %>%
      as.integer()
    curSums = curSums[curOrder]
    M = nrow(curGeno)
    curBins = floor(curSums / (M * binWidth)) + 1
    goodInds = unique(curBins)
    print(paste("The SNPs fall into", length(goodInds), "unique bins"))
    for (ind in goodInds) {
      curInds = (curBins == ind)
      curMat = curGeno[, curInds, drop = FALSE]
      curSNPs = allSNPs[curInds]
      if (ind >= numBins / 2) {
        curMat = 1 - curMat
      }
      curPos = which(curMat == 1, arr.ind = TRUE)
      curPos[,1] = curIDs [curPos[,1]]
      curPos[,2] = curSNPs[curPos[,2]]
      rownames(curPos) = NULL
      if (!is.null(allSNPs[[ind]])) {
        allBins[[ind]] %<>% rbind(curPos)
      } else {
        allBins[[ind]] =          curPos
      }
    }
  }
  setwd(initDir)
  for (ind in 1:numBins) {
    allSNPs = allBins[[ind]]
    save(allSNPs, file = paste0("BinnedGenotypes", ind, ifelse(!is.null(Subset), paste0("Size", length(Subset)), ""), ".RData"))
  }
  allBins
}

processAllSimFiles = function(scriptName = "fullScript.sh", goodDir = "../../TriageInputs") {
  initDir <- getwd()
  setwd(goodDir)
  goodFiles <- map_chr(list.files(pattern = 'TriageInput*'), ~{normalizePath(.)})
  setwd(initDir)
  Kvals       <- as.integer(str_sub(str_extract(goodFiles, "\\_K[0-9]+\\_"), 3, -2))
  Lvals       <- as.integer(str_sub(str_extract(goodFiles, "\\_L[0-9]+\\_"), 3, -2))
  numPatients <- as.integer(str_sub(str_extract(goodFiles, "\\_M[0-9]+\\_"), 3, -2))
  MY_TYPES = "DNF"
  N <- length(goodFiles)
  output <- NULL
  if (!is.null(scriptName)) {
    Lines <- rep("", N + 1)
    Lines[1]  = "#!/bin/bash"
    pos <- 2
  }
  for (type in MY_TYPES) {
    objective <- MY_OBJECTIVES
    extremeValue <- MAX_P
    myP <- extremeValue %>% as.character
    myF <- ""
    for (ind in 1:N) {
      K = Kvals[ind]
      L = Lvals[ind]
      N = numPatients[ind]
      filePath <- goodFiles[ind]
      if (!is.null(scriptName)) {
        Lines[pos] <- paste("Rscript", "ExtraDriver.R", "-f", filePath, "-t", type, "-k", K, "-l", L, "-o", objective, "-i", pos - 1)
        if (HPC) {
          Lines[pos] <- paste("bsub", "-q", "\"long\"", "-n", NUM_THREADS, "-R", '"broadwell || skylake"', 
                              "-R", paste0("\"span[ptile=", NUM_THREADS, "]\""), 
                              "-M", paste0(4 * TREE_MEMORY, "MB"), "-R", paste0("\"rusage[mem=", 4 * TREE_MEMORY, "MB]\""),
                              "-oo", paste0("/hpc/grid/scratch/chindl/CPLEXRunsSelected", pos, ".log"), 
                              "-eo", paste0("/hpc/grid/scratch/chindl/CPLEXRunsSelected", pos, ".err"),
                              paste0("'", Lines[pos], "'"))
        }
        pos <- pos + 1
        print(pos)
      }
    }
  }
  if (!is.null(scriptName)) {
    fn <- file(scriptName)
    writeLines(Lines, fn)
    close(fn)
  }
  output
}

analyzeSimResults = function(Dir = "SimulationResults/", phenoNames = PHENO_NAMES, pvalue = MAX_P) {
  initDir = getwd()
  setwd(Dir)
  LF = list.files(pattern = "RData") %>%
    str_remove_all(".RData")
  allSettings = str_extract_all(LF, "[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?", simplify = TRUE)
  mode(allSettings) = "double"
  allSettings %<>%
    set_colnames(c("Clauses", "Items", "Patients", "noise", "iteration", "SNPs")) %>%
    as_tibble() %>%
    mutate_at(c("Clauses", "Items", "Patients", "iteration", "SNPs"), as.integer) %>%
    mutate(fname = LF, .before = 1) %>%
    mutate(triageSuccess = NA, otherTriaged = NA, truePheno = NA, foundPheno = NA, logicEquiv = FALSE, functionEquiv = FALSE) %>%
    mutate(truePVal = NA, foundPVal = NA, ratio = NA, within2 = NA, trueMAF = NA, foundMAF = NA, TP = NA, FP = NA, FN = NA, TN = NA)
  for (index in 1:nrow(allSettings)) {
    if (index %% 10 == 0) { print(index) }
    curRow = allSettings %>%
      slice(index)
    curFile = curRow$fname
    altFile = str_replace(curFile, "TriageInput", "../SimulationDraws/TriageDraws") %>%
      str_replace("_OptimizationResults", ".RData")
    load(altFile)
    curTerms = matrix(PHENO_NAMES[allTerms], nrow = curRow$Clauses, ncol = curRow$Items)
    curRow$truePheno = paste0("(", paste(apply(curTerms, 1, function(x) {paste(x, collapse = " AND ")}), collapse = ") OR ("), ")")
    curRow$trueMAF  = sum(miniTab["TRUE", ])/sum(miniTab)
    curTP      = miniTab["TRUE", "TRUE"]
    curPos     = sum(miniTab["TRUE", ])
    curNeg     = sum(miniTab["FALSE", ])
    curPredPos = sum(miniTab[, "TRUE"])
    curRow$truePVal = phyper(curTP - 1, curPos, curNeg, curPredPos, lower.tail = FALSE, log.p = TRUE)
    load(paste0(curFile, ".RData"))
    curSummary = curRes$summary
    curPheno   = curRes$phenotype$truth
    curRow$foundMAF = sum(curPheno)/length(curPheno)
    curRow$foundPheno = curSummary$formula
    curRow$foundPVal = curSummary$LogFisherExactP
    curRow$ratio = exp(abs(curRow$foundPVal - curRow$truePVal))
    curRow$within2 = (curRow$ratio <= 2 && curRow$ratio >= 1/2)
    if (near(pvalue, 1e-3)) {
      extraFile = str_replace(curFile, "_OptimizationResults", 
        paste0("_DNF_K", curRow$Clauses, "_L", curRow$Items, "_P0.001_FInf_agreement_Triage_Results_Summary.csv"))
      extraTab = read_csv(paste0("../SimulationResultsP1e-3/", extraFile), guess_max = Inf, show_col_types = FALSE)
      curRow$triageSuccess = (extraTab %>% slice(1) %>% pull(LogFisherExactP) <= log(pvalue))
      curRow$otherTriaged = sum(extraTab$status[-1] == "optimal")
    } else {
      curRow$triageSuccess = (curRow$foundPVal <= log(pvalue))
      curRow$otherTriaged = 0
    }
    splitFormula = str_split(curRow$foundPheno, " OR ")[[1]]
    fullySplitFormula <- map(splitFormula, function(x) { 
      x %>% str_trim %>% str_remove("^\\(") %>% str_remove("\\)$") %>% str_split(" AND ") %>% extract2(1) })
    foundTerms = lapply(fullySplitFormula, function(x) {match(x, PHENO_NAMES)})
    maxLength  = max(sapply(foundTerms, length))
    foundTerms %<>%
      lapply(fillToLength, maxLength)
    foundTerms = matrix(unlist(foundTerms), ncol = maxLength, byrow = TRUE)
    allTerms %<>%
      putInCanonicalOrder()
    if (!all(is.na(foundTerms))) {
      foundTerms %<>%
        putInCanonicalOrder()
      curRow$logicEquiv = (all(dim(allTerms) == dim(foundTerms)) && all(allTerms == foundTerms))
      curRow$functionEquiv = curRow$logicEquiv
      if (curRow$logicEquiv) {
        curRow$TP = sum(miniTab["TRUE", ])
        curRow$TN = sum(miniTab["FALSE", ])
        curRow$FP = 0
        curRow$FN = 0
      }
    }
    allSettings[index, ] = curRow
  }
  allSettings %<>%
    mutate(inputFile = str_replace(fname, "_OptimizationResults", ".csv"))
  setwd(initDir)
  allSettings %<>%
    postprocessProblematic()
  write_csv(allSettings, paste0("SimulationSummary", pvalue, ".csv"))
  allSettings
}

putInCanonicalOrder = function(inputMatrix) {
  inputMatrix %<>%
    apply(1, function(x) {unique(sort(x))}, simplify = FALSE)
  maxValue = max(unlist(inputMatrix))
  altMatrix = matrix(FALSE, nrow = maxValue, ncol = length(inputMatrix))
  for (ind in 1:length(inputMatrix)) {
    altMatrix[inputMatrix[[ind]], ind] = TRUE
  }
  altMatrix %<>%
    extractMinimalColumns(keepDuplicates = FALSE)
  inputMatrix = apply(altMatrix, 2, which, simplify = FALSE)
  maxLength = max(sapply(inputMatrix, length))
  inputMatrix %<>%
    lapply(fillToLength, maxLength)
  inputMatrix = matrix(unlist(inputMatrix), ncol = maxLength, byrow = TRUE)
  colnames(inputMatrix) = paste0("L", 1:ncol(inputMatrix))
  bestOrder = do.call(order, as.data.frame(inputMatrix))
  inputMatrix = inputMatrix[bestOrder, , drop = FALSE]
  inputMatrix %<>%
    extractUniqueRows(repeats = FALSE) %>%
    extract2(1)
  rownames(inputMatrix) = paste0("K", 1:nrow(inputMatrix))
  inputMatrix
}

fillToLength = function(Vector, Length) {
  output = c(Vector, rep(tail(Vector, 1), Length - length(Vector)))
  output
}

postprocessProblematic = function(allSettings, Dir = "SimulationInputs/") {
  initDir = getwd()
  setwd(Dir)
  for (index in 1:nrow(allSettings)) {
    if (index %% 10 == 0) { print(index) }
    curRow = allSettings %>%
      slice(index)
    if (curRow$triageSuccess && !(curRow$logicEquiv) && !(curRow$foundPheno == "")) {
      curInputFile = curRow$inputFile
      inputTab   = read_csv(curInputFile, col_types = cols(ID = "c", .default = "l"))
      myMatrices = prepareMatrices(inputTab, numSNPs = 50, complement = 0)
      phenotypes = myMatrices$phenotypes
      truePhenoFull  = applyFormula(parseFormula(curRow$truePheno , type = "DNF"), phenotypes, type = "DNF")
      foundPhenoFull = applyFormula(parseFormula(curRow$foundPheno, type = "DNF"), phenotypes, type = "DNF")
      finalTab = matrix(0, 2, 2, dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE")))
      miniTab = table(truePhenoFull, foundPhenoFull)
      finalTab[rownames(miniTab), colnames(miniTab)] = miniTab
      curRow$TP = finalTab["TRUE",  "TRUE"]
      curRow$FP = finalTab["FALSE", "TRUE"]
      curRow$FN = finalTab["TRUE",  "FALSE"]
      curRow$TN = finalTab["FALSE", "FALSE"]
      curRow$functionEquiv = (curRow$FP == 0 && curRow$FN == 0)
      allSettings[index,] = curRow
    }
  }
  setwd(initDir)
  allSettings
}

postprocessSimResults = function(Dir = "../../SimulationResultsNew", single = FALSE) {
  initDir = getwd()
  setwd(Dir)
  LF = list.files()
  fullList = vector("list", length = length(LF)) 
  for (ind in 1:length(LF)) {
    if (ind %% 100 == 0) { print(ind) }
    fname = LF[ind]
    curTab = read_csv(fname, show_col_types = FALSE)
    stopifnot(nrow(curTab) == 1)
    fullList[[ind]] = curTab
  }
  numPatients <- as.integer(str_sub(str_extract(LF, "\\_M[0-9]+\\_"), 3, -2))
  noiseLevels <- as.double (str_sub(str_extract(LF, "\\_eps[0-9\\.]+\\_"), 5, -2))
  redFilename <- paste0(str_extract(LF, "Triage[A-Za-z0-9\\.\\_]+N50"), ".csv")
  Nums = sapply(fullList, ncol)
  inds0 = which(Nums == min(Nums))
  Tab0 = fullList %>% 
    magrittr::extract(inds0) %>%
    bind_rows() %>%
    mutate(fname = redFilename[inds0], Patients = numPatients[inds0], noise = noiseLevels[inds0], .before = 1)
  inds1 = setdiff(1:length(LF), inds0)
  Tab1 = fullList %>% 
    magrittr::extract(inds1) %>%
    bind_rows() %>%
    mutate(fname = redFilename[inds1], Patients = numPatients[inds1], noise = noiseLevels[inds1], .before = 1)
  badCols = setdiff(colnames(Tab0), colnames(Tab1))
  finalCNames = colnames(Tab0)
  finalCNames[match(badCols, finalCNames)] = paste0(finalCNames[match(badCols, finalCNames)], "_Discovery")
  stopifnot(all(finalCNames %in% colnames(Tab1)))
  colnames(Tab0) = finalCNames
  Tab = bind_rows(Tab0, Tab1) %>%
    rename(Clauses = K, Items = L)
  if (single) {
    Tab %<>%
      mutate(K = 1, L = 1)
  }
  write_csv(Tab, paste0("FullSimulationResults", ifelse(single, "Solo", ""), ".csv"))
  setwd(initDir)
  Tab
}