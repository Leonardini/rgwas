mainDriver = function(inputFile, type = "DNF", objective = "agreement", K = MY_K, L = MY_L, complement = FALSE, 
                      ext = str_sub(inputFile, -4), index = 1L, extremeValue = ifelse(type == "LCP", 0, log(MAX_P)), 
                      outputSummary = TRUE, outputPhenotype = FALSE, outputAssociations = FALSE, saveOutput = FALSE) {
  miniFile <- inputFile
  if (str_sub(inputFile, -3) == ".gz") {
    ext <- str_sub(inputFile, -7, -4)
    miniFile <- str_sub(inputFile, 1, -4)
  }
  readFunction <- ifelse(ext == '.tsv', read_tsv, read_csv)
  inputTab <- readFunction(inputFile, col_types = cols(ID = "c", id = "c", .default = ifelse(type == "LCP", "d", "l")))
  if ("id" %in% colnames(inputTab)) { inputTab <- inputTab %>% rename(ID = id) }
  numSNPs <- str_extract(miniFile, paste0("([0-9]+)", ext)) %>% str_remove(ext) %>% parse_integer
  numPheno <- (ncol(inputTab) - 1 - numSNPs) * (1 + (complement %in% c(-1, 2)))
  numCombos <- choose(choose(numPheno, L), K)
  myP <- extremeValue %>% as.character
  myF <- ifelse(!is.na(SINGLE_BEST_RATIO), paste0("F", SINGLE_BEST_RATIO), "")
  baseFilename <- str_remove(miniFile, ext)
  stage <- ifelse(OPTIMIZATION, "Optimization", "Triage")
  specString <- paste(c("", type, paste0("K", K), paste0("L", L), paste0(ifelse(type == "LCP", "R", "P"), myP), myF, objective, 
                        stage, "Results"), collapse = "_")
  outFileSum <- str_replace(miniFile, ext, str_c(specString, '_Summary.csv'))
  if (!(file.exists(outFileSum))) {
    output <- createAndSolveILPs(inputTab, numSNPs = numSNPs, type = type, objective = objective, K = K, L = L, index = index, 
                                 extremeValue = extremeValue, complement = complement, baseFilename = baseFilename, outputAssociations = FALSE, strata = TRUE)
    if (outputSummary && !is.null(output$summary)) {     
      write_csv(output$summary,      str_replace(miniFile, ext, str_c(specString,      '_Summary.csv')))
    }
    if (outputPhenotype && !is.null(output$phenotype))   { 
      write_csv(output$phenotype,    str_replace(miniFile, ext, str_c(specString,    '_Phenotype.csv')))
    }
    if (outputAssociations && !is.null(output$associations)) {
      write_csv(output$associations, str_replace(miniFile, ext, str_c(specString, '_Associations.csv')))
    }
    if (saveOutput)    {    
      saveRDS(object = output, file = str_replace(miniFile, ext, str_c(specString,              '.RDS')))
    }
    output
  }
}

mainDriverExtended = function(inputFile, type, objective, complement, KLpairs) {
  Dim <- ncol(KLpairs)
  for (index in 1:Dim) {
    curPair <- KLpairs[,index]
    curK <- curPair[1]
    curL <- curPair[2]
    mainDriver(inputFile, type = type, objective = objective, K = curK, L = curL, complement = complement, index = index)
  }
}

preOptimizingDriver = function(fname, type = "DNF", complement = FALSE, K = 3, L = 3, index = 0, pvalue = log(MAX_P)) {
  miniFile = str_replace(fname, "N50.csv", "truth_N1.csv")
  Tab = read_csv(fname, col_types = cols(ID = "c", id = "c", .default = "l"))
  miniTab = Tab %>% select(!starts_with("S"))
  write_csv(miniTab, miniFile)
  curRes = optimizingDriver(inputFile = miniFile, startPValue = pvalue, factor = 2, type = type, complement = complement, 
                            K = K, L = L, startIndex = index)
  save(curRes, file = str_replace(fname, ".csv", "_OptimizationResults.RData"))
  curRes
}
 
optimizingDriver = function(inputFile, startPValue = log(1e-10), factor = 2, type = "CNF", complement = -1, K = MY_K, L = MY_L, startIndex = 0) {
  stopifnot(startPValue <= 0)
  curPValue <- startPValue
  ind <- 0
  curOutput <- mainDriver(inputFile, extremeValue = curPValue, index = startIndex * THOUSAND + ind, complement = complement, 
                  type = type, K = K, L = L, outputSummary = FALSE, outputPhenotype = FALSE, outputAssociations = FALSE)
  prevOutput <- curOutput
  while (curOutput$summary$status[1] %in% c("feasible", "optimal")) {
    prevOutput <- curOutput
    newPValue <- prevOutput$summary$LogFisherExactP[1]
    if (is.na(newPValue)) {
      print("Warning: The p-value was NA, so using the previous value times factor as a stopgap measure!")
      newPValue <- curPValue
    }
    curPValue <- newPValue - log(factor)
    ind <- ind + 1
    curOutput <- mainDriver(inputFile, extremeValue = curPValue, index = startIndex * THOUSAND + ind, complement = complement, 
                    type = type, K = K, L = L, outputSummary = FALSE, outputPhenotype = FALSE, outputAssociations = FALSE)
    if (curOutput$summary$status[1] %in% c("timeout", "memoryout")) {
      print(paste("Warning: the optimization is terminated prematurely due to insufficient CPLEX resources (time or memory)"))
      break
    }
  }
  print(paste("The best feasible p-value to within a factor of", factor, "is", exp(prevOutput$summary$LogFisherExactP[1])))
  return(prevOutput)
}

processAllFiles = function(scriptName = "fullScript.sh", goodDir = paste0('/Users/chindl/Documents/GitHub/ReverseGWAS/AsaData/InputFiles/', MY_TYPES[1], 'Inputs/')) {
  initDir <- getwd()
  setwd(goodDir)
  goodFiles <- map_chr(list.files(pattern = paste0('CPLEXInput.*', EXT)), ~{normalizePath(.)})
  setwd(initDir)
  if (length(goodFiles) == 0) { goodFiles = "CPLEXInputs/CPLEXInput_SelfReported_EU_100_15.csv" }
  myDatabase  <- str_split(str_remove(goodFiles, EXT), "\\_")
  inputSize   <- parse_integer(sapply(myDatabase, function(x) {x[length(x)]}))
  numPatients <- parse_integer(sapply(myDatabase, function(x) {x[length(x) - 1]}))
  bestOrder <- order(numPatients, inputSize)
  goodFiles <- goodFiles[bestOrder]
  n1 <- length(MY_TYPES)
  n2 <- length(MY_K)
  n3 <- length(MY_L)
  stopifnot(n2 == n3)
  n4 <- length(MY_OBJECTIVES)
  n5 <- length(goodFiles)
  N <- n1 * n2 * n4 * n5
  if (!is.null(scriptName)) {
    output <- NULL
    Lines <- rep("", N + 1)
    Lines[1]  = "#!/bin/bash"
    pos <- 2
  } else {
    output <- vector("list", N)
    Names <- outer(outer(outer(outer(goodFiles, MY_TYPES, function(x,y) {paste0(x, "Type", y)}), MY_OBJECTIVES, function(x,y) {paste0(x, "Objective", y)}), 
                         MY_K, function(x,y) {paste0(x, "K", y)}), MY_L, function(x,y) {paste0(x, "L", y)})
    names(output) <- Names
  }
  for (type in MY_TYPES) {
    extremeValue <- ifelse(type == "LCP", 0, MAX_P)
    myP <- extremeValue %>% as.character
    myF <- ifelse(!is.na(SINGLE_BEST_RATIO), paste0("F", SINGLE_BEST_RATIO), "")
    for (ind in 1:n2) {
      K = MY_K[ind]
      L = MY_L[ind]
      for (objective in MY_OBJECTIVES) {
        for (curFile in goodFiles) {
          # print(paste("Currently processing", curFile))
          stage <- ifelse(OPTIMIZATION, "Optimization", "Triage")
          fullFile <- curFile
          curExt <- paste0("\\", str_sub(curFile, -nchar(EXT)))
          specString <- paste(c("", type, paste0("K", K), paste0("L", L), paste0(ifelse(type == "LCP", "R", "P"), myP), myF,
                                objective, stage, "Results"), collapse = "_")
          outFile  <- str_replace(fullFile, curExt, str_c(specString, '_Summary.csv'))
          if (!(file.exists(outFile))) {
            print(outFile)
            if (!is.null(scriptName)) {
              filePath <- fullFile
              Lines[pos] <- paste("Rscript", "ExtraDriver.R", "-f", filePath, "-t", type, "-k", K, "-l", L, "-o", 
                                  objective, "-i", pos - 1)
              if (HPC) {
                Lines[pos] <- paste("bsub", "-q", "\"at_risk\"", "-n", NUM_THREADS, "-R", '"broadwell || skylake"', 
                                    "-R", paste0("\"span[ptile=", NUM_THREADS, "]\""), 
                                    "-M", paste0(4 * TREE_MEMORY, "MB"), "-R", paste0("\"rusage[mem=", 4 * TREE_MEMORY, "MB]\""),
                                    "-oo", paste0("/hpc/grid/scratch/chindl/CPLEXRunsSelected", pos, ".log"), 
                                    "-eo", paste0("/hpc/grid/scratch/chindl/CPLEXRunsSelected", pos, ".err"),
                                    paste0("'", Lines[pos], "'"))
              } else {
                Lines[pos] <- paste(Lines[pos], rep(" &", (pos %% 2 == 0)))
              }
              pos <- pos + 1
            } else {
              print(paste("Processing", fullFile, "with objective", objective, "and type", type))
              curResult <- mainDriver(inputFile = fullFile, type = type, objective = objective, K = K, L = L, complement = 2, 
                                      extremeValue = extremeValue)
              output[[paste0(fullFile, "Type", type, "Objective", objective, "K", K, "L", L)]] <- curResult
            }
          }
        }
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
