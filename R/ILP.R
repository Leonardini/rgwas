srcFiles <- c("Geometry", "Logical", "Preprocess", "Postprocess", "Utilities", "Validation")
for (srcFile in srcFiles) { source(paste0(srcFile, ".R")) }

### This function produces instructions controlling the execution of a CPLEX optimization to write to a script file
makeControlLines = function(parMode = PAR_MODE, nThreads = NUM_THREADS, treeMemory = TREE_MEMORY, timeLimit = TIME_LIMIT, 
                            probeLevel = PROBE_LEVEL, emphasis = EMPHASIS, varSelect = VAR_SELECT, cutLevel = CUT_LEVEL, 
                            cutPasses = CUT_PASSES, repairTries = REPAIR_TRIES, heuristicFr = HEURISTIC_FR, nSolutions = NUM_SOL,
                            quadCPTol = QUAD_CP_TOL, quadCPStrat = QUAD_CP_STRAT, writeLvl = WRITE_LEVEL, boundValue = NA,
                            numEmphasis = NUM_EMPHASIS, intTolerance = INTTOLERANCE, API = FALSE) {
  if (API) {
    Control = list()
    if (!is.na(parMode))     { Control$parallel.mode = parMode   }
    if (!is.na(nThreads))    { Control$threads       = nThreads  }
    if (!is.na(timeLimit))   { Control$tilim         = timeLimit }
    if (!is.na(probeLevel))  { Control$probe         = probeLevel}
    if (!is.na(emphasis))    { Control$mipemphasis   = emphasis  }
    if (!is.na(varSelect))   { Control$varsel        = varSelect }
    if (!is.na(cutLevel))    { Control$disjcuts = cutLevel; Control$cliques = cutLevel; Control$flowcovers = cutLevel }
    if (!is.na(treeMemory))  { warnOption("tree memory")        }
    if (!is.na(cutPasses))   { warnOption("cut passes")         }
    if (!is.na(repairTries)) { warnOption("repair tries")       }
    if (!is.na(heuristicFr)) { warnOption("heuristic frequency")}
    if (!is.na(writeLvl))    { warnOption("write level")        }
    if (!is.na(boundValue))  { warnOption("upper cutoff")       }
    if (!is.na(numEmphasis)) { warnOption("numerical emphasis") }
    if (!is.na(intTolerance)){ warnOption("integral tolerance") }
    Control
  } else {
    Lines <- ""
    if (!is.na(parMode))     { Lines <- c(Lines, paste(c("set", "parallel", parMode), collapse = " "))                            }
    if (!is.na(nThreads))    { Lines <- c(Lines, paste(c("set", "threads", nThreads), collapse = " "))                            }
    if (!is.na(treeMemory))  { Lines <- c(Lines, paste(c("set", "mip", "limits", "treememory", treeMemory), collapse = " "))      }
    if (!is.na(timeLimit))   { Lines <- c(Lines, paste(c("set", "timelimit", timeLimit), collapse = " "))                         }
    if (!is.na(probeLevel))  { Lines <- c(Lines, paste(c("set", "mip", "strategy", "probe", probeLevel), collapse = " "))         }
    if (!is.na(emphasis))    { Lines <- c(Lines, paste(c("set", "emphasis", "mip", emphasis), collapse = " "))                    }
    if (!is.na(varSelect))   { Lines <- c(Lines, paste(c("set", "mip", "strategy", "variableselect", varSelect), collapse = " ")) }
    if (!is.na(cutLevel))    { Lines <- c(Lines, paste(c("set", "mip", "cuts", "all", cutLevel), collapse = " "))                 }
    if (!is.na(cutPasses))   { Lines <- c(Lines, paste(c("set", "mip", "limits", "cutpasses", cutPasses), collapse = " "))        }
    if (!is.na(repairTries)) { Lines <- c(Lines, paste(c("set", "mip", "limits", "repairtries", repairTries), collapse = " "))    }
    if (!is.na(heuristicFr)) { Lines <- c(Lines, paste(c("set", "mip", "strategy", "heuristicfreq", heuristicFr), collapse = " "))}
    ## if (!is.na(nSolutions))  { Lines <- c(Lines, paste(c("set", "mip", "limits", "solutions", nSolutions), collapse = " "))       }
    ## if (!is.na(quadCPTol))   { Lines <- c(Lines, paste(c("set", "barrier", "qcpconvergetol", quadCPTol), collapse = " "))         }
    ## if (!is.na(quadCPStrat)) { Lines <- c(Lines, paste(c("set", "mip", "strategy", "miqcpstrat", quadCPStrat), collapse = " "))   }
    if (!is.na(writeLvl))    { Lines <- c(Lines, paste(c("set", "output", "writelevel", writeLvl), collapse = " "))               }
    if (!is.na(boundValue))  { Lines <- c(Lines, paste(c("set", "mip", "tolerances", "uppercutoff", boundValue), collapse = " ")) }
    if (!is.na(numEmphasis)) { Lines <- c(Lines, paste(c("set", "emphasis", "numerical", "yes"), collapse = " "))                 }
    if (!is.na(intTolerance)){ Lines <- c(Lines, paste(c("set", "mip", "tolerances", "integrality", intTolerance),collapse = " "))}
    Lines
  }
}

warnOption = function(optName) {
  message = paste("Warning: the", optName, "option has not yet been implemented in the API; try setting it via a parameter file")
  message
}

### This function is a simplified version of that in Alternative.R in the ClassTR package published in PLoS Comp Bio
createScript = function(scriptName = "auxScript", filename = NULL, startSol = NULL, boundValue = NA, save = FALSE) {
  topLines  <- c("#!/bin/bash", "source ~/.bash_profile")
  initLines <- c(paste0("for file in ", filename, "; do"), "cat <<EOF | ./cplex | tee $file.log", "set log *")
  ctrlLines <- makeControlLines(boundValue = boundValue)
  botLines  <- rep("", 10)
  botLines[1] <- "read $file"
  if (!is.null(startSol)) { 
    botLines[2] <- paste(c("read", startSol), collapse = " ")
  }
  botLines[3:4] <- c("opt", "disp sol var *")
  if (save) {
    botLines[5] <- paste(c("write", paste0(filename, ".sol")), collapse = " ")
  }
  botLines[6:7] = c("EOF", "done")
  Lines <- c(topLines, initLines, ctrlLines, botLines)
  Lines <- Lines[Lines != ""]
  fw = file(scriptName, "w")
  writeLines(Lines, fw)
  close(fw)
  system2("chmod", args = c("+x", scriptName))
}

### This function is a simplified version of that in Alternative.R in the ClassTR package published in PLoS Comp Bio
runFile = function(inputFile, boundValue = NA, removeScript = TRUE, removeInput = FALSE, stdout = "", startSol = NULL) {
  scriptName <- paste0("auxScript", str_remove(tail(str_split(inputFile, "\\/")[[1]], 1), "\\.[A-Za-z]{2,3}"))
  createScript(scriptName = scriptName, filename = inputFile, boundValue = boundValue, startSol = startSol)
  system2("/bin/bash", paste0("./", scriptName), stdout = stdout)
  if (removeScript) {
    file.remove(scriptName)
  }
  if (removeInput) {
    file.remove(inputFile)
  }
}

### This function is a simplified version of that in Alternative.R in the ClassTR package published in PLoS Comp Bio
runAndParse = function(inputFile, boundValue = NA, cplexDir = CPLEX_DIR, prefix1 = PREFIX1, prefix2 = PREFIX2, 
                       removeInput = TRUE, removeLog = TRUE, startSol = NULL, optFactor = -1L) {
  inputFile %<>%
    normalizePath
  if (!is.null(startSol)) {
    startSol %<>%
      normalizePath
  }
  initDir <- getwd()
  setwd(CPLEX_DIR)
  runFile(inputFile, boundValue = boundValue, startSol = startSol)
  LOGfile <- str_c(inputFile, ".log")
  resultLines <- readLines(LOGfile, warn = FALSE)
  if (removeInput) {
    file.remove(inputFile)
  }
  if (removeLog) {
    file.remove(LOGfile)
  }
  setwd(initDir)
  output <- parseLines(resultLines, prefix1, prefix2, optFactor = optFactor)
  output
}

### This function creates and solves an ILP using a call to CPLEX and parses the solution into an appropriate format
createAndSolveProblem = function(coeffMatrix, objVector, type, K, L, filename, boundValue = NA, startSol = NULL, extraConstraints = NULL) {
  if (type == "CNF" || type == "DNF") {
    output <- createAndSolveILP(coeffMatrix, objVector, type = type, K = K, L = L, filename = filename, boundValue = boundValue,
                                startSol = startSol, extraConstraints = extraConstraints)
  } else { ### type = "LCP"
    output <- createAndSolveMIQCP(coeffMatrix, objVector, K = K, M = L, filename = filename, startSol = startSol, 
                                  boundValue = boundValue)
  }
  output
}

### This function creates and solves ILPs (integer linear programs) for maximizing an objective function over input phenotypes
### The first numSNPs columns of the input table are assumed to be the SNPs of interest; the rest are assumed to be phenotypes
### The value of objective can be "agreement", "covariance", "correlation", "OR" (odds ratio) or "Cohen" (Cohen's kappa)
### NOTE: only the first two implemented at the moment because all the others involve ratios and are thus non-linear functions!
### The value of type can be "CNF" (conjunctive normal form), "DNF" (disjunctive normal form), or "LCP" (a linear combination)
### The value of K determines the maximum allowed number of combinations (for CNF or DNF) or phenotypes (for LCP) in the output
### The value of L determines the maximum allowed number of elements in a combination (this is used for CNF or DNF, not for LCP)
### If complement =  1 or 2, the genotypes  are augmented with their negation (col name prefix "NOT") before processing starts.
### If complement = -1 or 2, the phenotypes are augmented with their negation (col name prefix "NOT") before processing starts.
### baseFilename is used for storing the boundary when p-value constraints are present; the index is used for distinguishing them
### If strata = TRUE, the Cochran-Mantel-Haenszel statistic is computed to correct for association with the single best phenotype
### extremeValue is a log p-value upper bound, above which any solution found is declared infeasible; this leads to default results
### The program assumes that there may be NAs in the genotypes, but NOT in the phenotypes; if that is not the case it will fail!
### The NAs in each individual genotype are removed, alongside the corresponding phenotype rows, during the preprocessing stage
### The return value is a list with two components; the first one is a table of statistics which includes the optimal formulas,
### and the second one is a matrix containing the complex phenotypes corresponding to each of these formulas, which may be empty
### If outputAssociations = TRUE, then the output additionally contains the matrix of single genotype - phenotype associations
createAndSolveILPs =  function(inputTab, numSNPs = 1, objective = "agreement", type = "CNF", K = 3, L = 3, complement = 2,
                                baseFilename = NULL, index = 0, outputAssociations = FALSE, strata = TRUE,
                                extremeValue = ifelse(type == "LCP", 0, log(MAX_P))) {
  stopifnot(type %in% c("CNF","DNF") && objective %in% c("agreement","covariance") || type == "LCP" && objective == "correlation")
  if (type == "LCP" && complement != 0) { warning("Setting complement to ", complement, " with type LCP will slow you down!") }
  stopifnot(type == "LCP" || extremeValue <= 0)
  myMatrices <- prepareMatrices(inputTab, numSNPs, complement)
  phenotypes <- myMatrices$phenotypes
  genotypes  <- myMatrices$genotypes
  print(paste("There are", ncol(phenotypes), "phenotypes over", nrow(phenotypes), "patients"))
  print(paste("There are", ncol(genotypes), "genotypes over", nrow(genotypes), "patients"))
  origPhenotypes <- phenotypes
  origGenotypes  <- genotypes
  output <- vector("list", ncol(genotypes))
  names(output) <- colnames(genotypes)
  complexPs <- matrix(NA, nrow(phenotypes), ncol(genotypes), dimnames = list(rownames(phenotypes), names(output)))
  sumYs  <- colSums(genotypes, na.rm = TRUE)
  if (type == "LCP") {
    ns      <- colSums(!is.na(genotypes))
    sumYSqs <- colSums(genotypes^2, na.rm = TRUE)
  } else {
    sumNYs  <- colSums(!genotypes, na.rm = TRUE)
    ns      <- sumYs + sumNYs
  }
  prepValues <- prepareCoeffsAndExtraScores(genotypes, objective, sumYs, ns)
  coeffsTrue <- prepValues$coeffsTrue
  coeffsAgree <- prepValues$coeffsAgree
  extraScores <- prepValues$extraScores
  M <- nrow(genotypes)
  N <- ncol(genotypes)
  objVectors  <- matrix(coeffsAgree, M, N, byrow = TRUE) * genotypes + matrix(coeffsTrue, M, N, byrow = TRUE)
  print("Reducing the input matrices")
  if (type != "LCP") {
    out <- reducePhenotypesAndObjectives(phenotypes, objVectors, columnsOnly = FALSE, genotypes = genotypes)
  } else {
    out <- reducePhenotypesAndObjectives(phenotypes, objVectors, columnsOnly = TRUE , genotypes = NULL)
  }
  phenotypes <- out$phenotypes
  objVectors <- out$objVectors
  extraScores %<>% add(out$extraScores)
  numSegments <- rep(NA, ncol(genotypes))
  if (type != "LCP" && !is.null(extremeValue)) { extraDef <- out$extraDef }
  print("Computing single associations")
  singleRes <- computeBestPValues(origGenotypes, origPhenotypes, cor = (type == "LCP"), fisher = (type != "LCP"))
  extremeValues <- prepareExtremeValues(extremeValue, type, singleRes$value, SINGLE_BEST_RATIO, USE_TIGHTER_BOUND)
  boundValues   <- prepareBoundValues(extremeValues, type, objective, ns, sumYs, coeffsTrue, coeffsAgree)
  MH <- rep(NA, ncol(genotypes))
  print(paste("There are", ncol(genotypes), "SNPs to process"))
  for (ind in 1:ncol(genotypes)) {
    extremeValue <- extremeValues[ind]
    print(paste("Currently processing SNP", colnames(genotypes)[ind]))
    curName <- colnames(genotypes)[ind]
    curObjVector <- objVectors[, ind, drop = FALSE]
    if (any(is.na(curObjVector))) { print(paste("Removing", sum(is.na(curObjVector)), "missing entries in", curName)) }
    goodInds <- which(!(is.na(curObjVector)))
    P <- phenotypes[goodInds, , drop = FALSE]
    curObjVector <- curObjVector[goodInds, , drop = FALSE]
    curExt <- NULL
    if (type != "LCP" && !is.null(extremeValue)) {
      extFilename = paste0(baseFilename, "T", type)
      curExt <- prepareExtras(extraDef, ind, goodInds, extremeValue, ns, sumYs, extFilename, index = index)
      numSegments[ind] <- ifelse(any(is.na(curExt$B)), 0, nrow(curExt$B))
    }
    fn <- paste0(ifelse(is.null(baseFilename), "Test", baseFilename), "C", index, "I", ind, "T", type, ".lp")
    if (type == "LCP") {
      G <- genotypes[goodInds, ind, drop = FALSE]
      initBound <- boundValues[ind]
      S <- sqrt(sumYSqs[ind] - sumYs[ind]^2 / ns[ind])
      curResult <- processLCP(G, P, S, ind, initBound, startSol = NULL, extremeValue, K = K, L = L, filename = fn)
      curSol    <- curResult$solution
      upperBound <- curResult$bound
      boundValues[ind] <- upperBound
      if (complement > 0 && ind < length(ns)/2) {
        boundValues[ind + length(ns)/2] <- upperBound
      }
    } else {
      curSol <- createAndSolveILP(phenotypes = P, objVector = curObjVector, type = type, K = K, L = L, filename = fn, 
                                boundValue = boundValues[ind], startSol = NULL, extraConstraints = curExt, API = TRUE)
    }
    curPheno <- postprocessSolution(origPhenotypes, curSol$usedVars, type = type, complement = complement)
    curPhenoMini <- curPheno$solution
    curGenoMini  <- origGenotypes[, ind, drop = TRUE]
    myTab <- table(curGenoMini, curPhenoMini)
    if (strata) {
      bestSingle <- origPhenotypes[, singleRes$phenotype[[ind]], drop = TRUE]
      if (sum(bestSingle) > 1 && sum(!bestSingle) > 1) { ### if the sum is 0, 1, N - 1 or N, a strata is too small and the test fails! 
        MH[ind] <- computeStratifiedPVal(genotype = curGenoMini, phenotypeC = curPhenoMini, phenotypeS = bestSingle)
      }
    }
    if (type != "LCP") {
      checkResult <- computeOptValue(myTab, objective, coeffsTrue, coeffsAgree, ind)
      checkOpt    <- checkResult$checkOpt
    } else {
      checkOpt    <- cor(curGenoMini, curPhenoMini[,1], use = "pairwise")
      if (!is.na(checkOpt) && checkOpt < 0) { ### implicitly multiply the coefficients by -1
        checkOpt     <- -checkOpt
        curPhenoMini <- -curPhenoMini
        curPheno$clause <- negateCoefficients(curPheno$clause)
      }
    }
    if (is.na(checkOpt)) { checkOpt <- 0 }
    if (is.na(curSol$optimum)) {
      curOpt <- ifelse(type == "LCP", NA, checkOpt)
    } else {
      curOpt <- curSol$optimum + extraScores[ind]
      stopifnot(near(checkOpt, curOpt, tol = .Machine$double.eps^0.5 * 2 * max(ns[ind], nrow(origGenotypes) - ns[ind])))
    }
    curOutput <- list(time = curSol$time, gap = curSol$gap, opt = curOpt, status = curSol$status, formula = curPheno$clause)
    if (type != "LCP") {
      curOutput %<>% 
        c(TP = checkResult$TP, FP = checkResult$FP, FN = checkResult$FN, TN = checkResult$TN)
    }
    output[[curName]] <- curOutput
    complexPs[, curName] <- curPhenoMini
  }
  fullOutput <- prepareOutput(output, complexPhenotypes = complexPs, SNP = colnames(genotypes), ID = rownames(origPhenotypes), 
                              type = type, numSegments = numSegments, objective = objective, singleRes = singleRes, 
                              boundValues = boundValues, ns = ns, MH = MH, outputAssociations = outputAssociations)
  fullOutput
}
