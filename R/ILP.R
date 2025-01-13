#' Produce instructions controlling the execution of a CPLEX optimization
#' @noRd
makeControlObject = function(timeLimit = TIME_LIMIT, probeLevel = PROBE_LEVEL, emphasis = EMPHASIS, varSelect = VAR_SELECT, cutLevel = CUT_LEVEL) {
  Control = list()
  if (!is.na(timeLimit))   { Control$tilim         = timeLimit }
  if (!is.na(probeLevel))  { Control$probe         = probeLevel}
  if (!is.na(emphasis))    { Control$mipemphasis   = emphasis  }
  if (!is.na(varSelect))   { Control$varsel        = varSelect }
  if (!is.na(cutLevel))    { Control$disjcuts = cutLevel; Control$cliques = cutLevel; Control$flowcovers = cutLevel }
  Control
}

#' Create and solve ILPs (integer linear programs) for maximizing an objective function over input phenotypes
#' The first numSNPs columns of the input table are assumed to be the SNPs of interest; the rest are assumed to be phenotypes
#' The value of objective can be "agreement" or "covariance"
#' The value of type can be "CNF" (conjunctive normal form) or "DNF" (disjunctive normal form)
#' The value of K determines the maximum allowed number of combinations in the output
#' The value of L determines the maximum allowed number of elements in a combination
#' If complement =  1 or 2, the genotypes  are augmented with their negation (col name prefix "NOT") before processing starts.
#' If complement = -1 or 2, the phenotypes are augmented with their negation (col name prefix "NOT") before processing starts.
#' baseFilename is used for storing the boundary when p-value constraints are present; the index is used for distinguishing them
#' If strata = TRUE, the Cochran-Mantel-Haenszel statistic is computed to correct for association with the single best phenotype
#' extremeValue is a log p-value upper bound, above which any solution found is declared infeasible; this leads to default results
#' The program assumes that there may be NAs in the genotypes, but NOT in the phenotypes; if that is not the case it will fail!
#' The NAs in each individual genotype are removed, alongside the corresponding phenotype rows, during the preprocessing stage
#' The return value is a list with two components; the first one is a table of statistics which includes the optimal formulas,
#' and the second one is a matrix containing the complex phenotypes corresponding to each of these formulas, which may be empty
#' If outputAssociations = TRUE, then the output additionally contains the matrix of single genotype - phenotype associations
#' @noRd
createAndSolveILPs =  function(inputTab, numSNPs = 1, objective = "agreement", type = "CNF", K = 3, L = 3, complement = 2,
                                baseFilename = NULL, index = 0, outputAssociations = FALSE, strata = TRUE, extremeValue = log(MAX_P)) {
  stopifnot(objective %in% c("agreement","covariance"))
  stopifnot(extremeValue <= 0)
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
  sumNYs  <- colSums(!genotypes, na.rm = TRUE)
  ns      <- sumYs + sumNYs
  prepValues <- prepareCoeffsAndExtraScores(genotypes, objective, sumYs, ns)
  coeffsTrue <- prepValues$coeffsTrue
  coeffsAgree <- prepValues$coeffsAgree
  extraScores <- prepValues$extraScores
  M <- nrow(genotypes)
  N <- ncol(genotypes)
  objVectors  <- matrix(coeffsAgree, M, N, byrow = TRUE) * genotypes + matrix(coeffsTrue, M, N, byrow = TRUE)
  print("Reducing the input matrices")
  out <- reducePhenotypesAndObjectives(phenotypes, objVectors, genotypes = genotypes)
  phenotypes <- out$phenotypes
  objVectors <- out$objVectors
  extraScores %<>% magrittr::subtract(out$extraScores) ### Replaced addition by subtraction, 11/11/2024
  numSegments <- rep(NA, ncol(genotypes))
  if (!is.null(extremeValue)) { extraDef <- out$extraDef }
  print("Computing single associations")
  singleRes <- computeBestPValues(origGenotypes, origPhenotypes, fisher = TRUE)
  extremeValues <- prepareExtremeValues(extremeValue, type, singleRes$value)
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
    if (!is.null(extremeValue)) {
      extFilename = paste0(baseFilename, "T", type)
      curExt <- prepareExtras(extraDef, ind, goodInds, extremeValue, ns, sumYs, extFilename, index = index)
      numSegments[ind] <- ifelse(any(is.na(curExt$B)), 0, nrow(curExt$B))
    }
    fn <- paste0(ifelse(is.null(baseFilename), "Test", baseFilename), "C", index, "I", ind, "T", type, ".lp")
    curSol <- createAndSolveILP(phenotypes = P, objVector = curObjVector, type = type, K = K, L = L, filename = fn,
                                boundValue = boundValues[ind], extraConstraints = curExt)
    curPheno <- postprocessSolution(origPhenotypes, curSol$usedVars, type = type, complement = complement)
    curPhenoMini <- as.vector(curPheno$solution)
    curGenoMini  <- origGenotypes[, ind, drop = TRUE]
    myTab <- table(curGenoMini, curPhenoMini)
    if (strata) {
      bestSingle <- origPhenotypes[, singleRes$phenotype[[ind]], drop = TRUE]
      if (sum(bestSingle) > 1 && sum(!bestSingle) > 1) { # if the sum is 0, 1, N - 1 or N, a strata is too small and the test fails!
        MH[ind] <- computeStratifiedPVal(genotype = curGenoMini, phenotypeC = curPhenoMini, phenotypeS = bestSingle)
      }
    }
    checkResult <- computeOptValue(myTab, objective, coeffsTrue, coeffsAgree, ind)
    checkOpt    <- checkResult$checkOpt
    if (is.na(checkOpt)) { checkOpt <- 0 }
    if (is.na(curSol$optimum)) {
      curOpt <- checkOpt
    } else {
      curOpt <- curSol$optimum + extraScores[ind]
      stopifnot(dplyr::near(checkOpt, curOpt, tol = .Machine$double.eps^0.5 * 2 * max(ns[ind], nrow(origGenotypes) - ns[ind])))
    }
    curOutput <- list(time = curSol$time, gap = curSol$gap, opt = curOpt, status = curSol$status, formula = curPheno$clause) %>%
        c(TP = checkResult$TP, FP = checkResult$FP, FN = checkResult$FN, TN = checkResult$TN)
    output[[curName]] <- curOutput
    complexPs[, curName] <- curPhenoMini
  }
  fullOutput <- prepareOutput(output, complexPhenotypes = complexPs, SNP = colnames(genotypes), ID = rownames(origPhenotypes),
                              type = type, numSegments = numSegments, objective = objective, singleRes = singleRes,
                              boundValues = boundValues, ns = ns, MH = MH, outputAssociations = outputAssociations)
  fullOutput
}
