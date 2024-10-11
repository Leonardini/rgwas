#' Compute the matrix of pairwise covariances between columns of matrix1 and matrix2, assumed to have the same number of rows
#' If scaleUp = TRUE, do not normalize by the number of non-missing pairs.
computeCovariances = function(matrix1, matrix2, scaleUp = TRUE) {
  stopifnot(nrow(matrix1) == nrow(matrix2))
  scaledCovariances <- cov(matrix1, matrix2, use = "pairwise.complete.obs")
  if (scaleUp) {
    numNonMissingPairs <- crossprod(!is.na(matrix1), !is.na(matrix2))
    scaledCovariances %<>%
      magrittr::multiply_by(numNonMissingPairs - 1)
  }
  scaledCovariances
}

#' Compute the column with best metric in matrix2 for each column of matrix1, returning both its value (arg1) and name (arg2).
#' The matrices are assumed (without checking) to have the same number of rows and to have binary values (0/1).
#' If fisher = TRUE use the Fisher exact p-value, otherwise chi-squared log p-value.
computeBestPValues = function(matrix1, matrix2, fisher = TRUE) {
  M <- nrow(matrix1)
  stopifnot(nrow(matrix2) == M)
  N1 <- ncol(matrix1)
  bestPs <- rep(NA, N1)
  bestps <- rep("", N1)
  N2 <- ncol(matrix2)
  rown <- colnames(matrix1)
  coln <- colnames(matrix2)
  print(paste("There are", N1, "columns to process"))
  pFun <- function(x) {
    ifelse(fisher, phyper(x[1] - 1, x[1] + x[2], x[3] + x[4], x[1] + x[3], lower.tail = FALSE, log.p = TRUE),
           pchisq((sum(x) * (x[1] * sum(x) - (x[1] + x[2]) * (x[1] + x[3]))^2) /
                    ((x[1] + x[2]) * (x[3] + x[4]) * (x[1] + x[3]) * (x[2] + x[4])), df = 1, lower.tail = FALSE, log.p = TRUE)) }
  allValues <- matrix(NA, N1, N2)
  for (ind1 in 1:N1) {
    if (ind1 %% 100 == 0) { print(ind1) }
    curCol <- matrix1[, ind1]
    curInds <- which(!is.na(curCol))
    curM <- length(curInds)
    curSum <- sum(curCol, na.rm = TRUE)
    curTPs <- colSums(matrix2 & curCol, na.rm = TRUE)
    curFPs <- curSum - curTPs
    cSums  <- colSums(matrix2[curInds, , drop = FALSE])
    curFNs <- cSums - curTPs
    curTNs <- curM - (curFPs + curFNs + curTPs)
    curMat <- cbind(curTPs, curFPs, curFNs, curTNs)
    curPVals <- apply(curMat, 1, pFun)
    allValues[ind1, ] <- curPVals
    bestValue    <- min(curPVals, na.rm = TRUE)
    bestPs[ind1] <- bestValue
    bestps[ind1] <- coln[which(dplyr::near(curPVals, bestValue, tol = .Machine$double.eps^0.5 * M))[1]]
  }
  dimnames(allValues) <- list(rown, coln)
  output <- list(value = bestPs, phenotype = bestps, allValues = allValues)
  output
}

#' Compute the matrix of pairwise agreements between columns of matrix1 and matrix2, assumed to have the same number of rows
computeAgreements = function(matrix1, matrix2) {
  stopifnot(nrow(matrix1) == nrow(matrix2))
  N1 <- ncol(matrix1)
  N2 <- ncol(matrix2)
  if (N2 < N1) {
    tempMat <- matrix1
    matrix1 <- matrix2
    matrix2 <- tempMat
    Nt <- N1
    N1 <- N2
    N2 <- Nt
  }
  output <- matrix(NA, N1, N2, dimnames = list(colnames(matrix1), colnames(matrix2)))
  for (ind in 1:N1) {
    curMat <- matrix(matrix1[, ind], nrow(matrix2), N2)
    agree <- (curMat == matrix2)
    disagree <- (curMat != matrix2)
    output[ind,] <- colSums(agree, na.rm = TRUE) - colSums(disagree, na.rm = TRUE)
  }
  output
}

#' This function computes, for a given N and nPlus, the set of critical value pairs (i, j) such that the 2x2 contingency table
#' with row sums nPlus and N - nPlus, and column sums i and N - i has log p-value below pMax when its top left entry is at least j.
computePValueBounds = function(N, nPlus, pMax) {
  print(paste("The log p-value being considered as the cutoff is", pMax))
  print(paste("Preparing p-value bounds with N =", N, "and nPlus =", nPlus))
  stopifnot(nPlus <= N)
  nMinus <- N - nPlus
  if (nPlus == 0) {
    return(tibble::tibble(Sum = integer(), Total = integer()))
  }
  iRange <- 1:N
  criticalValues <- rep(NA, N)
  prevValue <- 0
  for (i in iRange) {
    range <- max(prevValue - M_STEP, 0, na.rm = TRUE) : min(prevValue + M_STEP, i, na.rm = TRUE)
    if (nPlus > nMinus) {
      hyperVector <- phyper(nMinus - i + range, nMinus, nPlus , N - i, lower.tail = FALSE, log.p = TRUE)
    } else {
      hyperVector <- phyper(             range, nPlus , nMinus,     i, lower.tail = FALSE, log.p = TRUE)
    }
    prevValue <- range[min(which(hyperVector < pMax))]
    criticalValues[i] <- prevValue + 1 # the +1 is due to the offset of phyper with respect to fisher.test
    if (is.finite(prevValue) && prevValue >= nPlus) {
      break
    }
  }
  output <- tibble::tibble(Sum = iRange, Total = criticalValues) %>%
    dplyr::mutate_all(as.integer) %>%
    dplyr::filter(Sum >= Total & Total <= nPlus)
  output
}

#' This function prepares a genotype and a phenotype matrix based on an inputTab and a specified number of SNPs, numSNPs
#' If complement =  1 or 2, each genotype  is present twice, the second time in negated form.
#' If complement = -1 or 2, each phenotype is present twice, the second time in negated form.
prepareMatrices = function(inputTab, numSNPs, complement) {
  IDs <- inputTab %>%
    dplyr::pull("ID")
  inputTab %<>%
    dplyr::select(-ID) %>%
    as.matrix %>%
    magrittr::set_rownames(IDs)
  phenotypes <- inputTab[, -(1:numSNPs), drop = FALSE]
  stopifnot(all(!is.na(phenotypes)))
  # print(paste("There are", ncol(phenotypes), "phenotypes over", nrow(phenotypes), "patients"))
  genotypes <- inputTab[, 1:numSNPs, drop = FALSE]
  if (complement > 0) {
    nGenotypes <- matrix(!genotypes, nrow(genotypes), dimnames = list(rownames(genotypes), paste0("NOT_", colnames(genotypes))))
    genotypes %<>%
      cbind(nGenotypes)
  }
  if (complement > 1 || complement < 0) {
    nPhenotypes <- matrix(!phenotypes, nrow(phenotypes), dimnames = list(rownames(phenotypes), paste0("NOT_", colnames(phenotypes))))
    phenotypes %<>%
      cbind(nPhenotypes)
  }
  output <- list(phenotypes = phenotypes, genotypes = genotypes)
  output
}

#' This function prepares the coefficients and the extra scores based on a specified genotype matrix and objective function.
prepareCoeffsAndExtraScores = function(genotypes, objective, sumYs, ns, TPscore = 1, TNscore = 1, FPscore = -1, FNscore = -1) {
  coeffsTrue  <- rep(1, ncol(genotypes))
  coeffsAgree <- rep(1, ncol(genotypes))
  if (objective == "agreement") {
    coeffsTrue  <- rep(TNscore - FPscore, ncol(genotypes))
    coeffsAgree <- rep(FPscore + FNscore - TPscore - TNscore, ncol(genotypes))
  } else if (objective == "covariance") {
    coeffsTrue  <- sumYs
    coeffsAgree  <- -ns
  }
  extraScores <- rep(0, ncol(genotypes))
  if (objective == "agreement") {
    extraScores <- TNscore * ns + (FNscore - TNscore) * sumYs
  }
  output <- list(coeffsTrue = coeffsTrue, coeffsAgree = coeffsAgree, extraScores = extraScores)
  output
}

#' This function reduces an input consisting of a phenotype matrix and a matrix of objective function vectors (same # of rows)
#' Optionally, it can take a matrix of genotypes and if that input is not NULL, also return a matrix of multiplicities extraDef
#' The return value is a list containing:
#' - the reduced phenotype and objective matrices
#' - the status of each row and column: FALSE if it has been removed, TRUE if it has been preserved
#' - the vector of extra scores to be added to the objective function
#' - optionally, a matrix of extra definitions specifying the multiplicity and total weight of each surviving position
reducePhenotypesAndObjectives = function(phenotypes, objVectors, genotypes = NULL) {
  stopifnot(nrow(objVectors) == nrow(phenotypes))
  rowStatuses <- rep(TRUE, nrow(phenotypes))
  colStatuses <- rep(TRUE, ncol(phenotypes))
  extraScores <- rep(0       , ncol(objVectors))
  if (!is.null(genotypes)) {
    extraDefS <- matrix(0, ncol(objVectors), 1 + nrow(objVectors), dimnames = list(paste0("S", 1:ncol(objVectors)), c()))
    extraDefT <- matrix(0, ncol(objVectors), 1 + nrow(objVectors), dimnames = list(paste0("T", 1:ncol(objVectors)), c()))
  }
  rSums <- rowSums(phenotypes)
  allTrueRows  <- which(rSums == ncol(phenotypes))
  allFalseRows <- which(rSums == 0)
  if (!is.null(genotypes) && length(allTrueRows) > 0) {
    extraDefT[, 1] <- colSums( genotypes[allTrueRows, , drop = FALSE], na.rm = TRUE)
    extraDefS[, 1] <- colSums(!genotypes[allTrueRows, , drop = FALSE], na.rm = TRUE) + extraDefT[, 1]
    extraScores    <- colSums(objVectors[allTrueRows, , drop = FALSE], na.rm = TRUE)
  }
  allConstantRows <- c(allTrueRows, allFalseRows)
  if (length(allConstantRows) > 0) {
    rowStatuses[allConstantRows] <- FALSE
    objVectors  %<>% magrittr::extract(-allConstantRows, , drop = FALSE)
    phenotypes  %<>% magrittr::extract(-allConstantRows, , drop = FALSE)
    if (!is.null(genotypes)) {
      genotypes %<>% magrittr::extract(-allConstantRows, , drop = FALSE)
    }
  }
  cSums <- colSums(phenotypes)
  allTrueCols  <- which(cSums == nrow(phenotypes))
  allFalseCols <- which(cSums == 0)
  allConstantCols <- c(allTrueCols, allFalseCols)
  if (length(allConstantCols) > 0) {
    colStatuses[allConstantCols] <- FALSE
    phenotypes  %<>% magrittr::extract(, -allConstantCols, drop = FALSE)
  }
  activeCols <- which(colStatuses)
  activeRows <- which(rowStatuses)
  uniqueRows <- extractUniqueRows(phenotypes)
  uniqueGroupsR <- uniqueRows[[2]]
  phenotypes <- uniqueRows[[1]]
  for (ind in seq_along(uniqueGroupsR)) {
    curGroup <- uniqueGroupsR[[ind]]
    rowStatuses[activeRows[curGroup]] <- ind * c(1L, rep(-1L, length(curGroup) - 1))
    if (!is.null(genotypes)) {
      extraDefT[, ind + 1] = colSums( genotypes[curGroup, , drop = FALSE], na.rm = TRUE)
      extraDefS[, ind + 1] = colSums(!genotypes[curGroup, , drop = FALSE], na.rm = TRUE) + extraDefT[, ind + 1]
    }
  }
  if (!is.null(genotypes)) {
    usableColumns <- 1:(length(uniqueGroupsR) + 1)
    extraDefS %<>% magrittr::extract(, usableColumns, drop = FALSE)
    extraDefT %<>% magrittr::extract(, usableColumns, drop = FALSE)
    extraDef  <- rbind(extraDefS, extraDefT)
  }
  objVectors <- do.call(rbind, lapply(uniqueGroupsR, function(x) {
    cur <- objVectors[x, , drop = FALSE]; res = colSums(cur, na.rm = TRUE); res[colSums(is.na(cur)) == nrow(cur)] = NA; res
  }))
  uniqueCols <- extractUniqueRows(t(phenotypes))
  uniqueGroupsC <- uniqueCols[[2]]
  phenotypes <- t(uniqueCols[[1]])
  for (ind in seq_along(uniqueGroupsC)) {
    curGroup <- uniqueGroupsC[[ind]]
    colStatuses[activeCols[curGroup]] <- ind * c(1L, rep(-1L, length(curGroup) - 1))
  }
  print(paste("After reduction, there are", ncol(phenotypes), "phenotypes over", nrow(phenotypes), "patients"))
  print(paste("After reduction, there are", ncol(objVectors),  "genotypes over", nrow(objVectors), "patients"))
  output <- list(phenotypes = phenotypes, objVectors = objVectors, rowStatuses = rowStatuses, colStatuses = colStatuses,
                 extraScores = extraScores, extraDef = extraDef)
  output
}

#' NB: the extreme value is in log-space
prepareExtremeValues = function(extremeValue, type, singleValues, singleBestRatio = Inf, useTighterBound = FALSE) {
  extremeValues <- rep(extremeValue, length(singleValues))
  if (!is.null(extremeValue)) {
    if (!is.na(singleBestRatio)) {
      altBounds <- singleValues %>%
        magrittr::subtract(log(singleBestRatio))
      altFunction <- ifelse(!useTighterBound, pmax, pmin)
      extremeValues %<>%
        altFunction(altBounds, na.rm = TRUE)
    }
  }
  extremeValues
}

#' Prepares the upper cutoffs based on the input, type and objective function
prepareBoundValues = function(extremeValues, type, objective, ns, sumYs, coeffsTrue, coeffsAgree) {
  startBound <- NA
  if (any(is.null(extremeValues))) {
    startBound <- -DELTA # subtract DELTA to ensure that the value is at least trivial + DELTA
  }
  boundValues <- rep(startBound, length(ns))
  if (type == "CNF" && objective == "agreement") {
    boundValues %<>%
      magrittr::add(ns * coeffsTrue + sumYs * coeffsAgree)
  }
  boundValues
}

#' This function prepares the extra information relating to the boundary
prepareExtras = function(extraDefinitions, ind, goodInds, extremeValue, ns, sumYs, baseFilename, index = 0L) {
  stopifnot(extremeValue <= 0)
  curDefinitions <- extraDefinitions[c(paste0("S", ind), paste0("T", ind)), c(1, 1 + goodInds), drop = FALSE]
  curRFile <- paste0(ifelse(is.null(baseFilename), "Test", baseFilename), "I", ind, "C", index, "Boundary", extremeValue, ".csv")
  curPBounds    <- computePValueBounds(N = ns[ind], nPlus = sumYs[ind], pMax = extremeValue)
  curBoundaries <- callImaiIri(curPBounds, dryRun = FALSE, width = WIDTH, fname = curRFile)[[1]]
  if (file.exists(curRFile)) {
    file.remove(curRFile)
  }
  curExtras <- list(B = curBoundaries, D = curDefinitions)
  curExtras
}
