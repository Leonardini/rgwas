#' Compute the optimum value an output shoud have in order to compare it to the actual optimum value
computeOptValue = function(myTab, objective, coeffsTrue, coeffsAgree, ind, TPscore = 1, TNscore = 1, FPscore = -1, FNscore = -1) {
  initTab <- matrix(0, 2, 2, dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE")))
  initTab[rownames(myTab), colnames(myTab)] %<>%
    magrittr::add(myTab)
  TP <- initTab["TRUE" , "TRUE" ]
  FP <- initTab["FALSE", "TRUE" ]
  FN <- initTab["TRUE" , "FALSE"]
  TN <- initTab["FALSE", "FALSE"]
  checkOpt <- ifelse(objective == "agreement", sum(c(TP, FP, FN, TN) * c(TPscore, FPscore, FNscore, TNscore)),
                   -(coeffsTrue[ind] * (TP + FP) + coeffsAgree[ind] * TP))
  output <- list(checkOpt = checkOpt, TP = TP, FP = FP, FN = FN, TN = TN)
  output
}

#' Compute the p-value of the Cochran-Mantel-Haenszel test, estimating the association of genotype to phenotype
#' stratifying by a secondary phenotype; exact = TRUE is a Fisher exact test analog; correct = TRUE is a continuity correction
computeStratifiedPVal = function(genotype, phenotypeC, phenotypeS, exact = FALSE, correct = FALSE) {
  microTabs <- tibble::tibble(A = genotype, B = phenotypeC, C = phenotypeS) %>%
    dplyr::filter(!is.na(A) & !is.na(B) & !is.na(C)) %>%
    dplyr::group_by(A, B, C) %>%
    dplyr::summarize(N = dplyr::n(), .groups = "keep") %>%
    dplyr::ungroup()
  fullTabs <- array(0, c(2, 2, 2), dimnames = list(c("FALSE", "TRUE"), c("FALSE", "TRUE"), c("FALSE", "TRUE")))
  for (ind in seq_len(nrow(microTabs))) {
    curRow <- microTabs %>%
      dplyr::slice(ind)
    fullTabs[as.character(curRow$A), as.character(curRow$B), as.character(curRow$C)] <- curRow$N
  }
  extraStat <- log(mantelhaen.test(fullTabs, alternative = "g", exact = exact, correct = correct)$p.value)
  extraStat
}

#' Create and solve an ILP (integer linear program) for maximizing an objective function over input phenotypes
#' The arguments are similar to the function below, except that the genotype is not supplied and phenotypes must be a matrix
#' The objVector specifies the vector of objective function values to apply to each entry of the phenotype during optimization
#' The boundValue parameter specifies a bound on the objective function; not achieving this bound makes the problem infeasible
createAndSolveILP = function(phenotypes, objVector, type = "CNF", K = 3, L = 3, filename = "Test.lp", boundValue = NA,
                                extraConstraints = NULL) {
  if (length(extraConstraints) == 2 && any(is.na(extraConstraints[[1]]))) {
    return(list(usedVars = NA, optimum = NA, status = HOPELESS, time = 0, gap = NA))
  }
  n <- nrow(phenotypes)
  p <- ncol(phenotypes)
  goodEntries <- which(phenotypes == ifelse(type == "CNF", 1, 0), arr.ind = TRUE)
  goodEntries %<>%
    tibble::as_tibble() %>%
    dplyr::arrange(row, col)
  N <- nrow(goodEntries)
  numVar <- K * (n + p) + n
  numConst <- K * (N + 2 * n + 2) + n
  UInds <- matrix(1:(K * p), nrow = K)
  lastUInd <- K * p
  PInds <- matrix(lastUInd + (1:(K * n)), nrow = K)
  lastPInd <- lastUInd + K * n
  pInds <- lastPInd + (1:n)
  lastOrRow <- K * (n + N)
  andRows <- lastOrRow + (1:(n * (K + 1)))
  lastAndRow <- lastOrRow + n * (K + 1)
  LRows <- lastAndRow + (1:K)
  FRows <- lastAndRow + K + (1:K)
  coln <- c(paste0("U", as.vector(outer(1:K, 1:p, function(x, y) {paste(y, x, sep = "I")}))),
            paste0("P", as.vector(outer(1:K, 1:n, function(x, y) {paste(y, x, sep = "I")}))), paste0("p", 1:n))
  phenoSums <- rowSums(phenotypes * ifelse(type == "CNF", 1, -1) + ifelse(type == "CNF", 0, 1))
  numNonZeros <- K * (3 * N + 4 * n + 2 * p) + n
  if (!is.null(extraConstraints)) {
    numBoundaries <- nrow(extraConstraints[[1]])
    if (numBoundaries == 0) {
      return(list(usedVars = matrix(0, 0, 0), optimum = NA, status = INFEASIBLE, time = 0, gap = NA))
    }
    numSeg <- numBoundaries - 1
    stopifnot(nrow(extraConstraints[[2]]) == 2)
    numSTCoeffs <- ncol(extraConstraints[[2]])
    numExtraVar <- 2 * numSeg + 3
    numVar %<>%
      magrittr::add(numExtraVar)
    numExtraConst <- 2 * numSeg + 6
    numConst %<>%
      magrittr::add(numExtraConst)
    numExtraNonZeros <- 2 * numSTCoeffs + 8 * numSeg + 4
    numNonZeros %<>%
      magrittr::add(numExtraNonZeros)
    coln %<>%
      c("Sum", "Total", "Z", paste0(rep("Z", numSeg), seq_len(numSeg)), paste0(rep("S", numSeg), seq_len(numSeg)))
  }
  fullObjVector <- rep(0, numVar)
  fullObjVector[pInds] <- objVector
  rowInds <- rep(NA, numNonZeros)
  colInds <- rep(NA, numNonZeros)
  values  <- rep(NA, numNonZeros)
  pos <- 0
  rhs <- rep(0, numConst)
  for (iter in 1:K) {
    curRows <- N * (iter - 1) + (1:N)
    curRowsN <- (N * K) + n * (iter - 1) + (1:n)
    curRowsNext <- rep(curRowsN, phenoSums)
    curIndsP <- UInds[iter, goodEntries$col]
    curIndsN <- PInds[iter, goodEntries$row]
    curRange <- pos + (1:(3 * N + n))
    rowInds[curRange] <- c(curRows,   curRows,                         curRowsN,                        curRowsNext)
    colInds[curRange] <- c(curIndsP,  curIndsN,                        PInds[iter,],                    curIndsP)
    values [curRange] <- c(rep(1, N), rep(1 - 2 * (type == "CNF"), N), rep(2 * (type == "CNF") - 1, n), rep(-1, N))
    if (type == "DNF") {
      rhs[curRows]  <- 1
      rhs[curRowsN] <- -1
    }
    pos <- pos + 3 * N + n
  }
  endRange <- (pos + 1):(numNonZeros - ifelse(!is.null(extraConstraints), numExtraNonZeros, 0) - 2 * p * K)
  Sz <- n * K
  nextIndsP <- rep(pInds, K)
  nextIndsN <- as.vector(t(PInds))
  expAndRows <- rep(andRows[-(1:Sz)], K)
  rowInds[endRange] <- c(andRows[1:Sz], andRows[1:Sz], expAndRows, andRows[-(1:Sz)])
  colInds[endRange] <- c(nextIndsP,     nextIndsN,     nextIndsN,  pInds           )
  values [endRange] <- c(rep(1, Sz),    rep(-1, Sz),   rep(1, Sz), rep(-1, n)      )
  if (type == "CNF") {
    rhs[lastOrRow + n * K + (1:n)] <- K - 1
  }
  if (type == "DNF") {
    values[(3 * N + n) * K + (1:(3 * Sz + n))] %<>%
      magrittr::multiply_by(-1)
  }
  lastRange <- numNonZeros - ifelse(!is.null(extraConstraints), numExtraNonZeros, 0) - (2 * p * K) + 1:(p * K)
  rowInds[lastRange] <- rep(LRows, p)
  colInds[lastRange] <- as.vector(UInds)
  values [lastRange] <- rep(1, p * K)
  rhs[lastAndRow + (1:K)] <- L
  finalRange <- numNonZeros - ifelse(!is.null(extraConstraints), numExtraNonZeros, 0) - (p * K) + 1:(p * K)
  rowInds[finalRange] <- rep(FRows, p)
  colInds[finalRange] <- as.vector(UInds)
  values [finalRange] <- rep(-1, p * K)
  rhs[lastAndRow + K + (1:K)] <- -1
  if (!is.null(extraConstraints)) {
    extraRange <- K * (3 * N + 4 * n + 2 * p) + n + 1:(numExtraNonZeros)
    rowInds[extraRange] <- c(rep(lastAndRow + 2 * K + (1:2), each = numSTCoeffs),
                             rep(lastAndRow + 2 * K + 3, 2),
                             rep(lastAndRow + 2 * K + 4, numSeg),
                             rep(lastAndRow + 2 * K + 4 + seq_len(numSeg), each = 2),
                             rep(lastAndRow + 2 * K + 4 + numSeg + (1:2), each = 2 * numSeg + 1),
                             rep(lastAndRow + 2 * K + 4 + numSeg + 2 + seq_len(numSeg)))
    indS <- numVar - numExtraVar + 1
    indT <- numVar - numExtraVar + 2
    indZ <- numVar - numExtraVar + 3
    indZi <- numVar - numExtraVar + 3 + seq_len(numSeg)
    indSi <- numVar - numExtraVar + 3 + numSeg + seq_len(numSeg)
    colInds[extraRange] <- c(pInds, indS,
                             pInds, indT,
                             indS,  indZ,
                             indZi,
                             as.vector(matrix(rbind(indSi, indZi), ncol = 2)), # this "interlaces" S and Z indices like zip
                             indT, indZi, indSi, # x-value constraint
                             indZ, indZi, indSi, # y-value constraint
                             indSi)
    Xvalues <- extraConstraints[[1]] %>%
      dplyr::pull(1)
    deltaXs <- diff(Xvalues)
    Yvalues <- extraConstraints[[1]] %>%
      dplyr::pull(2)
    deltaYs <- diff(Yvalues)
    values [extraRange] <- c(-extraConstraints[[2]][1, -1], 1,
                             -extraConstraints[[2]][2, -1], 1,
                             1, -1,
                             rep(1, numSeg),
                             rep(c(1, -1), numSeg),
                             c(-1, Xvalues %>% head(-1), deltaXs),
                             c(-1, Yvalues %>% head(-1), deltaYs),
                             rep(-1, numSeg))
    rhs[lastAndRow + 2 * K + (1:4)] <- c(extraConstraints[[2]][1, 1], extraConstraints[[2]][2, 1], 0, as.double(numSeg > 0))
    if (numSeg == 0) {
      rhs[lastAndRow + 2 * K + (5:6)] <- c(-Xvalues[1], -Yvalues[1])
    }
  }
  Control = makeControlLines()
  if (!is.na(boundValue)) {
    numConst = numConst + 1
    rowInds  = c(rowInds, rep(numConst, numVar))
    colInds  = c(colInds, 1:numVar)
    values   = c(values, fullObjVector)
    rhs      = c(rhs, boundValue)
  }
  Mat = Matrix::sparseMatrix(i = rowInds, j = colInds, x = values, dims = c(numConst, numVar), index1 = TRUE)
  varTypes = rep("B", numVar)
  constDir = rep("L", numConst)
  if (!is.null(extraConstraints)) {
    varTypes[(numVar - numExtraVar + 1):numVar] = c(rep("I", 2), "C", rep("B", numSeg), rep("C", numSeg))
    constDir[numConst - (numExtraConst + !is.na(boundValue)) + c(1, 2, 4, numSeg + 5, numSeg + 6)] = "E"
  }
  UB = as.integer(ifelse(varTypes == "I", n, 1))
  solution = Rcplex::Rcplex(cvec = fullObjVector, Amat = Mat, bvec = rhs, ub = UB, control = Control, objsense = "min", sense = constDir, vtype = varTypes)
  solution$xopt %<>%
    magrittr::set_names(coln)
  Rcplex::Rcplex.close()
  output = extractSolution(solution)
  if (!all(is.na(output$usedVars))) {
    rownames(output$usedVars) <- colnames(phenotypes)[seq_len(nrow(output$usedVars))]
  }
  output
}
