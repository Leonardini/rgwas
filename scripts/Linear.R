### Specific options for LCP
QUAD_CP_STRAT= ifelse(OPTIMIZATION,  NA, NA)
QUAD_CP_TOL  = ifelse(OPTIMIZATION,  NA, NA) 

### This function creates and solves a MIQCP (mixed integer quadratically constrained program) to maximize an objective function
### The rMatrix indicates the constraint matrix, the objVector specifies the vector of objective function values, S is the scale
### Note: In the latest version, rMatrix is an upper triangular matrix R such that the "true" coefficient matrix is t(R) %*% R!
### If K is not NULL, there is an additional restriction on the sparsity (0-norm) of the output vector: it must not exceed K
### If M is not NULL, the binary variables defining sparsity use its value as the big-M (upper bound on the output coefficients)
### The filename used to write the program down may be specified, and an optional starting solution may be specified in startSol
### If boundValue is not NA, it is used to lower bound the found solution
createAndSolveMIQCP = function(rMatrix, objVector, S = 1, K = 3, M = 10, filename = "Test.lp", startSol = NULL, boundValue = NA) {
  Dim <- nrow(rMatrix)
  stopifnot(length(objVector) == Dim)
  numLines <- 8 + 6 * Dim
  Lines <- rep("", numLines)
  Lines[1] <- "Minimize"
  Lines[2] <- paste0(paste(-objVector/S, paste0("U", 1:Dim)), collapse = " +")
  Lines[3] <- "Subject to"
  Lines[4] <- paste0("[ ", paste0("Y", 1:Dim, "^2", collapse = " +"), " ] <= ", 1) ### and now, we express the Y's via the U's
  for (ind in 1:Dim) {
    curRange <- ind:Dim
    curRow <- rMatrix[ind, curRange]
    Lines[4 + ind] <- paste0(paste0(curRow, " U", curRange, collapse = " +"), " - Y", ind, " = 0")
  }
  start <- 4 + Dim
  Lines[start + (1:Dim)]       <- paste0("U", 1:Dim, " - ", M, " Z", 1:Dim, " <= 0")
  Lines[start + Dim + (1:Dim)] <- paste0("-U", 1:Dim, " - ", M, " Z", 1:Dim, " <= 0")
  start2 <- start + 2 * Dim + 1
  Lines[start2] <- paste0(paste0("Z", 1:Dim, collapse = " + "), " <= ", K)
  Lines[start2 + 1] <- "Bounds"
  Lines[start2 + 1 + (1:Dim)] <- paste0("U", 1:Dim, " free")
  start3 <- start2 + 1 + Dim
  Lines[start3 + (1:Dim)] <- paste0("Y", 1:Dim, " free")
  start4 <- start3 + Dim
  Lines[start4 + 1] <- "Binary"
  Lines[start4 + 1 + (1:Dim)] <- paste0("Z", 1:Dim)
  Lines[length(Lines)] <- "End"
  Lines %<>%
    str_replace_all("\\+-", "-")
  if (!is.null(startSol)) {
    solFname <- str_replace(filename, ".lp", ".mst")
    usedVars <- startSol$usedVars
    numSols <- nrow(usedVars)
    coln <- colnames(rMatrix)
    transMap <- match(rownames(usedVars), coln)
    trueUInds <- usedVars[,1]
    names(trueUInds) <- transMap
    trueZInds <- rep(1L, length(transMap))
    names(trueZInds) <- Dim + transMap
    trueInds <- c(trueUInds, trueZInds)
    fullColn <- as.vector(t(outer(c("U","Z"), 1:Dim, paste0)))
    makeSolutionFile(fullColn, trueInds, optValue = startSol$optimum, filename = solFname, type = "LCP")
  } else {
    solFname <- NULL
  }
  f = file(filename, 'w')
  writeLines(Lines, f)
  close(f)
  output <- runAndParse(filename, prefix2 = NULL, boundValue = boundValue)
  if (!all(is.na(output$usedVars))) {
    rownames(output$usedVars) <- colnames(rMatrix)[seq_len(nrow(output$usedVars))]
  }
  output
}

### Auxiliary function to process an LCP problem given an input genotype G, phenotype P, scaling factor S, and other parameters.
processLCP = function(G, P, S, ind, initBound, startSol, extremeValue, K, L, filename) {
  Pc <- scale(P, center = TRUE, scale = FALSE)
  Gc <- scale(G, center = TRUE, scale = FALSE)
  QR <- qr(Pc)
  rk <- QR$rank
  print(paste("There are", rk, "linearly independent columns in the phenotype matrix"))
  rMatrix <- qr.R(QR)[1:rk, 1:rk]
  goodCols <- QR$pivot[1:rk]
  objVector <- t(Gc) %*% Pc[, goodCols]
  if (is.na(initBound)) {
    intermedResult <- backsolve(rMatrix, forwardsolve(rMatrix, t(objVector), upper.tri = TRUE, transpose = TRUE))
    upperBound <- sqrt(objVector %*% intermedResult/sum(Gc^2))
  } else {
    upperBound <- initBound
  }
  print(paste("The maximum value we could possibly hope for (without sparsity) is", upperBound))
  if (!is.na(upperBound) && upperBound >= extremeValue) {
    if (K >= ncol(Pc)) { 
      print(paste("Special case: we may use all", ncol(Pc), "phenotypes, so CCA suffices for the optimization"))
      Z <- cancor(Gc, Pc, xcenter = FALSE, ycenter = FALSE)
      stopifnot(near(Z$cor, upperBound))
      usedVars <- Z$ycoef[colnames(Pc), 1, drop = FALSE]
      usedVars <- usedVars/max(abs(usedVars))
      curSol <- list(usedVars = usedVars, optimum = upperBound, status = OPTIMAL, time = 0, gap = 0)
    } else {
      curSol <- createAndSolveMIQCP(rMatrix, objVector, S = S, K = K, M = L, filename = filename, startSol = startSol, 
                                    boundValue = -extremeValue) ### use min instead of max
    }
  } else {
    print(paste("The current problem has no shot at getting a correlation at or above", extremeValue))
    curSol <- list(usedVars = NA, optimum = NA, status = INFEASIBLE, time = 0, gap = NA)
  }
  output <- list(solution = curSol, bound = upperBound)
  output
}