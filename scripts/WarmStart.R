### This function creates an MST file containing the solution to a problem containing only binary decision variables
### based on the list of variable names, the list of true indices, and the (optional) objective function value.
### If type == "LCP", both binary and general variables are allowed
makeSolutionFile = function(colnames, trueVars, optValue = 0, filename = "Test.mst", type = "CNF") {
  Lines <- rep("", 26 + length(colnames))
  Lines[1]  <- '<?xml version = "1.0" standalone="yes"?>'
  Lines[2]  <- '<CPLEXSolution version="1.2">'
  Lines[3]  <- ' <header'
  Lines[4]  <- paste0('problemName="', str_replace(filename, '.mst', '.lp'), '"')
  Lines[5]  <- 'solutionName="incumbent"'
  Lines[6]  <- 'solutionIndex="-1"'
  Lines[7]  <- paste0('objectiveValue="', optValue, '"')
  Lines[8]  <- 'solutionTypeValue="3"'
  Lines[9]  <- 'solutionStatusValue="101"'
  Lines[10] <- 'solutionStatusString="integer optimal solution"'
  Lines[11] <- 'solutionMethodString="mip"'
  Lines[12] <- 'primalFeasible="1"'
  Lines[13] <- 'dualFeasible="1"'
  Lines[14] <- 'MIPNodes="219"'
  Lines[15] <- 'MIPIterations="65222"'
  Lines[16] <- 'writeLevel="1"/>'
  Lines[17] <- '<quality'
  Lines[18] <- 'epInt="1e-05"'
  Lines[19] <- 'epRHS="1e-06"'
  Lines[20] <- 'maxIntInfeas="0"'
  Lines[21] <- 'maxPrimalInfeas="0"'
  Lines[22] <- 'maxX="40"'
  Lines[23] <- 'maxSlack="2"/>'
  Lines[24] <- '<variables>'
  numVars <- length(colnames)
  values <- rep(0L, numVars)
  if (type == "LCP") {
    values[as.integer(names(trueVars))] <- trueVars
  } else {
    values[trueVars] <- 1L
  }
  Lines[25:(length(Lines) - 2)] <- paste0('<variable name="', colnames, '" index="', 0:(numVars - 1), '" value="', values, '"/>')
  Lines[length(Lines) - 1] <- '</variables>'
  Lines[length(Lines)] <- '</CPLEXSolution>'
  f = file(filename)
  writeLines(Lines, f)
  close(f)
}

### This auxiliary function prepares a file with an initial solution to be used for a logical formula (CNF or DNF) with an ILP.
makeInitialSolution = function(startSol, filename, type, K, phenotypes) {
  n <- nrow(phenotypes)
  p <- ncol(phenotypes)
  UInds <- matrix(1:(K * p), nrow = K)
  lastUInd <- K * p
  PInds <- matrix(lastUInd + (1:(K * n)), nrow = K)
  lastPInd <- lastUInd + K * n
  pInds <- lastPInd + (1:n)
  coln <- c(paste0("U", as.vector(outer(1:K, 1:p, pasteI))), paste0("P", as.vector(outer(1:K, 1:n, pasteI))), paste0("p", 1:n))
  usedVars <- startSol$usedVars
  numSols <- ncol(usedVars)
  if (numSols < K) { ### complete the matrix by duplicating the last column
    usedVars %<>%
      cbind(usedVars[, rep(numSols, K - numSols)])
  }
  transMap <- match(rownames(usedVars), colnames(phenotypes))
  trueUInds <- c()
  truePInds <- c()
  fullPheno <- rep(ifelse(type == "CNF", TRUE, FALSE), n)
  aggregateFunction <- ifelse(type == "CNF", magrittr::and, magrittr::or)
  for (index in 1:K) {
    curInds <- which(usedVars[, index] == 1)
    transInds <- transMap[curInds]
    trueUInds %<>%
      c(UInds[index, transInds])
    curPheno <- apply(phenotypes[, transInds, drop = FALSE], 1, ifelse(type == "CNF", any, all))
    truePInds %<>%
      c(PInds[index, curPheno])
    fullPheno %<>%
      aggregateFunction(curPheno)
  }
  truepInds <- pInds[fullPheno]
  allTrueVars <- c(trueUInds, truePInds, truepInds)
  solFname <- str_replace(filename, ".lp", ".mst")
  output <- makeSolutionFile(coln, allTrueVars, optValue = startSol$optimum, filename = solFname, type = type)
  output
}
