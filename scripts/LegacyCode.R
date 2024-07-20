### Packages that are no longer used
library(expm)
library(FRACTION)
library(geometry)
library(Hmisc)
library(igraph)
library(magic)
library(pracma)

### SOME FUNCTIONS FROM UTILITIES.R
getFreqs = function(inputFile) {
  ext <- str_sub(inputFile, -4)
  readFunction <- ifelse(ext == '.tsv', read_tsv, read_csv)
  writeFunction <- ifelse(ext == ".tsv", write_tsv, write_csv)
  inputTab <- readFunction(inputFile, col_types = cols(ID = "c", .default = "l")) %>%
    select(-ID)
  outputTab <- bind_cols(name = colnames(inputTab), Sum = as.integer(colSums(inputTab)), N = nrow(inputTab)) %>%
    mutate(frequency = round(Sum/N, 3))
  outputFile <- str_replace(inputFile, ext, paste0("Frequencies", ext))
  writeFunction(outputTab, outputFile)
  outputTab
}

deleteZeroRows = function(inputMatrix) {
  inputMatrix %<>% 
    subset(rowSums(inputMatrix) > 0)
  inputMatrix
}

deleteFullColumns = function(inputMatrix) {
  goodCols <- (colSums(inputMatrix) < nrow(inputMatrix))
  inputMatrix %<>% 
    subset(select = goodCols)
  inputMatrix
}

findSubsetColumns = function(inputMatrix, inputVector) {
  M <- nrow(inputMatrix)
  N <- ncol(inputMatrix)
  L <- length(inputVector)
  stopifnot(L == M)
  sumVector <- sum(inputVector)
  goodCols <- (colSums(inputMatrix | matrix(inputVector, M, N)) == sumVector)
  goodCols
}

findSupersetColumns = function(inputMatrix, inputVector) {
  M <- nrow(inputMatrix)
  N <- ncol(inputMatrix)
  L <- length(inputVector)
  stopifnot(L == M)
  sumVector <- sum(inputVector)
  goodCols <- (colSums(inputMatrix & matrix(inputVector, M, N)) == sumVector)
  goodCols
}

reduceGroupsByInclusion = function(groups, inputMatrix) {
  M <- nrow(inputMatrix)
  N <- ncol(inputMatrix)
  setNames <- colnames(inputMatrix)
  stopifnot(!is.null(setNames))
  L <- length(groups)
  redGroups <- vector("list", L)
  for (ind in seq_len(L)) {
    curGroup <- groups[[ind]]
    curNames <- setNames[curGroup]
    curSubmatrix <- t(inputMatrix[, curGroup, drop = FALSE])
    redSubmatrix <- t(extractMinimalRows(curSubmatrix))
    redGroups[[ind]] <- curGroup[curNames %in% colnames(redSubmatrix)]
  }
  redGroups
}

### Set J is dominated if there is an I with setsP[J] <= setsP[I], setsN[J] >= setsN[I], and at least one is strict.
### Note that duplicates are automatically removed as well!
removeDominated = function(setsP, setsN) {
  N <- ncol(setsP)
  origInds <- 1:N
  stopifnot(ncol(setsN) == N)
  M1 <- nrow(setsP)
  ### start by pruning for equality in setsP, with strict inequality in setsN
  U <- extractUniqueRows(t(setsP), repeats = TRUE)
  groups <- U[[2]]
  groups <- reduceGroupsByInclusion(groups, setsN)
  keepCols <- unlist(groups)
  origInds <- keepCols
  setsP <- subset(setsP, select = keepCols)
  setsN <- subset(setsN, select = keepCols)
  redN <- ncol(setsP)
  status <- rep(FALSE, redN)
  for (ind in 1:redN) {
    curPCol <- setsP[,ind]
    curNCol <- setsN[,ind]
    curCands <- which(findSubsetColumns(setsP, curPCol))
    curConfirmed <- setdiff(curCands[findSupersetColumns(setsN[, curCands, drop = FALSE], curNCol)], ind)
    status[curConfirmed] <- TRUE
  }
  origInds <- origInds[!status]
  setsP <- subset(setsP, select = !status)
  setsN <- subset(setsN, select = !status)
  output <- list(setsP, setsN, origInds)
  output
}

### Auxiliary function for creating a named list with the given names of the named objects inside the caller function 
makeNamedList = function(listOfNames) {
  output <- mget(listOfNames, envir = parent.frame())
  names(output) <- listOfNames
  output
}

executeArgs = function(String, splitBy = ",") {
  altStrings <- String %>%
    str_replace_all("=", "<-") %>%
    str_split(splitBy) %>%
    extract2(1)
  for (ind in seq_along(altStrings)) {
    part <- altStrings[[ind]] %>%
      str_trim
    if (str_detect(part, "<-")) {
      eval.parent(str2lang(part))
    }
  }
}

executeLoop = function(String) {
  altString <- String %>%
    str_remove("for") %>%
    str_remove("\\(") %>%
    str_remove("\\)") %>%
    str_remove("\\{") %>%
    str_trim(side = "both") %>%
    str_replace("in", "<-") %>%
    str_c("[1]")
  eval.parent(str2lang(altString))
}

negateCoefficients = function(Formula, joinString = " \\+ ") {
  splitFormula <- str_split(Formula, joinString)[[1]]
  for (ind in seq_along(splitFormula)) {
    part <- splitFormula[[ind]]
    if (str_detect(part, '\\(-')) {
      splitFormula[[ind]] <- str_replace(part, '\\(-', '(')
    } else {
      splitFormula[[ind]] <- str_replace(part, '\\(', '(-')
    }
  }
  negFormula <- str_c(splitFormula, collapse = str_remove(joinString, '\\\\'))
  negFormula
}

makeRowColumnGraph = function(inputMatrix, directed = FALSE) {
  edgeList <- which(inputMatrix, arr.ind = TRUE)
  edgeList[,1] <- paste0("R", edgeList[,1])
  edgeList[,2] <- paste0("C", edgeList[,2])
  initGraph <- graph.data.frame(edgeList, directed = directed)
  initGraph
}

getConnectedComponents = function(inputMatrix, orderIncreasing = FALSE) {
  initGraph <- makeRowColumnGraph(inputMatrix)
  dimnames(inputMatrix) <- list(1:nrow(inputMatrix), 1:ncol(inputMatrix))
  clusters <- components(initGraph)
  numClusters <- clusters$no
  sizeClusters <- clusters$csize
  output <- vector("list", numClusters)
  for (ind in 1:numClusters) {
    curVerts <- V(initGraph)[clusters$membership == ind]$name
    rowVerts <- (substr(curVerts, 1, 1) == "R")
    curRows <- as.integer(substr(curVerts[rowVerts], 2, nchar(curVerts[rowVerts])))
    curCols <- as.integer(substr(curVerts[!rowVerts], 2, nchar(curVerts[!rowVerts])))
    output[[ind]] <- inputMatrix[curRows, curCols, drop = FALSE]
  }
  if (orderIncreasing) {
    myOrder <- order(sizeClusters)
    output <- output[myOrder]
  }
  output
}

### SOME FUNCTIONS FROM SIMULATION.R

### Assumptions: the file only contains the constraints on the phenotypes and choices, and the constraint on T, plus formatting
parseCPLEXInput = function(inputFile = "TrickyTest.rtf") {
  Lines <- readLines(inputFile)
  Lines <- tail(Lines, -8)
  Lines <- str_remove(Lines, "R[0-9]+\\:")
  Lines <- str_remove(Lines, "\\\\")
  Lines <- str_trim(Lines)
  lastIneqLine <- max(which(str_detect(Lines, "<= 0")))
  genoLines <- Lines[(lastIneqLine + 1):length(Lines)]
  genotype <- str_extract_all(genoLines, "p([0-9]+)", FALSE) %>%
    unlist %>%
    str_remove("p") %>%
    parse_integer
  phenoLines <- Lines[str_detect(Lines, paste0(PREFIX2, 1, " <= ", 0))]
  phenoLines %<>%
    str_remove("<= 0") %>%
    str_remove_all(paste0(PREFIX2, "1")) %>%
    str_remove_all("\\-") %>%
    str_remove_all("\\+") %>%
    str_squish
  goodLines <- which(
    phenoLines %>%
      str_count(" ")
    == 1)
  phenoLines %<>%
    magrittr::extract(goodLines)
  allInds <- str_split_fixed(phenoLines, " ", n = 2) %>%
    as_tibble(.name_repair = "minimal")
  colnames(allInds) <- c("col", "row")
  allInds %<>%
    mutate_all(~parse_integer(str_sub(., start = 2))) %>%
    distinct
  M <- max(allInds$row)
  N <- max(allInds$col)
  phenoMatrix <- matrix(FALSE, M, N)
  phenoMatrix[cbind(allInds$row, allInds$col)] <- TRUE
  genoVector <- matrix(FALSE, M, 1)
  genoVector[genotype, ] <- TRUE
  output <- list(SNP = genoVector, phenotype = phenoMatrix)
  output
}

#### Alternative test: see how much the running time improves on TrickyTest.mps
checkImprovement = function() {
  Z <- parseCPLEXInput()
  genotype0 <- Z[[1]]
  phenotypes0 <- Z[[2]]
  n0 <- nrow(phenotypes0)
  p0 <- ncol(phenotypes0)
  testTab0 <- cbind(genotype0, phenotypes0)
  colnames(testTab0) <- c("SNP1", paste0("p", 1:p0))
  fullTab0 <- as_tibble(testTab0) %>% bind_cols(tibble(ID = paste0("Patient", 1:n0))) %>% select(ID, everything())
  allSolutions <- vector("list", 2)
  allSolutions[[1]] <- createAndSolveILPs(fullTab0, numSNPs = 1, type = "CNF", complement = FALSE, K = 3, objective = "agreement")
  allSolutions[[2]] <- createAndSolveILPs(fullTab0, numSNPs = 1, type = "CNF", complement = FALSE, K = 3, objective = "covariance")
  allSolutions
}

### SOME FUNCTIONS FROM PREPROCESS.R
createPairs = function(inputFile ="CPLEXInputs/CPLEXInput_SelfReportedReduced_EU_100000_15.csv", numSNPs = 15, type = "CNF") {
  ext = str_sub(inputFile, -4)
  readFunction = ifelse(ext == ".csv", read_csv, read_tsv)
  Tab = readFunction(inputFile)
  allSNPs = Tab[, 1:(1 + numSNPs)]
  phenos = Tab %>%
    select(-(1:(1 + numSNPs))) %>%
    as.matrix()
  M = ncol(phenos)
  myFun = ifelse(type == "CNF", or, and)
  pairedP = combn(M, 2, function(x) { myFun(phenos[,x[1]], phenos[,x[2]]) } )
  colnames(pairedP) = combn(colnames(phenos), 2, function(x) {paste0(x[1], ifelse(type == "CNF", " + ", " * "), x[2])})
  fullPhenos = cbind(phenos, pairedP) %>%
    as_tibble()
  extTab = bind_cols(allSNPs, fullPhenos)
  write_csv(extTab, file = str_replace(inputFile, "ed", "edPaired"))
}

reduceFeatures = function(inputFile ="CPLEXInputs/CPLEXInput_SelfReported_EU_10000_15.csv", numSNPs = 15, SNPpos = 1) {
  ext = str_sub(inputFile, -4)
  readFunction = ifelse(ext == ".csv", read_csv, read_tsv)
  Tab = readFunction(inputFile)
  allSNPs = colnames(Tab)[1 + (1:numSNPs)]
  curSNPName = allSNPs[SNPpos]
  curSNP = Tab %>%
    select(all_of(curSNPName))
  phenotypes = Tab %>%
    select(-(1:(1 + numSNPs)))
  M = nrow(phenotypes)
  cSums = unlist(colSums(phenotypes))
  phenotypes = phenotypes %>%
    select(which(cSums >= 0.01 * M))
  print(paste("There are", ncol(phenotypes), "interesting phenotypes left"))
  posMat = phenotypes %>%
    filter(curSNP == 1)
  negMat = phenotypes %>%
    filter(curSNP == 0)
  reduced = removeDominated(posMat, negMat)
  if (length(reduced[[3]]) < ncol(phenotypes)) {
    print(paste("Reduction was effective for SNP number", SNPpos))
  }
  return(ncol(phenotypes) - length(reduced[[3]]))
}

### This function prepares a multiplicity specification for the genotypes and phenotypes based on a reduction record
prepareMultiplicitySpecification = function(genotypes, origPhenotypes, out) {
  constantTrueRows <- which(origPhenotypes[out$rowStatuses == CONSTANT, 1] == 1)
  offsetS  <- length(constantTrueRows)
  offsetsT <- colSums(genotypes[constantTrueRows, , drop = FALSE], na.rm = TRUE)
  posRowStatuses <- out$rowStatuses[out$rowStatuses != CONSTANT]
  rowMultsS <- table(abs(posRowStatuses))
  stopifnot(all(names(rowMultsS) == as.character(seq_along(rowMultsS))))
  ### NOTE: these lines implicitly use the fact that the index for CONSTANT is 0 and duplicates are negative
  redGenotypes <- genotypes[out$rowStatuses != CONSTANT, , drop = FALSE]
  rowMultsT <- do.call(cbind, lapply(split(redGenotypes, abs(posRowStatuses)), function(x) {
    colSums(matrix(x, ncol = ncol(genotypes)), na.rm = TRUE)
  }))
  extraDefinitions <- rbind(c(offsetS, rowMultsS), cbind(offsetsT, rowMultsT))
  rownames(extraDefinitions) <- c("S", paste0("T", seq_along(offsetsT)))
  extraDefinitions
}

### This function creates and solves an ILP (integer linear program) for maximizing an objective function over input phenotypes
### The arguments are similar to the function below, except that the genotype is not supplied and phenotypes must be a matrix
### The objVector specifies the vector of objective function values to apply to each entry of the phenotype during optimization
### The boundValue parameter specifies a bound on the objective function; not achieving this bound makes the problem infeasible
createAndSolveILPOld = function(phenotypes, objVector, type = "CNF", K = 3, L = 3, filename = "Test.lp", boundValue = NA, 
                                startSol = NULL, extraConstraints = NULL) {
  n <- nrow(phenotypes)
  p <- ncol(phenotypes)
  goodEntries <- which(phenotypes == ifelse(type == "CNF", 1, 0), arr.ind = TRUE)
  goodEntries %<>%
    as_tibble %>%
    arrange(row, col)
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
  coln <- c(paste0("U", as.vector(outer(1:K, 1:p, pasteI))), paste0("P", as.vector(outer(1:K, 1:n, pasteI))), paste0("p", 1:n))
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
      add(numExtraVar)
    numExtraConst <- 2 * numSeg + 6
    numConst %<>%
      add(numExtraConst)
    numExtraNonZeros <- 2 * numSTCoeffs + 8 * numSeg + 4
    numNonZeros %<>%
      add(numExtraNonZeros)
    coln %<>%
      c("Sum", "Total", "Z", paste0(rep("Z",numSeg), seq_len(numSeg)), paste0(rep("S",numSeg), seq_len(numSeg)))
  }
  fullObjVector <- rep(0, numVar)
  fullObjVector[pInds] <- objVector
  rowInds <- rep(NA, numNonZeros)
  colInds <- rep(NA, numNonZeros)
  values  <- rep(NA, numNonZeros)
  pos <- 0
  lhs <- rep(0, numConst)
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
      lhs[curRows]  <- 1
      lhs[curRowsN] <- -1
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
    lhs[lastOrRow + n * K + (1:n)] <- K - 1
  }
  if (type == "DNF") {
    values[(3 * N + n) * K + (1:(3 * Sz + n))] %<>%
      multiply_by(-1)
  }
  lastRange <- numNonZeros - ifelse(!is.null(extraConstraints), numExtraNonZeros, 0) - (2 * p * K) + 1:(p * K)
  rowInds[lastRange] <- rep(LRows, p)
  colInds[lastRange] <- as.vector(UInds)
  values [lastRange] <- rep(1, p * K)
  lhs[lastAndRow + (1:K)] <- L
  finalRange <- numNonZeros - ifelse(!is.null(extraConstraints), numExtraNonZeros, 0) - (p * K) + 1:(p * K)
  rowInds[finalRange] <- rep(FRows, p)
  colInds[finalRange] <- as.vector(UInds)
  values [finalRange] <- rep(-1, p * K)
  lhs[lastAndRow + K + (1:K)] <- -1
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
                             as.vector(matrix(rbind(indSi, indZi), ncol = 2)), ### this "interlaces" S and Z indices like zip
                             indT, indZi, indSi, ### x-value constraint
                             indZ, indZi, indSi, ### y-value constraint
                             indSi)
    Xvalues <- extraConstraints[[1]] %>% 
      pull(1)
    deltaXs <- diff(Xvalues)
    Yvalues <- extraConstraints[[1]] %>% 
      pull(2)
    deltaYs <- diff(Yvalues)
    values [extraRange] <- c(-extraConstraints[[2]][1, -1], 1, 
                             -extraConstraints[[2]][2, -1], 1,
                             1, -1,
                             rep(1, numSeg),
                             rep(c(1, -1), numSeg),
                             c(-1, Xvalues %>% head(-1), deltaXs),
                             c(-1, Yvalues %>% head(-1), deltaYs),
                             rep(-1, numSeg))
    lhs[lastAndRow + 2 * K + (1:4)] <- c(extraConstraints[[2]][1, 1], extraConstraints[[2]][2, 1], 0, as.double(numSeg > 0))
    if (numSeg == 0) {
      lhs[lastAndRow + 2 * K + (5:6)] <- c(-Xvalues[1], -Yvalues[1])
    }
  }
  model <- glpkAPI::initProbGLPK()
  glpkAPI::setProbNameGLPK(model, paste(type, "for maximum agreement or covariance with genotype"))
  glpkAPI::setObjDirGLPK(model, glpkAPI::GLP_MIN)
  glpkAPI::addColsGLPK(model, ncols = numVar)
  glpkAPI::setColsBndsObjCoefsGLPK(model, j = seq_len(numVar), lb = NULL, ub = NULL, obj_coef = fullObjVector,
                                   type = rep(glpkAPI::GLP_FR, numVar))
  colTypeVector <- rep(glpkAPI::GLP_BV, numVar)
  if (!is.null(extraConstraints)) {
    colTypeVector[numVar - numExtraVar + (1:numExtraVar)] <- 
      c(rep(glpkAPI::GLP_IV, 2), rep(glpkAPI::GLP_CV, 1), rep(glpkAPI::GLP_BV, numSeg), rep(glpkAPI::GLP_CV, numSeg))
  }
  glpkAPI::setColsKindGLPK(model, j = seq_len(numVar), kind = colTypeVector)
  glpkAPI::addRowsGLPK(model, nrows = numConst)
  glpkAPI::loadMatrixGLPK(model, ne = numNonZeros, ia = rowInds, ja = colInds, ra = values)
  rowTypeVector = rep(glpkAPI::GLP_UP, numConst)
  if (!is.null(extraConstraints)) {
    rowTypeVector[numConst - numExtraConst + c(1, 2, 4, numSeg + 5, numSeg + 6)] <- glpkAPI::GLP_FX
  }
  lowerBounds <- rep(0, numConst)
  lowerBounds[rowTypeVector == glpkAPI::GLP_FX] <- lhs[rowTypeVector == glpkAPI::GLP_FX]
  glpkAPI::setRowsBndsGLPK(model, i = seq_len(numConst), lb = lowerBounds, ub = lhs, type = rowTypeVector)
  glpkAPI::setColsNamesGLPK(model, j = seq_len(numVar), cnames = coln)
  glpkAPI::writeLPGLPK(model, fname = filename)
  if (!is.null(startSol)) {
    solFname <- makeInitialSolution(startSol, filename, type, K, phenotypes)
  } else {
    solFname <- NULL
  }
  # note: the objective function is negated during the processing step
  output <- runAndParse(filename, boundValue = boundValue, startSol = solFname)
  if (!all(is.na(output$usedVars))) {
    rownames(output$usedVars) <- colnames(phenotypes)[seq_len(nrow(output$usedVars))]
  }
  output
}

### Slow version of the program currently in Randomizer.R - the new version works directly on the matrix via fast index searches
### This function produces a random member of the set of binary matrices with the same size, row and column sums as the given one
### If keepRownames = TRUE, the original row    names are used; otherwise, R1, R2, ..., Rm are the row    names where m = # rows
### If keepColnames = TRUE, the original column names are used; otherwise, C1, C2, ..., Cn are the column names where n = # cols
randomizeBinaryMatrixOld = function(binaryMatrix, keepRownames = FALSE, keepColnames = FALSE) {
  m <- nrow(binaryMatrix)
  n <- ncol(binaryMatrix)
  nnz <- sum(binaryMatrix)
  print(paste("Processing a", m, "by", n, "input matrix with", nnz, "nonzeros"))
  initGraph <- graph.empty(m + n, directed = TRUE)
  V(initGraph)$name <- c(paste0("R", 1:m), paste0("C", 1:n))
  initGraph %<>%
    graph.union(makeRowColumnGraph(matrix(as.logical(binaryMatrix), ncol = ncol(binaryMatrix)), directed = TRUE), byname = TRUE)
  randGraph <- initGraph
  iter <- 1
  M <- ecount(initGraph)
  maxIter <- KAPPA * M
  print(paste("There are", maxIter, "edge switches to perform"))
  while (iter < maxIter) {
    curInds  <- sample(E(randGraph), 2)
    curEdges <- ends(randGraph, curInds, names = FALSE)
    A <- curEdges[1,1]
    B <- curEdges[1,2]
    C <- curEdges[2,1]
    D <- curEdges[2,2]
    if (A == C || B == D || are_adjacent(randGraph, A, D) || are_adjacent(randGraph, C, B)) {
      next
    }
    randGraph <- delete_edges(randGraph, curInds)
    randGraph <- add_edges(randGraph, edges = c(A, D, C, B))
    iter <- iter + 1
    if (iter %% 100 == 0) {
      print(iter)
    }
  }
  outMatrix <- get.adjacency(randGraph, type = "both", names = FALSE, sparse = TRUE)
  outMatrix <- outMatrix[1:m, m + (1:n), drop = FALSE]
  outMatrix <- as.matrix(outMatrix)
  if (keepRownames) {
    rownames(outMatrix) <- rownames(binaryMatrix)
  }
  if (keepColnames) {
    colnames(outMatrix) <- colnames(binaryMatrix)
  }
  outMatrix
}

### From ForDmitri.R - a bad formulation for the TSP (way too slow!)
solveBestOrderILP = function(coeffMat, filename = "BestOrder.lp") {
  N <- nrow(coeffMat)
  stopifnot(ncol(coeffMat) == N)
  numConst <- N^3 + 2*N
  Lines <- rep("", numConst + 6)
  Lines[1] <- "Minimize"
  Lines[2] <- paste0(coeffMat, " F", row(coeffMat), "C", col(coeffMat), collapse = " + ")
  Lines[3] <- "Subject to"
  pos <- 3
  KRange <- 1:(N - 1)
  nextK <- KRange + 1
  for (i1 in 1:N) {
    for (i2 in 1:N) {
      range <- pos + KRange
      Lines[range] <- paste0("X", i1, "C", KRange, " + X", i2, "C", nextK, " - F", i1, "C", i2, " <= 1")
      pos <- tail(range, 1)
    }
  }
  for (i in 1:N) {
    Lines[pos + i]     <- paste0(paste0("X", i, "C", 1:N, collapse = " + "), " = 1")
    Lines[pos + N + i] <- paste0(paste0("X", 1:N, "C", i, collapse = " + "), " = 1")
  }
  pos <- pos + 2 * N + 1
  Lines[pos] <- "X1C1 = 1"
  Lines[pos + 1] <- "Binary"
  Lines[pos + 1 + (1:N^2)] <- as.vector(outer(1:N, 1:N, function(u, v) { paste0("X", u, "C", v) }))
  Lines[length(Lines)] <- "End"
  f = file(filename)
  writeLines(Lines, f)
  close(f)
}

### This obsolete function prepares a starting solution using the topN best phenotypes from each end and other input parameters.
prepareStartingSolution = function(topN, objValues, ind,  P, curObjVector, type, K, L, filename) {
  curOrder <- order(objValues[ind,])
  highestN <- head(curOrder, topN)
  lowestN  <- tail(curOrder, topN)
  redPheno <- P[, sort(unique(c(lowestN, highestN)))]
  print(paste("The original phenotype matrix is", nrow(P), "by", ncol(P)))
  print(paste("The reduced phenotype matrix is", nrow(redPheno), "by", ncol(redPheno)))
  if (ncol(redPheno) < ncol(P)) {
    microP   <- reducePhenotypesAndObjectives(redPheno, curObjVector, columnsOnly = (type == "LCP"))
    startSol <- createAndSolveProblem(microP$phenotype, microP$objVectors, type = type, K = K, L = L, filename = filename)
  } else {
    startSol <- NULL
  }
  startSol
}

### Alternative way to generate the coefficient matrix (full version)
alternativeCreateAndSolveILP = function(breakSymmetry = TRUE, byClause = FALSE, extraConstraints = NULL,
                                        TPscore = 1, TNscore = 1, FPscore = -1, FNscore = -1, topN = TOP_N, smallerK = FALSE) {
  rown <- c(paste0("OrLower", as.vector(outer(1:N, 1:K, pasteI))), paste0("OrUpper", as.vector(outer(1:n, 1:K, pasteI))))
  rown <- c(rown, paste0("AndLower", c(as.vector(outer(1:n, 1:K, pasteI)))), paste0("AndUpper", 1:n), "TotalSum", "PartialSum")
  coeffMat <- matrix(0, numConst, numVar, dimnames = list(rown, coln))
  coeffMat[cbind(curRows, curIndsP)] <- 1
  coeffMat[cbind(curRows, curIndsN)] <- -1
  coeffMat[cbind(curRowsNext, curIndsN)] <- 1
  coeffMat[cbind(curRowsNext, curIndsP)] <- -1
  ### If breakSymmetry = TRUE additional constraints are added to the problem in order to simplify it (avoid symmetry of clauses)
  ### If byClause = FALSE then symmetry is broken based on the partial phenotypes instead; it involves more variables for large n
  if (breakSymmetry) {
    Dim <- ifelse(byClause, p, n)
  }
  numConst <- K * (N + 2 * n + 2) + n + ifelse(breakSymmetry, Dim * (K - 1), 0)
  if (breakSymmetry) {
    BRows <- lastAndRow + 2 * K + (1:(Dim * (K - 1)))
  }
  numNonZeros <- K * (3 * N + 4 * n + 2 * p) + n + ifelse(breakSymmetry, (K - 1) * (Dim * (Dim + 3)/2), 0)
  if (!is.null(extraConstraints)) {
    extraRows <- lastAndRow + 2 * K + ifelse(breakSymmetry, Dim * (K - 1), 0) + (1:(numExtraConst + 2))
  }
  endRange <- (pos + 1):(numNonZeros - 2 * p * K - ifelse(breakSymmetry, (K - 1) * (Dim * (Dim + 3)/2), 0))
  lastRange <- (numNonZeros - 2 * p * K - ifelse(breakSymmetry, (K - 1) * (Dim * (Dim + 3)/2), 0)) + (1:(p * K))
  finalRange <- numNonZeros - (p * K) - ifelse(breakSymmetry, (K - 1) * (Dim * (Dim + 3)/2), 0) + (1:(p * K))
  if (breakSymmetry) {
    symRange <- numNonZeros - (K - 1) * (Dim * (Dim + 3)/2) + (1:((K - 1) * (Dim * (Dim + 3)/2)))
    rowInds[symRange] <- c(BRows                , rep(BRows, rep(1:Dim, each = K - 1)))
    if (byClause) {
      colInds[symRange] <- as.vector(c(UInds[-1,] , do.call(c, lapply(1:Dim, function(x) {t(UInds[-K, 1:x])}))))
    } else {
      colInds[symRange] <- as.vector(c(PInds[-1,] , do.call(c, lapply(1:Dim, function(x) {t(PInds[-K, 1:x])}))))
    }
    values [symRange] <- c(rep(1, (K - 1) * Dim), rep(-1, (K - 1) * Dim * (Dim + 1)/2))
  }
}  
### Alternative that is no longer used: specifying priorties (see below) 
alternativeCreateAndSolveILPs = function(breakSymmetry = FALSE) {
  ### If ord = TRUE, then there is a specification for the solver to branch on the usage variables; priority given by the index k
  if (ord) {
    ordFN <- "Test.ord"
    makeOrdFile(ordFN, ncol(phenotypes), K)
  } else {
    ordFN <- NULL
  }
  ### was used previously, but not anymore
  usedRows <- match(seq_len(max(out$rowStatuses)), out$rowStatuses)
  if (!is.null(topN) && ncol(phenotypes) > 2 * topN) {
    valueFunction <- ifelse(type == "LCP", computeCorr, ifelse(objective == "covariance", computeCovariances, computeAgreements))
    objValues <- valueFunction(origGenotypes, origPhenotypes) %>%
      magrittr::extract(, usedCols, drop = FALSE)
  }
  ### was used previously, but has been replaced by a PWL constraint
  if (!is.null(maxPValue)) {
    curDefinitions <- curDefinitions[, c(1, 1 + goodInds), drop = FALSE]
    curConstraints <- constrainPValues(N = ns[ind], nPlus = sumYs[ind], pMax = maxPValue, convex = TRUE)
    curExtras <- list(C = curConstraints, D = curDefinitions) ### TODO: Consider adding them to the warm start problem as well
  } else {
    curExtras <- NULL
  }
  ### obsolete solutions
  if (!is.null(topN)) { ### construct a warm start by using N phenotypes with the highest and lowest objective w.r.t. genotype
    startSol <- prepareStartingSolution(topN, singleRes$allValues, ind, P, curObjVector, type, K, L, 
                                        filename = paste0("Micro", fn))
  } else {
    if (smallerK) { ### construct a warm start by using a lower value of K
      startSol <- createAndSolveProblem(P, curObjVector, type = type, K = K - 1, L = L, filename = paste0("Mini", fn), 
                                        boundValue = ifelse(type == "LCP", NA, boundValues[ind]), startSol = NULL)
    } else {
      startSol <- NULL
    }
  }
  ### Excerpt from a past post-processing step
  if (type == "LCP") {
    myTab %<>%
      mutate_at("opt", ~divide_by(., sqrt(sumYs * (ns - sumYs) / ns)))
  }
}

### No longer used; function to create an order on the choice variables
makeOrdFile = function(filename, numCol, numIter, prefix1 = PREFIX1, prefix2 = PREFIX2) {
  numVars <- numCol * numIter
  Lines <- rep("", numVars + 2)
  Lines[1] <- "NAME"
  allLines <- outer(1:numCol, 1:numIter, function(x, y) { paste0(" ", " ", " ", " ", prefix1, x, prefix2, y, " ", y) } )
  Lines[2:(numVars + 1)] <- as.vector(allLines)
  Lines[length(Lines)] <- "ENDATA"
  f = file(filename, "w")
  writeLines(Lines, f)
  close(f)
}

### No longer used due to incorrectness for the application
findUpperEnvelope = function(Z0) {
  Z <- Z0 %>% 
    as_tibble() %>%
    group_by(Total) %>%
    summarize(Sum = max(Sum)) %>%
    ungroup %>%
    as.matrix
  N <- nrow(Z)
  ind <- 1
  envelope <- c(ind)
  while (ind < N) {
    j <- ind + 1
    while (testCorrectPiece(Z, ind, j)) {
      j %<>%
        add(1)
      if (j > N) {
        break
      }
    }
    j %<>%
      subtract(1)
    envelope %<>%
      c(j)
    ind <- j
  }
  output <- Z[envelope,]
  output
}

### Obsolete functions from Simulation.R

### This function carries out a random (seeded) simulation on a specified number of patients and phenotypes with a known formula
### It is useful for testing the performance of the logical code in situations where an approximately correct solution is known. 
conductSimulation = function(Q = 4, p = 10, n = 10^Q, K = 3, L = 3, TYPE = "CNF", numPerturb = 0, FULL = FALSE) {
  set.seed(MY_SEED)
  innerF <- ifelse(TYPE == "CNF", or, and)
  outerF <- ifelse(TYPE == "CNF", and, or)
  phenotypes <- matrix(as.logical(round(runif(n * p))), nrow = n, dimnames = list(paste0("Patient", 1:n), paste0("p", 1:p)))
  if (FULL) {
    genotype <- outerF(outerF(innerF(innerF(phenotypes[,1], phenotypes[,2]), phenotypes[,3]), 
                              innerF(innerF(phenotypes[,4], phenotypes[,5]), phenotypes[,6])),
                       innerF(innerF(phenotypes[,7], phenotypes[,8]), phenotypes[,9]))
  } else {
    if (p >= 9) {
      genotype <- outerF(outerF(innerF(phenotypes[,1], phenotypes[,2]), innerF(phenotypes[,4], phenotypes[,5])),
                         innerF(phenotypes[,3], phenotypes[,9]))
    } else {
      genotype <- outerF(outerF(innerF(phenotypes[,1], phenotypes[,2]), innerF(phenotypes[,4], phenotypes[,5])), (phenotypes[,3]))
    }
  }
  genotype[sample(1:n, n / 10, replace = FALSE)] <- NA # add 10% random NAs to the genotype 
  perturbInds <- sample(1:n, numPerturb, replace = FALSE)
  genotype[perturbInds] <- !(genotype[perturbInds])
  curFilename <- paste0("TestProblem", TYPE, Q, ".lp")
  objVector <- rep(0, length(genotype))
  objVector[genotype == TRUE]  <- -2
  objVector[genotype == FALSE] <- +2
  nPhenotypes <- matrix(!phenotypes, nrow(phenotypes), dimnames = list(rownames(phenotypes), paste0("NOT_", colnames(phenotypes))))
  fullPhenotypes <- cbind(phenotypes, nPhenotypes)
  solution <- createAndSolveILP(fullPhenotypes, objVector, type = TYPE, K = K, L = L, filename = curFilename)
  if (TRUE) {return(solution)}
  foundGenotype <- postprocessSolution(fullPhenotypes, solution$usedVars, type = TYPE, complement = -1)
  numEqual <- sum(foundGenotype$solution == genotype, na.rm = TRUE)
  print(paste("There are", numEqual, "matching entries"))
  negFilename <- paste0("TestProblemNegated", TYPE, Q, ".lp")
  negSolution <- createAndSolveILP(fullPhenotypes, -objVector, type = TYPE, K = K, L = L, filename = negFilename)
  foundNegGenotype <- postprocessSolution(fullPhenotypes, negSolution$usedVars, type = TYPE, complement = -1)
  numNegEqual <- sum(foundNegGenotype$solution == !genotype, na.rm = TRUE)
  print(paste("There are", numNegEqual, "matching entries after negation"))
  testTab <- cbind(SNP1 = genotype, phenotypes)
  fullTab <- as_tibble(testTab) %>%
    mutate_all(as.numeric) %>%
    bind_cols(tibble(ID = paste0("Patient", 1:n))) %>% 
    select(ID, everything())
  write_tsv(fullTab, paste0("TestProblem", TYPE, Q, "K", K, "L", L, ".tsv"))
  extraTab <- cbind(SNP1 = genotype, fullPhenotypes)
  extraFullTab <- as_tibble(extraTab) %>%
    mutate_all(as.numeric) %>%
    bind_cols(tibble(ID = paste0("Patient", 1:n))) %>% 
    select(ID, everything())
  write_tsv(extraFullTab, paste0("TestProblem", TYPE, Q, "K", K, "L", L, "Extended.tsv"))
  fullTabLogical <- fullTab %>%
    mutate_if(is.numeric, as.logical)
  EV <- log(MAX_P)
  altSolution <- createAndSolveILPs(fullTabLogical, type = TYPE, K = K, L = L, complement = 2, extremeValue = EV, index = 1)
  XTYPE <- ifelse(TYPE == "CNF", "DNF", "CNF")
  crossSolution <- createAndSolveILPs(fullTabLogical, type = XTYPE, K = K, L = L, complement = 2, extremeValue = EV, index = 2)
  output <- list(solution, negSolution, altSolution, crossSolution, genotype, foundGenotype, foundNegGenotype, fullPhenotypes)
  output
}

### This function produces an input file for a multi-SNP test of the brute-force enumeration and the CPLEX-based pipelines.
### p, n and N are the numbers of phenotypes, patients and SNPs, respectively; K is the number of clauses; L is the clause size.
### TYPE is CNF or DNF, according to the formula to generate; fPerturb is the fractional perturbation to use (= numPerturb/n).
### numPerturb is the number of perturbations per combination; fRand is the fraction of random SNPs (out of N).
simulateForDmitri = function(p = 10, n = 10^3, N = 10^3, K = 3, L = 3, TYPE = "CNF", fPerturb = 0.01, numPerturb = 10, fRand = 0.5, sparsity = 0.2) {
  set.seed(MY_SEED)
  innerF <- ifelse(TYPE == "CNF", or, and)
  outerF <- ifelse(TYPE == "CNF", and, or)
  patientNames <- paste0("P", 0:(n - 1))
  phenoNames <- paste0("p", 0:(p - 1))
  phenotypes <- matrix(as.logical(runif(n * p) >= sparsity), nrow = n, dimnames = list(patientNames, phenoNames))
  numRandom <- round(N * fRand)
  numNonRandom <- N - numRandom
  numCombos <- round(numNonRandom/numPerturb)
  genotypes <- matrix(as.logical(runif(n * N) >= 0.5), nrow = n, dimnames = list(patientNames, c()))
  allCombos <- combn(1:p, L)
  pos <- 1
  perturbSize <- round(fPerturb * n)
  records <- vector("list", N)
  for (ind in 1:numCombos) {
    curClauses <- allCombos[, sample(1:ncol(allCombos), K, replace = FALSE), drop = FALSE]
    initGenotype <- rep(ifelse(TYPE == "CNF", TRUE, FALSE), n)
    for (index in 1:K) {
      miniGenotype <- rep(ifelse(TYPE == "CNF", FALSE, TRUE), n)
      curCol <- curClauses[, index]
      for (index2 in 1:L) {
        miniGenotype <- innerF(phenotypes[, curCol[index2]], miniGenotype)
      }
      initGenotype <- outerF(initGenotype, miniGenotype)
    }
    for (count in 1:numPerturb) {
      finalGenotype <- initGenotype
      perturbInds <- sample(1:n, perturbSize, replace = FALSE)
      finalGenotype[perturbInds] <- !(finalGenotype[perturbInds])
      genotypes[, pos] <- finalGenotype
      records[[pos]] <- curClauses
      pos <- pos + 1
    }
  }
  randomOrder <- sample(1:N)
  records <- records[randomOrder]
  genotypes <- genotypes[, randomOrder]
  SNPNames <- paste0("SNP", 0:(N - 1))
  colnames(genotypes) <- SNPNames
  names(records) <- SNPNames
  testTab <- cbind(genotypes, phenotypes)
  fullTab <- as_tibble(testTab) %>%
    mutate_all(as.numeric) %>%
    bind_cols(tibble(ID = paste0("Patient", 1:n))) %>% 
    select(ID, everything())
  write_tsv(fullTab, paste0("TestProblem", TYPE, "n", n, "p", p, "N", N, "K", K, "L", L, ".tsv"))
  save(records, file = paste0("TestSolutions", TYPE, "n", n, "p", p, "N", N, "K", K, "L", L, ".RData"))
  output <- list(fullTab, records)
  output
}

testResults = function(resFile = "ForDmitri/TestProblemCNFn1000p10N1000K3L3.scores_3x3.csv") {
  recordFile <- str_replace(str_replace(resFile, "Problem", "Solutions"), ".scores_3x3.csv", ".RData")
  testRes <- read_csv(resFile)
  load(recordFile)
  allNames <- names(records)
  for (Name in allNames) {
    curRes <- testRes %>%
      filter(SNP == Name) %>%
      slice(1) %>%
      select(c("status", "formula"))
    curCombo <- curRes %>%
      pull(formula) %>%
      str_split(" AND ") %>%
      unlist %>%
      str_remove("\\(") %>%
      str_remove("\\)") %>%
      str_remove_all("OR") %>%
      str_remove_all("p") %>%
      str_split("[ ]+") %>%
      lapply(function(x) {as.integer(x) + 1})
    curRec <- records[[Name]]
    if (!is.null(curRec)) {
      curTest <- curRec %>%
        t() %>%
        as.data.frame() %>%
        magrittr:extract(do.call(order, .), )
      for (ind in 1:nrow(curTest)) {
        stopifnot(sort(unique(curTest[ind,])) == curCombo[[ind]])
      }
    }
  }
  return(TRUE)
}

### This function creates symmetry-breaking constraints for (K,L)-[C/D]NF; if an LP file specified in the input, adds them to it
createSBConstraints = function(K, p, inputFile = NULL) {
  fVars <- paste0("F", as.vector(outer(1:K, 1:p, pasteI)))
  uVars <- paste0("U", as.vector(outer(1:K, 1:p, pasteI)))
  indsMat <- matrix(1:(K * p), nrow = K)
  ConstEq <- paste(fVars[indsMat[,1]], "-", uVars[indsMat[,1]], " = 0")
  ConstIneq1 <- paste(fVars[indsMat[,-1]], "-", uVars[indsMat[,-1]], " >= 0")
  ConstIneq2 <- paste(fVars[indsMat[,-1]], "-", fVars[indsMat[,-p]], ">= 0")
  ConstIneq3 <- paste(uVars[indsMat[,-1]], "+", fVars[indsMat[,-p]], "-", fVars[indsMat[,-1]], ">= 0")
  ConstIneq4 <- paste(fVars[indsMat[-1,]], "-", fVars[indsMat[-K,]], " <= 0")
  allConst <- c(ConstEq, ConstIneq1, ConstIneq2, ConstIneq3, ConstIneq4)
  allBounds <- paste("0 <= ", fVars, " <= 1")
  if (!is.null(inputFile)) {
    initLines <- readLines(inputFile)
    startConst <- which(initLines %in% c("Subject To", "Subject to"))
    startBounds <- which(initLines == "Bounds")
    startGenerals <- which(initLines == "Generals")
    stopifnot(startConst < startBounds)
    stopifnot(startBounds < startGenerals)
    finalLines <- c(initLines[1:startConst], allConst, initLines[(startConst + 1):startBounds], allBounds, 
                    initLines[(startBounds + 1):startGenerals], fVars, initLines[(startGenerals + 1):length(initLines)])
    outputFile <- str_replace(inputFile, ".lp", "Augmented.lp")
    fw <- file(outputFile, "w")
    writeLines(finalLines, fw)
    close(fw)
  }
  output <- list(generals = fVars, constraints = allConst, bounds = allBounds)
  output
}

### Obsolete functions from Postprocess.R
simplifyFormulas = function(inputFile = "UKBB_Diabetes_300kSNPs_feasible_discovery_all_chr_v3_june_full_sweep.csv") {
  Tab = read_csv(inputFile)
  stopifnot("formula" %in% colnames(Tab))
  N = nrow(Tab)
  Tab = Tab %>%
    mutate(simpleFormula = map_chr(formula, ~{simplifyFormula(., type = "CNF")}), .after = formula)
  write_csv(Tab, str_replace(inputFile, ".csv", "_extended.csv"))
  Tab
}

### Obsolete functions from Postprocess.R (now renamed to Preprocess.R) - used primarily for versions of the Imai-Iri algorithm
getPoint = function(allPoints, sign, index) {
  point <- allPoints[[sign]] %>%
    slice(index) %>%
    select(1:2)
  point
}

getPrevPoint = function(allPoints, sign, index) {
  prevPoint <- allPoints[[sign]] %>% 
    slice(index) %>%
    select(PrevX, PrevY)
  prevPoint
}

getNextPoint = function(allPoints, sign, index) {
  nextPoint <- allPoints[[sign]] %>% 
    slice(index) %>%
    select(NextX, NextY)
  nextPoint
}

setNextPoint = function(allPoints, sign, index, point) {
  allPoints[[sign]][index, c("NextX", "NextY")] <- point
  allPoints
}

setPoint = function(allPoints, sign, index, point) {
  allPoints[[sign]][index1, 1:2] <- point
  allPoints
}

setPrevPoint = function(allPoints, sign, index, point) {
  allPoints[[sign]][index, c("PrevX", "PrevY")] <- point
  allPoints
}

setNextPoint = function(allPoints, sign, index, point) {
  allPoints[[sign]][index, c("NextX", "NextY")] <- point
  allPoints
}

### Check to see how many extra points the line segment between points i and j includes relative to the specified point set Z
### Assumes that the point set has two columns; "Sum", which can vary from 1 to nrow(Z), and "Total", which has consecutive 
### If below = FALSE assume that the segment lies above the point set and computes the number of extra integer points below it
### If below = TRUE, assume that the segment lies below the point set and computes the number of extra integer points above it
testCorrectPiece = function(Z, startPoint, endPoint, below = TRUE) {
  i <- which(Z[,"Total"] == startPoint["Total"])
  j <- which(Z[,"Total"] == endPoint  ["Total"])
  if (j - i == 1) {
    return(0)
  }
  roundFunction <- ifelse(below, floor, ceiling)
  curDelta <- endPoint - startPoint
  curLength <- j - i - 1
  curSlope <- curDelta["Sum"]/curDelta["Total"]
  curPreds <- roundFunction((1:curLength) * curSlope) + startPoint["Sum"]
  curTruth <- Z[(i + 1):(j - 1), "Sum"]
  return(sum(curPreds - curTruth))
}

getDirections = function(convexHull, pointSet, tol = TOL) {
  convexHull %<>%
    apply(1, sort) %>%
    t
  M <- nrow(convexHull)
  N <- nrow(pointSet)
  dirs <- rep(NA, M)
  for (ind in 1:M) {
    curRow <- sort(convexHull[ind,])
    curStart <- curRow[1]
    curEnd  <- curRow[2]
    curInds <- setdiff(1:N, c(curStart, curEnd))
    iter <- 1
    while (iter <= N - 2) {
      curInd <- curInds[iter]
      turnAngle <- getAngle(pointSet[c(curStart, curInd, curEnd), ])
      if (abs(turnAngle) <= tol) {
        iter <- iter + 1
      } else if (turnAngle > tol) {
        dirs[ind] <- BELOW
        break
      } else { ### turnAngle < -tol
        dirs[ind] <- ABOVE
        break
      }
    }
  }
  convexHull %<>%
    set_colnames(c("start", "end")) %>%
    as_tibble %>%
    bind_cols(dir = dirs)
  convexHull
}

prepareConstraints = function(hull, below = TRUE) {
  if (nrow(hull) == 1) { ### special case: single extreme point
    return(tibble(Sum = c(1L, -1L, 0L), Total = c(0L, 0L, -1L), lhs = c(hull$Sum[1], -hull$Sum[1], -hull$Total[1])))
  }
  Q <- hull %>%
    as.matrix
  deltaX <- diff(Q[,1])
  deltaY <- diff(Q[,2])
  xCoeffs <- deltaY
  yCoeffs <- -deltaX
  uBounds <- deltaY * Q[-1, 1] - deltaX * Q[-1, 2]  
  constMat <- cbind(xCoeffs, yCoeffs, uBounds)
  colnames(constMat) <- c(colnames(hull), "lhs")
  constMat %<>%
    as_tibble %>%
    mutate_all(as.integer)
  ### Previously: simplified the line equations via:
  # constMat %<>%
  #   as_tibble %>%
  #   mutate(gcd = pracma::gcd(Sum, Total)) %>%
  #   mutate_all(as.integer) %>%
  #   mutate_all(~divide_by_int(., gcd)) %>%
  #   select(-gcd)
  constMat
}

constrainPValues = function(N, nPlus, pMax = 0.05, convex = TRUE, below = TRUE) {
  Z <- computePValueBounds(N = N, nPlus = nPlus, pMax = pMax, logSpace = TRUE, fast = TRUE)
  if (convex) {
    hull <- getConvexHull(Z, checkTightness = FALSE, below = below)
    constraints <- prepareConstraints(hull, below = below)
    print(paste("Prepared", nrow(constraints), "additional constraints to be added"))
    return(constraints)
  } else {
    print("Warning: no alternatives to convex hulls are currently available!")
  }
  return(Z)
}

getConvexHull = function(Z, checkTightness = TRUE, below = TRUE) {
  if (nrow(Z) == 1) {
    return(Z)
  }
  if (nrow(Z) == max(Z$Total) || nrow(Z) == 2) {
    hull <- Z %>%
      slice(1, nrow(Z)) %>%
      as.matrix
    return(hull)
  }
  CH <- convhulln(Z)
  CH %<>%
    getDirections(Z)
  hull <- CH %>%
    filter(dir == ifelse(below, BELOW, ABOVE)) %>%
    select(-dir) %>%
    arrange(start) %>%
    as.matrix
  Z %<>%
    as.matrix
  N <- nrow(hull)
  stopifnot(all(hull[-1, 1] == hull[-N, 2]))
  goodPoints <- as.vector(c(hull[ ,1], hull[N, 2]))
  hull <- Z[goodPoints, ]
  if (checkTightness) {
    selectorFunction <- ifelse(below, max, min)
    Z1 <- Z %>%
      as_tibble %>%
      group_by(Total) %>%
      summarize(Sum = selectorFunction(Sum)) %>%
      ungroup %>%
      select(Sum, Total) %>%
      as.matrix
    totalOK <- sum(Z[,"Sum"]) - sum(Z[,"Total"]) + nrow(Z)
    numExtras <- 0
    for (ind in seq_len(N)) {
      startPoint <- hull[ind, ]
      endPoint   <- hull[ind + 1, ]
      numExtras %<>%
        add(testCorrectPiece(Z1, startPoint, endPoint, below = below))
    }
    if (numExtras > 0) {
      print(paste("Warning:", numExtras, "extra integer points appear alongside", totalOK, "correct ones in the convex hull"))
    }
  }
  hull
}

makePlotForDaniel = function(exponents = 1:6, fracPositive = 0.5, pMax = 0.05, plot = TRUE) {
  size <- length(exponents)
  allNs <- rep(NA, size)
  allPluses <- rep(NA, size)
  inputSizes <- rep(NA, size)
  outputSizes <- rep(NA, size)
  for (ind in 1:size) {
    print(ind)
    curN <- 10^(exponents[ind])
    allNs[ind] <- curN
    curPlus <- fracPositive * curN
    allPluses[ind] <- curPlus 
    curZ <- computePValueBounds(curN, curPlus, pMax)
    curOutput <- callImaiIri(curZ)
    inputSizes[ind] <- curOutput[[2]]
    outputSizes[ind] <- nrow(curOutput[[1]])
  }
  Tab <- tibble(N = allNs, nPlus = allPluses, pMax = pMax, numNonCollinear = inputSizes, numBreakpoints = outputSizes)
  if (plot) {
    G <- ggplot(data = Tab, mapping = aes(x = N, y = numBreakpoints)) + 
      geom_point() +
      scale_x_log10()
    ggsave(filename = "ForDaniel.pdf")
  }
  Tab
}

sanityCheck = function(Nrange = 10:60) {
  for (N in Nrange) { 
    for (n in 1:N) {
      print(c(N,n))
      Z <- computePValueBounds(N, n, 0.05)
      Q <- callImaiIri(Z)
    }
  }
}

checkConvexPointSet = function(pointSetBdry) {
  N <- nrow(pointSetBdry)
  if (N <= 2) {
    return(TRUE)
  }
  pointSetBdry %<>%
    as.matrix
  slopes <- diff(pointSetBdry[,2])/diff(pointSetBdry[,1])
  return(all(diff(slopes) > 0))
}

### NB: assumes without checking that no triplet of consecutive input points is collinear!
ImaiIriRestrictedAlgorithm = function(inputPoints) {
  inputPoints %<>%
    as.matrix
  n <- nrow(inputPoints)
  stopifnot(ncol(inputPoints) == 2)
  if (n <= 2) {
    return(inputPoints)
  }
  allAngles <- rep(NA, n)
  for (ind in 2:(n - 1)) {
    allAngles[ind] <- sign(getAngle(inputPoints[(ind - 1):(ind + 1), ]))
  }
  starts <- rep(0, choose(n, 2))
  starts[1:(n - 1)] <- 1:(n - 1)
  ends   <- rep(0, choose(n, 2))
  ends  [1:(n - 1)] <- 2:n
  pos <- n
  print(paste("There are", n, "points to process"))
  for (ind1 in 1:(n - 2)) {
    if (ind1 %% 500 == 0) {
      print(paste("Processed", ind1, "points so far"))
    }
    curSign <- allAngles[ind1 + 1]
    curMax <- ind1 + 2
    while (curMax < n && allAngles[curMax] == curSign) {
      curMax <- curMax + 1
    }
    firstInd  <- ind1
    secondInd <- ind1 + 1
    thirdInd <- ind1 + 2
    while (thirdInd <= curMax) {
      if (checkEmptyTriangle(inputPoints, firstInd, secondInd, thirdInd, checkDiagonal = TRUE)) {
        starts[pos] <- firstInd
        ends  [pos] <- thirdInd
        pos <- pos + 1
        secondInd <- thirdInd
        thirdInd <- thirdInd + 1
      } else {
        break
      }
    }
  }
  G <- graph_from_edgelist(cbind(starts[1:(pos - 1)], ends[1:(pos - 1)]), directed = TRUE)
  myShortestPath <- shortest_paths(G, from = 1, to = n, output = "vpath")$vpath[[1]]
  boundaryPoints <- inputPoints[as.vector(myShortestPath), , drop = FALSE] %>%
    as_tibble
  boundaryPoints
}

checkEmptyTriangle = function(pointSet, pos1, pos2, pos3, checkDiagonal = TRUE) {
  point1 <- pointSet[pos1,]
  point2 <- pointSet[pos2,]
  point3 <- pointSet[pos3,]
  delta1 <- point2 - point1
  x1 <- delta1[1]
  y1 <- delta1[2]
  delta2 <- point3 - point2
  x2 <- delta2[1]
  y2 <- delta2[2]
  numPointsAB <- pracma::gcd(x1, y1)
  numPointsBC <- pracma::gcd(x2, y2)
  numPointsCA <- pracma::gcd(x1 + x2, y1 + y2)
  if (checkDiagonal && numPointsCA > 1) {
    return(FALSE)
  }
  Area <- abs(x1 * y2 - x2 * y1)/2
  perimeterPoints <- numPointsAB + numPointsBC + numPointsCA
  numIntegerPoints <- Area + 1 - perimeterPoints/2
  return(numIntegerPoints == 0)
}

### Obsolete function for data mining: looking at how often solutions come from top/bottom N
rankResults = function(group = "SelfReported", numPat = 1000, numSNPs = 15, type = "CNF", K = 3, L = 3, objective = "agreement", 
                       ext = ".csv", N = 10) {
  fname <- paste(c("CPLEXInput", group, "EU", numPat, numSNPs, type, paste0("K", K), paste0("L", L), objective), collapse = "_")
  fname <- paste(c(fname, "Results", "Summary.csv"), collapse = "_")
  Tab1 <- read_delim(fname, delim = ",")
  Tab1S <- Tab1 %>% 
    as.matrix %>%
    apply(1, function(x) {
      unlist(str_split(x,",")[[1]])}) %>%
    t
  colnames(Tab1S) <- colnames(Tab1) %>%
    str_split(",") %>%
    extract2(1) %>%
    unlist
  Tab1S %<>%
    as_tibble %>%
    filter(formula != "") %>%
    select(SNP, formula) %>%
    separate(formula, sep = ifelse(type == "CNF", "AND", "OR"), into = paste0("Clause", 1:K))
  for (ind in 1:K) {
    Tab1S %<>%
      separate(paste0("Clause", ind), sep = ifelse(type == "CNF", "OR", "AND"), into = paste0("Clause", ind, "Term", 1:L))
  }
  Tab1S %<>%
    mutate_all(~trimws(.)) %>%
    mutate_all(~str_remove(str_remove(., "^\\("), "\\)$")) %>%
    separate("SNP", sep = "_", into = c("direction", "variant", "letter"), fill = "left") %>%
    replace_na(list(direction = "")) %>%
    unite("SNP", variant:letter)
  inputFile <- paste(c("CPLEXInput", group, "EU", numPat, paste0(numSNPs, str_sub(ext, 2))), collapse = "_")
  readFunction <- ifelse(ext == '.tsv', read_tsv, read_csv)
  inputTab <- readFunction(inputFile, col_types = cols(ID = "c", .default = "l"))
  SNPs <- inputTab %>%
    select(Tab1S$SNP)
  phenotypes <- inputTab %>%
    select(-(1:(numSNPs + 1)))
  valueFunction <- ifelse(objective == "agreement", computeAgreements, computeCovariances)
  values <- valueFunction(as.matrix(SNPs), as.matrix(phenotypes))
  rankValues <- t(apply(values, 1, rank))
  myRanks <- Tab1S %>%
    select(Clause1Term1:paste0("Clause", K, "Term", L))
  for (ind in 1:ncol(myRanks)) {
    curCol <- myRanks[[ind]]
    goodPos <- which(!is.na(curCol))
    curCol[goodPos] <- rankValues[cbind(match(Tab1S$SNP[goodPos], colnames(SNPs)), match(curCol[goodPos], colnames(phenotypes)))]
    myRanks[[ind]] <- curCol
  }
  numTopBottomN <- (myRanks <= N | myRanks >= ncol(phenotypes) - N) %>%
    t %>%
    split(rep(1:K, each = L)) %>%
    lapply(function(x) {rowSums(matrix(x, ncol = L, byrow = TRUE))})
  numTopBottomN <- do.call(cbind, numTopBottomN)
  rownames(numTopBottomN) <- Tab1S$SNP
  output <- list(Tab = Tab1S, ranks = myRanks, topN = numTopBottomN)
  output
}

### Option from computeStats
computeStatsOld = function(Tibble, pvals = TRUE, logSpace = FALSE) {
  if (pvals) {
    ChiSquaredP  <- apply(myMat, 1, function(x) {  chisq.test(matrix(x, 2)                   )$p.value} )
    FisherExactP <- apply(myMat, 1, function(x) { fisher.test(matrix(x, 2), alternative = "g")$p.value} )
    Tibble %<>% mutate(ChiSquaredP = ChiSquaredP, FisherExactP = FisherExactP)
  }
  Tibble
}

### Option from CallImaiIri
callImaiIri = function(extremePoints, below = TRUE, verify = TRUE, width = WIDTH, dryRun = FALSE) {
  if (restricted) {
    breakpoints <- ImaiIriRestrictedAlgorithm(pointsP)
  } else {
    ... 
  }
}

### Plot the sizes of the boundaries obtained by the Imai-Iri algorithm at different widths, with n = 1e6, nPlus = 1e5, p = 5e-8
plotSegmentsVsEpsilon = function(plotFile = "PlotForDaniel.pdf") {
  Tab <- tibble(
    width =       c(0.6  , 0.7  , 0.8  , 0.9 , 0.96, 0.97, 0.98, 0.99, 0.995, 1.00), 
    numSegments = c(4667L, 3029L, 1441L, 884L, 343L, 259L, 193L, 150L, 147L,  121L)
  )
  Tab %<>%
    mutate(epsilon = 1 - width)
  ggplot(data = Tab, mapping = aes(x = epsilon, y = numSegments)) + geom_point()
  ggsave(filename = plotFile)
  Tab
}

testEmptyTriangle = function() {
  pointSet1 <- rbind(c(0, 0), c(2, 2), c(3, 4))
  stopifnot(checkEmptyTriangle(pointSet1, 1, 2, 3, checkDiagonal = TRUE))
  pointSet2 <- rbind(c(0, 0), c(1, 1), c(2, 4))
  stopifnot(checkEmptyTriangle(pointSet2, 1, 2, 3, checkDiagonal = FALSE))
  stopifnot(!checkEmptyTriangle(pointSet2, 1, 2, 3, checkDiagonal = TRUE))
  pointSet3 <- rbind(c(0, 0), c(3, 4), c(4, 6))
  stopifnot(!checkEmptyTriangle(pointSet3, 1, 2, 3, checkDiagonal = TRUE))
}

countIntegral = function(pointSet) {
  pointSet %<>%
    mutate(intSum = near(Sum, round(Sum)), intTotal = near(Total, round(Total)), intBoth = intSum & intTotal)
  print(paste("There are", sum(pointSet$intBoth), "integral points out of", nrow(pointSet)))
}

### This function computes, for a given N and nPlus, the set of critical value pairs (i, j) such that the 2x2 contingency table
### with row sums nPlus and N - nPlus, and column sums i and N - i has p-value below pMax when its top left entry is at least j. 
computePValueBoundsOld = function(N, nPlus, pMax, logSpace = FALSE, fast = TRUE) {
  print(paste("The p-value being considered as the cutoff is", pMax))
  print(paste("Preparing p-value bounds with N =", N, "and nPlus =", nPlus))
  stopifnot(nPlus <= N)
  nMinus <- N - nPlus
  if (nPlus > nMinus) {
    fast <- FALSE ### The conditions for speeding up the computation may not be met
  }
  stopifnot(0 <= pMax && pMax <= 1)
  if (nPlus == 0) {
    return(tibble(Sum = integer(), Total = integer()))
  }
  iRange <- 1:nPlus
  limit <- ifelse(logSpace, log(pMax), pMax)
  criticalValues <- rep(NA, nPlus)
  prevValue <- 0
  for (i in iRange) {
    range <- ifelse(fast, max(prevValue - M_STEP, 0, na.rm = TRUE), 0):ifelse(fast, min(prevValue + M_STEP, i, na.rm = TRUE), i)
    hyperVector <- phyper(range, nPlus, nMinus, i, lower.tail = FALSE, log.p = logSpace)
    prevValue <- range[min(which(hyperVector < limit))]
    criticalValues[i] <- prevValue + 1 ### the reason for the +1 is the offset of phyper with respect to fisher.test with or = 1
  }
  output <- tibble(Sum = iRange, Total = criticalValues) %>%
    filter(Sum >= Total)
  output
}
### End obsolete functions from Postprocess.R (now renamed to Preprocess.R)

### Compute the covariance between two given input matrices; if binary = TRUE, scale and round the result to the nearest integer
computeDotProducts = function(inputMatrix1, inputMatrix2, binary = FALSE) {
  n <- nrow(inputMatrix1)
  stopifnot(nrow(inputMatrix2) == n)
  output <- cov(inputMatrix1, inputMatrix2) * (n - 1)
  if (binary) {
    output %<>%
      multiply_by(n) %>%
      round
  }
  output
}

### From an older version of MIQCP, where there was no QR factorization of the coefficient matrix prepared ahead of time
createAndSolveMIQCPOld = function(qMatrix, objVector, S = 1, K = 3, M = 10, filename = "Test.lp", startSol = NULL, positive = FALSE,
                               boundValue = NA) {
  diagonal <- cbind(1:Dim, 1:Dim)
  upperHalf <- combn2(1:Dim)
  quadFunction1 <- paste0(qMatrix[upperHalf] * 2, " U", upperHalf[,1], " * ", "U", upperHalf[,2], collapse = " +")
  quadFunction2 <- paste0(qMatrix[diagonal], " U", diagonal[,1], "^2", collapse = " +")
  Lines[4] <- paste0("[ ", quadFunction1, " +", quadFunction2, " ] <= ", 1)
  ### obsolete option
  if (positive) {
    Lines[start2 + 1 + (1:Dim)] <- paste0("U", 1:Dim, " >= 0")
  }
}

### FROM Geometry.R

INT_BOUNDARY = FALSE ### If TRUE, makes the boundary consist only of integral points

### This function turns an input 2-dimensional vector into a unit vector in the same direction
normalizeVector = function(vector) {
  len <- pracma::hypot(vector[1], vector[2])
  if (near(len, 0)) {
    return(vector)
  }
  normVector <- vector/sqrt(len)
  normVector
}

### Determines the (sign of the) angle by which we need to turn to get from point 1 to point 3 via point 2 (ignores magnitudes!)
### The return value is positive if the turn is a right turn, negative if the turn is a left turn, 0 if the points are collinear
getAngle = function(threePointSet) {
  stopifnot(all(dim(threePointSet) == c(3, 2)))
  threePointSet %<>%
    as.matrix
  point1 <- threePointSet[1, ]
  point2 <- threePointSet[2, ]
  point3 <- threePointSet[3, ]
  if (comparePoints(point1, point2) || comparePoints(point2, point3)) { ### degenerate situation!
    return(0)
  }
  vector1 <- normalizeVector(point2 - point1)
  vector2 <- normalizeVector(point2 - point3)
  angle <- vector1[1] * vector2[2] - vector2[1] * vector1[2]
  angle
}

### Remove all but the first and the last point of any group of three or more consecutive input points that are collinear.
### Note that this function is idempotent (at least in exact arithmetic), i.e. running it twice is the same as running it once. 
removeCollinear = function(pointSet, tol = TOL) {
  N <- nrow(pointSet)
  if (N <= 2) {
    return(pointSet)
  }
  badInds <- rep(FALSE, N)
  for (ind in 2:(N - 1)) {
    if (abs(getAngle(pointSet[(ind - 1):(ind + 1), ])) < tol) {
      badInds[ind] <- TRUE
    }
  }
  keepPoints <- which(!badInds)
  pointSet %<>%
    slice(keepPoints)
  pointSet
}

### Find the slope and y-intercept of a line given by a 2 x 2 matrix, whose rows are points and columns are x and y coordinates
### Note: special cases are defined in the dependent functions for horizontal lines (slope = 0) and vertical lines (slope = Inf)
findSlopeAndIntercept = function(endPoints) {
  stopifnot(nrow(endPoints) == 2)
  endPoints %<>%
    as.matrix
  delta <- endPoints[2, ] - endPoints[1, ]
  slope <- delta[2]/delta[1]
  if (is.finite(slope)) {
    intercept <- endPoints[1, 2] - slope * endPoints[1, 1]
  } else if (is.na(slope)) {
    stop("The two endpoints are the same!")
  } else {
    intercept <- endPoints[1, 1]
  }
  output <- c(slope, intercept)
  output
}

### Find the intersection between two lines, each specified by a pair of points (2 x 2 matrices; row = point, columns = coords)
### Note: does not check that the intersection point lies between the starting and ending points of each segment!
findIntersection = function(endPoints1, endPoints2) {
  line1 <- findSlopeAndIntercept(endPoints1)
  line2 <- findSlopeAndIntercept(endPoints2)
  if (is.finite(line1[1]) && is.finite(line2[1])) {
    delta <- line2 - line1
    if (near(delta[1], 0)) {
      stop("The lines are parallel or near parallel!")
    }
    x0 <- -delta[2]/delta[1]
    y0 <- line1[1] * x0 + line1[2]
  } else if (is.finite(line1[1])) {
    x0 <- line2[2]
    y0 <- line1[1] * x0 + line1[2]
  } else if (is.finite(line2[1])) {
    x0 <- line1[2]
    y0 <- line2[1] * x0 + line2[2]
  } else {
    stop("The lines are both vertical!")
  }
  output <- c(x0, y0)
  output
}

### Compares two points numerically, up to a specified number of dimensions Dim (which can be smaller than the full dimension)
### Returns TRUE if they are numerically identical, FALSE otherwise. Does not check that the points have the same dimensions!
comparePoints = function(point1, point2, Dim = 2) {
  comp <- TRUE
  for (ind in 1:Dim) {
    if (!near(point1[ind], point2[ind], tol = TOL)) {
      comp <- FALSE
    }
  }
  comp
}

### Auxiliary function for the Imai-Iri algorithm
computeIndex <- function(index, sign, N) {
  allSigns <- c("positive", "negative", "computed", "other")
  return(index + (which(allSigns == sign) - 1) * N)
}

### Auxiliary function for the Imai-Iri algorithm
getIndex <- function(point) {
  return(point %>% pull(rowid))
}

### Auxiliary function for the Imai-Iri algorithm
getPoint <- function(pointSet, index) {
  return(pointSet %>% slice(index))
}

### Auxiliary function for the Imai-Iri algorithm
getPrevPoint <- function(pointSet, prevMap, point) {
  return(getPoint(pointSet, prevMap[getIndex(point)]))
}

### Auxiliary function for the Imai-Iri algorithm
getNextPoint <- function(pointSet, nextMap, point) {
  return(getPoint(pointSet, nextMap[getIndex(point)]))
}

### Auxiliary function for the Imai-Iri algorithm
setPoint <- function(pointSet, index, point) {
  pointSet[index, ] <- point
  pointSet
}

### Auxiliary function for the Imai-Iri algorithm
setPrevPoint <- function(prevMap, curPoint, prevPoint) {
  prevMap[getIndex(curPoint)] <- getIndex(prevPoint)
  prevMap
}

### Auxiliary function for the Imai-Iri algorithm
setNextPoint <- function(nextMap, curPoint, nextPoint) {
  nextMap[getIndex(curPoint)] <- getIndex(nextPoint)
  nextMap
}

### Auxiliary function for the Imai-Iri algorithm
makePoint <- function(point, index) {
  output <- tibble(Total = point[1], Sum = point[2], rowid = index)
  output
}

### Auxiliary function for the Imai-Iri algorithm
getSignedAngle <- function(point1, point2, point3, sign) {
  myPoints <- bind_rows(point1, point2, point3) %>%
    select(-rowid) %>%
    as.matrix
  myAngle <- getAngle(myPoints) * SIGNS[sign]
}

### Auxiliary function for the Imai-Iri algorithm
checkAcute <- function(point1, point2, point3, sign, tol = TOL) {
  myAngle <- getSignedAngle(point1, point2, point3, sign)
  return(myAngle > tol)
}

### Auxiliary function for the Imai-Iri algorithm
checkObtuse <- function(point1, point2, point3, sign, tol = TOL) {
  myAngle <- getSignedAngle(point1, point2, point3, sign)
  return(myAngle < -tol)
}

### Implementation of the Imai-Iri algorithm for finding the shortest piecewise linear path lying between pointsP and pointsM.
### Assumes, without checking, that the x-coordinates of pointsP and pointsM are identical, and ordered from smallest to largest.
### Also assumes, without checking, that the y-coordinates of pointsP are pointwise larger than those of pointsM.
### This implementation is a corrected version of the pseudocode in Sabine Neubauer's student thesis (University of Karlsruhe).
ImaiIriAlgorithm = function(pointsP, pointsM) {
  n <- nrow(pointsP)
  stopifnot(nrow(pointsM) == n)
  stopifnot(ncol(pointsP) == 2)
  stopifnot(ncol(pointsM) == 2)
  if (n <= 2) {
    return(pointsM)
  }
  # Global numbering scheme: the positive points are 1:n; the negative points are (n+1):2n; the computed Q points are (2n+1):4n
  allQ <- matrix(NA, n, 2, dimnames = list(c(), colnames(pointsP))) %>%
    as_tibble
  allPoints <- bind_rows(pointsP, pointsM, allQ, allQ) %>%
    rowid_to_column %>%
    select(-rowid, rowid)
  nextMap <- rep(NA, nrow(allPoints))
  names(nextMap) <- 1:nrow(allPoints)
  prevMap <- rep(NA, nrow(allPoints))
  names(prevMap) <- 1:nrow(allPoints)
  allP <- vector("list")
  allL <- vector("list")
  allR <- vector("list")
  for (sign in c("positive", "negative")) {
    nextMap %<>%
      setNextPoint(getPoint(allPoints, computeIndex(1, sign, n)), getPoint(allPoints, computeIndex(2, sign, n)))
    prevMap %<>%
      setPrevPoint(getPoint(allPoints, computeIndex(2, sign, n)), getPoint(allPoints, computeIndex(1, sign, n)))
    allP[[sign]] <- getPoint(allPoints, computeIndex(1, sign, n))
    allL[[sign]] <- getPoint(allPoints, computeIndex(1, sign, n))
    allR[[sign]] <- getPoint(allPoints, computeIndex(2, sign, n))
  }
  j <- 1
  print(paste("There are", n, "boundary points to process"))
  for (i in 3:n) {
    if (i %% 500 == 0) {
      print(paste("Processed", i, "boundary points so far"))
    }
    nextWindow <- FALSE
    for (sign in c("positive", "negative")) {
      curPrev  <- getPoint(allPoints, computeIndex(i - 1, sign, n))
      curPoint <- getPoint(allPoints, computeIndex(i, sign, n))
      curP     <- allP[[sign]]
      while (!comparePoints(curPrev, curP) && checkObtuse(curPoint, curPrev, getPrevPoint(allPoints, prevMap, curPrev), sign)) {
        curPrev <- getPrevPoint(allPoints, prevMap, curPrev)
      }
      nextMap %<>%
        setNextPoint(curPrev, curPoint)
      prevMap %<>% 
        setPrevPoint(curPoint, curPrev)
    }
    for (index in 1:2) {
      signStar <- ifelse(index == 1, "positive", "negative")
      signDiam <- ifelse(index == 1, "negative", "positive")
      curPoint <- getPoint(allPoints, computeIndex(i, signStar, n))
      curL <- allL[[signStar]]
      curR <- allR[[signDiam]]
      if (!nextWindow && checkAcute(curPoint, curL, curR, signStar)) {
        indQ <- computeIndex(j, "computed", n)
        myIntersect <- makePoint(findIntersection(bind_rows(curL, curR), bind_rows(allP[[signStar]], allP[[signDiam]])), indQ)
        allPoints %<>%
          setPoint(indQ, myIntersect)
        allP[[signDiam]][1,] <- allR[[signDiam]][1,]
        prevP <- getPoint(allPoints, computeIndex(i - 1, signStar, n))
        curP  <- getPoint(allPoints, computeIndex(i, signStar, n))
        indR <- computeIndex(j, "other", n)
        otherIntersect <- makePoint(findIntersection(bind_rows(curL, curR), bind_rows(prevP, curP)), indR)
        allPoints %<>%
          setPoint(indR, otherIntersect)
        allP[[signStar]][1,] <- otherIntersect
        j %<>%
          add(1)
        nextMap %<>%
          setNextPoint(otherIntersect, curP)
        prevMap %<>%
          setPrevPoint(curP, otherIntersect)
        allR[[signStar]][1,] <- getPoint(allPoints, computeIndex(i, signStar, n))
        allR[[signDiam]][1,] <- getPoint(allPoints, computeIndex(i, signDiam, n))
        allL <- allP
        while (checkAcute(allL[[signDiam]], allR[[signStar]], getNextPoint(allPoints, nextMap, allL[[signDiam]]), signDiam)) {
          allL[[signDiam]] <- getNextPoint(allPoints, nextMap, allL[[signDiam]])
        }
        nextWindow <- TRUE
      }
    }
    if (!nextWindow) {
      for (index in 1:2) {
        signStar <- ifelse(index == 1, "positive", "negative")
        signDiam <- ifelse(index == 1, "negative", "positive")
        curP <- getPoint(allPoints, computeIndex(i, signStar, n))
        if (checkAcute(curP, allL[[signDiam]], allR[[signStar]], signStar)) {
          allR[[signStar]][1,] <- getPoint(allPoints, computeIndex(i, signStar, n))
          while (checkAcute(curP, allL[[signDiam]], getNextPoint(allPoints, nextMap, allL[[signDiam]]), signStar)) {
            allL[[signDiam]][1,] <- getNextPoint(allPoints, nextMap, allL[[signDiam]])
          }
        }
      }
    }
  }
  m <- j + 1
  if (comparePoints(allL[["positive"]], allP[["positive"]])) {
    supportLine <- rbind(allL[["positive"]], allR[["negative"]])
  } else {
    supportLine <- rbind(allL[["negative"]], allR[["positive"]])
  }
  indQ <- computeIndex(m - 1, "computed", n)
  allPoints %<>%
    setPoint(indQ, makePoint(findIntersection(supportLine, rbind(allP[[1]], allP[[2]])), indQ))
  indP <- computeIndex(n, "positive", n)
  indM <- computeIndex(n, "negative", n)
  indR <- computeIndex(m, "computed", n)
  allPoints %<>%
    setPoint(indR, makePoint(findIntersection(supportLine, rbind(getPoint(allPoints, indP), getPoint(allPoints, indM))), indR))
  outputInds <- computeIndex(1, "computed", n):indR
  output <- allPoints %>%
    slice(outputInds) %>%
    select(-rowid)
  output
}

### This function finds all the possible non-zero partial sums of input listOfSizes that do not exceed the specified upperBound
findPartialSums = function(listOfSizes, upperBound) {
  fullPos <- rep(FALSE, upperBound)
  L <- length(listOfSizes)
  if (L > 0) {
    for (ind in 1:L) {
      curValue <- listOfSizes[ind]
      altPos <- c(rep(FALSE, curValue), head(fullPos, -curValue)) # result of adding the current value to previous partial sums
      fullPos %<>%
        magrittr::or(altPos)
    }
  }
  return(which(fullPos))
}

### This function computes, for a given N and nPlus, the set of critical value pairs (i, j) such that the 2x2 contingency table
### with row sums nPlus and N - nPlus, and column sums i and N - i has p-value below pMax when its top left entry is at least j. 
computePValueBounds = function(N, nPlus, pMax, logSpace = FALSE) {
  print(paste("The", rep("log", logSpace), "p-value being considered as the cutoff is", pMax))
  print(paste("Preparing p-value bounds with N =", N, "and nPlus =", nPlus))
  stopifnot(nPlus <= N)
  stopifnot((logSpace && is.finite(pMax) && pMax <= 0) || (!logSpace && 0 <= pMax && pMax <= 1))
  nMinus <- N - nPlus
  fast <- TRUE
  if (nPlus > nMinus) {
    fast <- FALSE ### The conditions for speeding up the computation may not be met
  }
  if (nPlus == 0) {
    return(tibble(Sum = integer(), Total = integer()))
  }
  iRange <- 1:nPlus
  criticalValues <- rep(NA, nPlus)
  prevValue <- 0
  for (i in iRange) {
    if (!fast) {
      myFunction <- function(x) { phyper(x, nPlus, nMinus, i, lower.tail = FALSE, log.p = logSpace) }
      result <- gtools::binsearch(myFunction, range = c(max(0, i - nMinus), min(nPlus, i)), target = pMax)
      if (result$flag == "Upper Boundary") { ### unattainable p-value
        prevValue <- NA
      } else if (result$flag == "Lower Boundary") { ### p-value attained at the lower end of the range
        prevValue <- max(0, i - nMinus)
      } else {
        prevValue <- result$where[1] + 1 ### placing the correct value where it belongs - note indexing is 0-based
      }
    } else {
      range <- max(prevValue - M_STEP, 0, na.rm = TRUE) : min(prevValue + M_STEP, i, na.rm = TRUE) # 0:nPlus
      hyperVector <- phyper(range, nPlus, nMinus, i, lower.tail = FALSE, log.p = logSpace)
      prevValue <- range[min(which(hyperVector < pMax))]
    }
    criticalValues[i] <- prevValue + 1 ### the +1 is due to the offset of phyper with respect to fisher.test
  }
  output <- tibble(Sum = iRange, Total = criticalValues) %>%
    filter(Sum >= Total)
  output
}

unitTestImaiIri = function(completeTest = TRUE) {
  N <- c(8529, 99505, 9996, 9988, 9993, 9922, 1e4, 9736, 9998, 9997, 99589, 99818, 99165, 99839, 99547, 99210, 99173, 99983) 
  P <- c(1954, 18002, 153,  784,  593 , 792, 1745, 1016, 701 , 918 , 2758 , 29936, 5216 , 1741 , 22856, 3180,  16016, 5852 )
  N <- c(N,  99888, 1e5, 97905, 9983, 98400, 99898)
  P <- c(P,  3704 , 2835, 8675, 1414, 36894, 23070)
  tCases <- tibble(N = N, nPlus = P, pMax = 5e-8) %>%
    arrange(N)
  allRatios <- rep(0, nrow(tCases))
  for (ind in 1:nrow(tCases)) {
    curQ <- computePValueBounds(N = tCases$N[[ind]], nPlus = tCases$nPlus[[ind]], pMax = tCases$pMax[[ind]])
    if (completeTest) {
      A <- Sys.time()
      curR   <- callImaiIri(curQ, CPP = FALSE)
      B <- Sys.time()
      curAlt <- callImaiIri(curQ, CPP = TRUE)
      C <- Sys.time()
      t1 <- as.double(B - A, units = "secs")
      t2 <- as.double(C - B, units = "secs")
      ratio <- t1/t2
      print(paste("For input", ind, "the timings are", signif(t1, 3), "and", signif(t2, 3), "for a ratio of", round(ratio)))
      allRatios[ind] <- ratio
    } else {
      curR <- callImaiIri(curQ, dryRun = TRUE)[[3]]
      write_csv(curR, paste0("TestCase", ind, ".csv"))
    }
  }
  if (completeTest) {
    return(allRatios)
  }
}

testPythonCode = function(startFile = "Test.csv") {
  fullStart <- normalizePath(startFile)
  outputFile <- str_replace(fullStart, ".csv", "OutputC.csv")
  a <- Sys.time()
  cmdName <- "/Users/lchindelevitch/Downloads/ReverseGWAS/rgwas/plfoptq"
  system2(command = cmdName, args = c(fullStart, WIDTH/2, ">", outputFile))
  b <- Sys.time() - a
  f <- file(outputFile, 'r')
  Lines <- readLines(f)
  close(f)
  badLine <- ifelse(any(Lines == "()"), max(which(Lines == "()")), 0)
  Lines %<>%
    magrittr::extract((badLine[1] + 1):length(Lines)) %>%
    magrittr::extract(nchar(.) > 0) %>%
    str_remove_all("\\[") %>%
    str_remove_all("\\]") %>%
    str_trim() %>%
    # str_split_fixed("[ ]+", n = 2)
    str_split_fixed(",", n = 2)
  # Lines <- read_csv(outputFile)
  points <- matrix(as.numeric(Lines), ncol = 2) %>%
    set_colnames(c("Total", "Sum")) %>%
    as_tibble
  points <- points %>%
    mutate_at("Sum", ~{add(., WIDTH/2)})
  Tab <- read_csv(startFile)
  TabP <- Tab %>%
    mutate_at("Sum", ~{add(., WIDTH)})
  TabM <- Tab
  Tab %<>%
    removeCollinear
  TabPR <- Tab %>%
    mutate_at("Sum", ~{add(., WIDTH)})
  TabMR <- Tab
  A <- Sys.time()
  goldStandard <- ImaiIriAlgorithm(TabPR, TabMR)
  B <- Sys.time() - A
  ratio <- as.double(B, units = "secs")/as.double(b, units = "secs")
  print(paste("On input file", startFile, "the speed-up factor is", round(ratio)))
  redPoints   <- head(points, -2)
  redStandard <- head(goldStandard, -2)
  print(max(abs(redPoints$Sum - redStandard$Sum)))
  print(max(abs(redPoints$Total - redStandard$Total)))
  checkPosition(TabP, points, below = TRUE,  reportCoincident = TRUE)
  checkPosition(TabM, points, below = FALSE, reportCoincident = FALSE)
  ratio
}

wrapperPythonTest = function(indMin = 1, indMax = 21) {
  # unitTestImaiIri(completeTest = FALSE)
  allRatios <- c()
  for (ind in indMin:indMax) {
    curFile <- paste0("TestCase", ind, ".csv")
    allRatios %<>%
      c(testPythonCode(curFile))
  }
  allRatios
}

### OBSOLETE FUNCTIONS FROM VALIDATION.R

checkNesting = function(inputFile) {
  Tab <- read_csv(inputFile)
  if ("status_Validation" %in% colnames(Tab)) {
    allPairs <- Tab %>%
      mutate(survived = !is.na(status_Validation)) %>%
      select(c(survived, K, L))
    survivorInds <- which(allPairs$survived)
    nonSurvivorInds <- setdiff(1:nrow(allPairs), survivorInds)
    for (ind1 in survivorInds) {
      row1 <- allPairs %>% 
        slice(ind1)
      for (ind2 in nonSurvivorInds) {
        row2 <- allPairs %>% 
          slice(ind2)
        if (row2$K >= row1$K && row2$L >= row1$L) {
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

checkPVals = function(inputFile, minValue = log(5e-8), ratio = 1.1, stringent = TRUE) {
  Tab <- read_csv(inputFile)
  if ("status_Validation" %in% colnames(Tab)) {
    allRows <- Tab %>%
      mutate(survived = !is.na(status_Validation)) %>%
      select(c(survived, LogFisherExactP_Discovery, bestSingleStat_Discovery)) %>%
      filter(survived)
    altBound <- allRows %>%
      pull(bestSingleStat_Discovery)
    if (!all(near(altBound, altBound[1]))) {
      return(FALSE)
    }
    curBound <- ifelse(stringent, min(minValue, altBound[1] - log(ratio)), max(minValue, altBound[1] - log(ratio)))
    foundVals <- allRows %>%
      pull(LogFisherExactP_Discovery)
    return(all(foundVals <= curBound))
  }
  return(TRUE)
}

checkFormulas = function(inputFile) {
  Tab <- read_csv(inputFile)
  if ("status_Validation" %in% colnames(Tab)) {
    allRows <- Tab %>%
      mutate(survived = !is.na(status_Validation)) %>%
      select(c(survived, formula_Discovery)) %>%
      filter(survived)
    allFormulas <- allRows %>%
      pull(formula_Discovery)
    goodFormulas <- str_detect(allFormulas, "AND") | str_detect(allFormulas, "OR")
    return(all(goodFormulas))
  }
  return(TRUE)
}

checkAll = function(testDir, checkFunction = checkNesting) {
  failures <- c()
  initDir <- getwd()
  setwd(testDir)
  allFiles <- list.files()
  for (File in allFiles) {
    curResult <- checkFunction(File)
    if (!curResult) {
      failures %<>%
        c(File)
    }
  }
  setwd(initDir)
  return(failures)
}

### Obsolete functions rom Driver.R

checkResultConsistency = function(resultDir1 = "AsaData/InputFiles/", resultDir2 = "AsaData/ProcessedFilesFastComplementPhenotypesBiobankSize/") {
  initDir <- getwd()
  List1 <- list.files(path = resultDir1) %>%
    str_subset("Summary")
  subList1 <- List1 %>%
    str_sub(end = str_locate(., "agreement")[,2])
  List2 <- list.files(path = resultDir2) %>%
    str_subset("Summary")
  subList2 <- List2 %>%
    str_sub(end = str_locate(., "agreement")[,2])
  myMatches <- match(subList2, subList1)
  for (ind in 1:length(List2)) {
    print(ind)
    curTab1 <- read_csv(paste0(resultDir2, List2[[ind]]))
    curTab2 <- read_csv(paste0(resultDir1, List1[[myMatches[ind]]]))
    stopifnot(all(curTab1$SNP == curTab2$SNP))
    status1 <- curTab1$status %>%
      na_if("timeout")
    status2 <- curTab2$status %>%
      recode("feasible" = "optimal")
    stopifnot(all(status1==status2, na.rm=TRUE))
  }
  setwd(initDir)
  return(TRUE)
}

optimizingMetaDriver = function(inputDir = "AsaData/InputFiles/") {
  initDir <- getwd()
  setwd(inputDir)
  List <- list.files() %>%
    str_subset("Summary")
  fullResults <- vector("list", length(List)) %>%
    set_names(outer(c("CNF", "DNF"), 1:(length(List)/2), paste0))
  for (ind in 1:length(List)) {
    curFile <- List[[ind]]
    curType <- ifelse(str_count(curFile, "CNF") > 0, "CNF", "DNF")
    curName <- str_extract(curFile, "_chr[0-9]+_") %>%
      str_sub(5, -2) %>%
      str_c(curType, ., sep = "")
    altFile <- paste0(str_sub(curFile, end = str_locate(curFile, "NF_")[,2] - 5), ".csv")
    print(paste("Currently processing", altFile))
    curTab <- read_csv(List[[ind]]) %>%
      filter(status %in% c("feasible", "optimal")) %>%
      select(c(SNP, LogFisherExactP))
    numSurvivors <- nrow(curTab)
    allResults <- vector("list", numSurvivors) %>%
      set_names(curTab$SNP)
    if (numSurvivors > 0) {
      print(paste("There are", numSurvivors, "SNPs to process!"))
      for (index in 1:numSurvivors) {
        print(index)
        curRow <- curTab %>% 
          slice(index)
        curSNP <- curRow$SNP
        curPValue <- curRow$LogFisherExactP
        allResults[[curSNP]] <- selectOptimizingDriver(inputFile = altFile, SNPs = curSNP, startPValue = curPValue, type = curType, Log = TRUE)
      }
    }
    fullResults[[curName]] <- allResults
  }
  setwd(initDir)
  save(fullResults, file = "AllOptimizationResults.RData")
  fullResults
}

selectOptimizingDriver = function(inputFile = "CPLEXInputUKBB_autoimmune_11_EUR_chr1_all_unrelated_281591_122.csv", SNPs = c("rs6679677_A"), startPValue = 5e-8, type = "CNF", Log = FALSE) {
  ext <- str_sub(inputFile, start = -4)
  readFunction <- ifelse(ext == '.tsv', read_tsv, read_csv)
  inputTab <- readFunction(inputFile, col_types = cols(ID = "c", .default = "l"))
  numSNPs <- str_extract(inputFile, paste0("([0-9]+)", ext)) %>%
    str_remove(ext) %>%
    parse_integer
  phenoNames <- colnames(inputTab)[-(1:(numSNPs + 1))]
  redTab <- inputTab %>%
    magrittr::extract(, c("ID", SNPs, phenoNames))
  redFilename <- str_replace(inputFile, ext, paste0(SNPs, "_", 1, ext))
  write_csv(redTab, redFilename)
  result <- optimizingDriver(redFilename, type = type, complement = -1, startPValue = startPValue, Log = Log)
  result
}

sampleOptimizingDriver = function() {
  source("MetaDriver.R")
  inFile <- subsampleFile("CPLEXInput_AutoImmune_EU_100000_15.tsv", nPheno = 4, nPat = 100000, Snp = "rs2413583_C", Dir = FALSE)
  output <- optimizingDriver(inFile)
  output
}
