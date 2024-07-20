options(warn = 0)
library(expm)
library(magrittr)
library(Matrix)
library(tidyverse)
OPTIMAL = "optimal"      ### status descriptor
OPTIMIZATION = TRUE
SINGLE_BEST_RATIO = 1.1  ### target ratio of the combined phenotype's summary statistic to that of the best singleton phenotype

mainDriver = function(inputFile, type = "LCP", objective = "correlation", K = Inf, complement = 0, Log = FALSE, ext = ".csv", index = 1L, 
                      extremeValue = 0, outputSummary = TRUE, outputPhenotype = FALSE, outputAssociations = FALSE, saveOutput = FALSE) {
  inputTab <- read_csv(inputFile, col_types = cols(ID = "c", id = "c", .default = "d"))
  if ("id" %in% colnames(inputTab)) { 
    inputTab %<>% 
      rename(ID = id) 
  }
  numSNPs <- str_extract(inputFile, paste0("([0-9]+)", ext)) %>% 
    str_remove(ext) %>% 
    parse_integer()
  myP <- extremeValue %>% 
    as.character
  myF <- ifelse(!is.na(SINGLE_BEST_RATIO), paste0("F", SINGLE_BEST_RATIO), "")
  baseFilename <- str_remove(inputFile, ext)
  stage <- ifelse(OPTIMIZATION, "Optimization", "Triage")
  specString <- paste(c("", type, paste0("K", K), paste0(ifelse(type == "LCP", "R", "P"), myP), myF, objective, 
                        stage, "Results"), collapse = "_")
  output <- createAndSolveILPs(inputTab, numSNPs = numSNPs, type = type, objective = objective, K = K, index = index, 
                               extremeValue = extremeValue, complement = complement, baseFilename = baseFilename, outputAssociations = FALSE)
  if (outputSummary)      { write_csv(output$summary,       str_replace(inputFile, ext, str_c(specString,     '_Summary.csv')))}
  if (outputPhenotype)    { write_csv(output$phenotype,     str_replace(inputFile, ext, str_c(specString,   '_Phenotype.csv')))}
  if (outputAssociations) { write_csv(output$associations,  str_replace(inpuFile, ext, str_c(specString, '_Associations.csv')))}
  if (saveOutput)         { saveRDS(object = output, file = str_replace(inputFile, ext, str_c(specString,             '.RDS')))}
  output
}

### This function creates and solves ILPs (integer linear programs) for maximizing an objective function over input phenotypes
### The first numSNPs columns of the input table are assumed to be the SNPs of interest; the rest are assumed to be phenotypes
### The value of objective can be "agreement", "covariance", or "correlation"
### The value of type can be "CNF" (conjunctive normal form), "DNF" (disjunctive normal form), or "LCP" (a linear combination)
### The value of K determines the maximum allowed number of combinations (for CNF or DNF) or phenotypes (for LCP) in the output
### If complement =  1 or 2, the genotypes  are augmented with their negation (col name prefix "NOT") before processing starts.
### If complement = -1 or 2, the phenotypes are augmented with their negation (col name prefix "NOT") before processing starts.
### baseFilename is used for storing the boundary when p-value constraints are present; the index is used for distinguishing them
### extremeValue is a p-value upper bound, above which any solutions found are declared infeasible; this leads to default results
### The program assumes that there may be NAs in the genotypes, but NOT in the phenotypes; if that is not the case it will fail!
### The NAs in each individual genotype are removed, alongside the corresponding phenotype rows, during the preprocessing stage
### The return value is a list with two components; the first one is a table of statistics which includes the optimal formulas,
### and the second one is a matrix containing the complex phenotypes corresponding to each of these formulas, which may be empty
### If outputAssociations = TRUE, then the output additionally contains the matrix of single genotype - phenotype associations
createAndSolveILPs =  function(inputTab, numSNPs, objective = "correlation", type = "LCP", K = Inf, complement = 0,
                               baseFilename = NULL, index = 0, outputAssociations = FALSE, extremeValue = 0) {
  stopifnot(type == "LCP" && objective == "correlation" && complement == 0)
  IDs <- inputTab %>% 
    pull("ID")
  inputTab %<>% 
    select(-ID) %>%
    as.matrix %>%
    set_rownames(IDs)
  phenotypes <- inputTab[, -(1:numSNPs), drop = FALSE]
  stopifnot(all(!is.na(phenotypes)))
  genotypes <- inputTab[, 1:numSNPs, drop = FALSE] 
  print(paste("There are", ncol(phenotypes), "phenotypes over", nrow(phenotypes), "patients"))
  print(paste("There are", ncol(genotypes), "genotypes over", nrow(genotypes), "patients"))
  origPhenotypes <- phenotypes
  origGenotypes  <- genotypes
  output <- vector("list", ncol(genotypes))
  names(output) <- colnames(genotypes)
  complexPs <- matrix(NA, nrow(phenotypes), ncol(genotypes), dimnames = list(rownames(phenotypes), names(output)))
  coeffsTrue  <- rep(1, ncol(genotypes)) 
  coeffsAgree <- rep(1, ncol(genotypes))
  M <- nrow(genotypes)
  N <- ncol(genotypes)
  objVectors  <- matrix(coeffsAgree, M, N, byrow = TRUE) * genotypes + matrix(coeffsTrue, M, N, byrow = TRUE)
  print("Reducing the input matrices")
  out <- reducePhenotypesAndObjectives(phenotypes, objVectors)
  phenotypes <- out$phenotypes
  objVectors <- out$objVectors
  print("Computing single associations")
  sRes <- computeBestPValues(origGenotypes, origPhenotypes, cor = TRUE, fisher = FALSE)
  boundValues   <- rep(NA, N)
  genoNames     <- colnames(genotypes)
  print(paste("There are", N, "SNPs to process"))
  for (ind in 1:N) {
    print(paste("Currently processing SNP", colnames(genotypes)[ind]))
    curName <- genoNames[ind]
    curObjVector <- objVectors[, ind, drop = FALSE]
    if (any(is.na(curObjVector))) { 
      print(paste("Removing", sum(is.na(curObjVector)), "missing entries in", curName)) 
    }
    goodInds <- which(!(is.na(curObjVector)))
    P <- phenotypes[goodInds, , drop = FALSE]
    G <- genotypes[goodInds, ind, drop = FALSE]
    initBound <- boundValues[ind]
    curResult  <- processLCP(G, P, K)
    curSol     <- curResult$solution
    upperBound <- curResult$bound
    boundValues[ind] <- upperBound
    curPheno   <- postprocessSolution(origPhenotypes, curSol$usedVars, type = type, complement = complement)
    curPhenoMini <- curPheno$solution
    curGenoMini  <- origGenotypes[, ind, drop = TRUE]
    checkOpt     <- cor(curGenoMini, curPhenoMini[,1], use = "pairwise")
    if (!is.na(checkOpt) && checkOpt < 0) { ### implicitly multiply the coefficients by -1
      checkOpt     <- -checkOpt
      curPhenoMini <- -curPhenoMini
      curPheno$clause <- negateCoefficients(curPheno$clause)
    }
    curOpt    <- curSol$optimum
    curOutput <- list(time = curSol$time, gap = curSol$gap, opt = curOpt, status = curSol$status, formula = curPheno$clause)
    output[[curName]]    <- curOutput
    complexPs[, curName] <- curPhenoMini
  }
  fullOutput <- prepareOutput(output, complexPs, genoNames, IDs, type, objective, sRes, boundValues, outputAssociations)
  fullOutput
}

### This function reduces an input consisting of a phenotype matrix and a matrix of objective function vectors (same # of rows)
### The return value is a list containing:
### - the reduced phenotype and objective matrices
### - the status of each row and column: FALSE if it has been removed, TRUE if it has been preserved
reducePhenotypesAndObjectives = function(phenotypes, objVectors) {
  stopifnot(nrow(objVectors) == nrow(phenotypes))
  rowStatuses <- rep(TRUE, nrow(phenotypes))
  colStatuses <- rep(TRUE, ncol(phenotypes))
  cSums <- colSums(phenotypes)
  allTrueCols  <- which(cSums == nrow(phenotypes))
  allFalseCols <- which(cSums == 0)
  allConstantCols <- c(allTrueCols, allFalseCols)
  if (length(allConstantCols) > 0) {
    colStatuses[allConstantCols] <- FALSE
    phenotypes  %<>% 
      magrittr::extract(, -allConstantCols, drop = FALSE)
  }
  activeCols <- which(colStatuses)
  uniqueCols <- extractUniqueRows(t(phenotypes))
  uniqueGroupsC <- uniqueCols[[2]]
  phenotypes    <- t(uniqueCols[[1]])
  for (ind in seq_along(uniqueGroupsC)) {
    curGroup <- uniqueGroupsC[[ind]]
    colStatuses[activeCols[curGroup]] <- ind * c(1L, rep(-1L, length(curGroup) - 1))
  }
  print(paste("After reduction, there are", ncol(phenotypes), "phenotypes over", nrow(phenotypes), "patients"))
  print(paste("After reduction, there are", ncol(objVectors),  "genotypes over", nrow(objVectors), "patients"))
  output <- list(phenotypes = phenotypes, objVectors = objVectors, rowStatuses = rowStatuses, colStatuses = colStatuses)
  output
}

### Partly based on the solution at r.789695.n4.nabble.com/Count-unique-rows-columns-in-a-matrix-td844731.html
extractUniqueRows = function(Table, repeats = TRUE) {
  Q <- nrow(Table)
  if (Q <= 1) {
    output <- list(Table)
    if (repeats) {
      output <- c(output, list(Q))
    }
    return(output)
  }
  myOrder <- do.call(order, as.data.frame(Table))
  Table <- Table[myOrder, , drop = FALSE]
  equalToPrevious <- rep(FALSE, Q - 1)
  prevRow <- Table[1,]
  index <- 1
  while (index < Q) {
    curRow <- Table[index + 1, ]
    curComp <- (all(is.na(prevRow) == is.na(curRow)) && all(prevRow[!is.na(prevRow)] == curRow[!is.na(curRow)]))
    if (curComp) {
      equalToPrevious[index] <- TRUE
    }
    prevRow <- curRow
    index <- index + 1
  }
  lowerBounds <- c(0, which(!equalToPrevious)) + 1
  output <- list(Table[lowerBounds, , drop = FALSE])
  if (repeats) {
    N <- nrow(output[[1]])
    repList <- vector("list", N)
    upperBounds <- c(lowerBounds[-1] - 1, Q)
    repList <- lapply(1:N, function(x) { 
      myOrder[lowerBounds[x]:upperBounds[x]] 
    })
    output <- c(output, list(repList))
  }
  output
}

### Compute the column with best metric in matrix2 for each column of matrix1, returning both its value (arg1) and name (arg2).
### The matrices are assumed (without checking) to have the same number of rows and, unless cor = TRUE, to also be binary (0/1).
### If cor = TRUE use abs(correlation) as the metric; else if fisher = TRUE use the Fisher exact, otherwise chi-squared p-value.
computeBestPValues = function(matrix1, matrix2, fisher = FALSE, cor = TRUE) {
  M <- nrow(matrix1)
  stopifnot(nrow(matrix2) == M)
  N1 <- ncol(matrix1)
  bestPs <- rep(NA, N1)
  bestps <- rep("", N1)
  bestQs <- rep(NA, N1)
  N2 <- ncol(matrix2)
  rown <- colnames(matrix1)
  coln <- colnames(matrix2)
  selectorFun <- ifelse(cor, max, min)
  print(paste("There are", N1, "columns to process"))
  allValues <- abs(cor(matrix1, matrix2, use = "pairwise.complete.obs", method = "pearson"))
  binTypes  <- map_lgl(1:N2, function(x) { curCol <- matrix2[,x]; curSize <- n_distinct(curCol); curSize <= 2 } )
  for (ind1 in 1:N1) {
    if (ind1 %% 100 == 0) { print(ind1) }
    curPVals     <- allValues[ind1, ]
    bestValue    <- selectorFun(curPVals, na.rm = TRUE)
    bestPs[ind1] <- bestValue
    curBest <- which(near(curPVals, bestValue, tol = .Machine$double.eps^0.5 * M))[1] ### improvement: fixed tol
    bestps[ind1] <- coln[curBest] 
    if (binTypes[curBest]) {
      bestQs[ind1] <- chisq.test(table(matrix1[, ind1], matrix2[, curBest]))$p.value
    } else {
      bestQs[ind1] <- cor.test(matrix1[, ind1], matrix2[, curBest], method = "pearson", alternative = "g")$p.value
    }
  }
  dimnames(allValues) <- list(rown, coln)
  output <- list(value = bestPs, phenotype = bestps, allValues = allValues, significance = bestQs)
  output
}

### Auxiliary function to process an LCP problem given an input genotype G, phenotype P, and maximum number of phenotypes K.
processLCP = function(G, P, K) {
  stopifnot(!is.null(colnames(P)))
  Pc <- scale(P, center = TRUE, scale = FALSE)
  Gc <- scale(G, center = TRUE, scale = FALSE)
  stopifnot(K >= ncol(Pc))
  Z <- cancor(Gc, Pc, xcenter = FALSE, ycenter = FALSE)
  usedVars <- Z$ycoef[colnames(Pc), 1, drop = FALSE]
  usedVars <- usedVars/max(abs(usedVars))
  curSol <- list(usedVars = usedVars, optimum = Z$cor, status = OPTIMAL, time = 0, gap = 0)
  output <- list(solution = curSol, bound = Z$cor)
  output
}

### This function postprocesses a given solution of an ILP/MIQCP to put it into a convenient format: the phenotype and the clause
postprocessSolution = function(origPhenotypes, solutionMatrix, type, complement = 0) {
  if (any(is.na(solutionMatrix)) || all(near(solutionMatrix, 0))) {
    return(list(solution = matrix(0, nrow(origPhenotypes), 1), clause = ""))
  }
  goodPos <- which(!(near(solutionMatrix, 0)))
  curSolution <- origPhenotypes[, rownames(solutionMatrix), drop = FALSE] %*% solutionMatrix
  solutionMatrix <- solutionMatrix[goodPos, 1, drop = FALSE]
  solOrder <- order(abs(solutionMatrix[,1]), decreasing = TRUE)
  clause <- paste0("(", solutionMatrix[solOrder, 1], ") * ", rownames(solutionMatrix)[solOrder], collapse = " + ")
  output <- list(solution = curSolution, clause = clause)
  output
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

### This function prepares the final output from a sequence of ILPs, one per SNP, in the form of 2 tables: summary and phenotype 
prepareOutput = function(output, complexPhenotypes, SNP, ID, type, objective, sRes, boundValues, log, outputAssociations) {
  myTab <- matrix(unlist(output), nrow = length(output), dimnames = list(names(output), names(output[[1]])), byrow = TRUE)
  myTab %<>%
    as_tibble %>%
    mutate(SNP = SNP) %>%
    select(SNP, everything())
  myTab %<>%
    mutate_at(c("time", "gap", "opt"), parse_double)
  myTab %<>%
    mutate(bestSingle = sRes$phenotype, bestSingleStat = sRes$value, computedBounds = boundValues, bestPValue = sRes$significance)
  if (!is.null(complexPhenotypes)) {
    complexPhenotypes %<>%
      as_tibble %>%
      mutate_all(as.numeric) %>%
      mutate(ID = ID) %>%
      select(ID, everything())
  }
  fullOutput <- list(summary = myTab, phenotype = complexPhenotypes)
  if (outputAssociations) { 
    fullOutput %<>% 
      c(list(associations = sRes$allValues)) 
  }
  fullOutput
}

parseFormula = function(Formula) {
  outerString   <- " \\+ "
  innerString   <- " \\* "
  splitFormula <- str_split(Formula, outerString)[[1]]
  fullySplitFormula <- map(splitFormula, function(x) { 
    x %>% 
      str_trim %>% 
      str_remove("^\\(") %>% 
      str_remove("\\)$") %>% 
      str_split(innerString) %>% 
      extract2(1) 
  })
  extractFunction <- function(x) {x[[2]]}
  extractCoeffs   <- function(x) {x[[1]] %>% str_remove("\\)$") %>% parse_double} 
  usedNames  <- map(fullySplitFormula, extractFunction)
  usedValues <- map(fullySplitFormula, extractCoeffs) 
  allUsedNames <- sort(unique(unlist(usedNames)))
  M <- length(allUsedNames)
  parsedFormula <- matrix(0, M, 1, dimnames = list(allUsedNames, 1))
  parsedFormula[unlist(usedNames), 1] <- unlist(usedValues)
  parsedFormula
}

testLCP = function(numPatients = 100000, numSNPs = 3, coeffs = 1:numSNPs, signs = (-1)^(coeffs - 1), extraSNPs = 3) {
  totalSNPs = numSNPs + extraSNPs
  inputMat  = cbind(matrix(1, numPatients), matrix(runif(numPatients * totalSNPs), ncol = totalSNPs))
  QR        = Matrix::qr(inputMat)
  transMat  = Matrix::qr.Q(QR)
  fullPheno = transMat[, 1 + (1:totalSNPs), drop = FALSE] %>%
    set_colnames(paste0("p", 1:totalSNPs))
  fullGeno  = fullPheno[, 1:numSNPs, drop = FALSE] %*% (coeffs * signs)
  IDvector  = paste0("P", 1:numPatients)
  fullInput = bind_cols(ID = IDvector, SNP1 = as.vector(fullGeno), fullPheno %>% as_tibble())
  inputFile = "TestLCP_1.csv"
  write_csv(fullInput, path = inputFile)
  result    = mainDriver(inputFile = inputFile)
  formula   = parseFormula(result$summary$formula)[,1]
  stopifnot(all(near(formula/min(abs(formula)), coeffs * signs/min(abs(coeffs * signs)))))
  stopifnot(all(sort(names(formula)) == sort(paste0("p", 1:numSNPs))))
  result
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
