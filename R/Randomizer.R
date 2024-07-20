KAPPA = 10  ### Expected number of swaps per edge in the simulation; larger values make it run slower

### This function produces a random member of the set of binary matrices with the same size, row and column sums as the given one
### If keepRownames = TRUE, the original row    names are used; otherwise, R1, R2, ..., Rm are the row    names where m = # rows
### If keepColnames = TRUE, the original column names are used; otherwise, C1, C2, ..., Cn are the column names where n = # cols

randomizeBinaryMatrix = function(binaryMatrix, keepRownames = FALSE, keepColnames = FALSE) {
  m <- nrow(binaryMatrix)
  n <- ncol(binaryMatrix)
  nnz <- sum(binaryMatrix)
  print(paste("Processing a", m, "by", n, "input matrix with", nnz, "nonzeros"))
  rSums <- rowSums(binaryMatrix)
  zeroRows <- which(rSums == 0)
  goodRows <- setdiff(1:m, zeroRows)
  randMatrix <- binaryMatrix[goodRows, , drop = FALSE]
  rSums <- rSums[goodRows]
  m0 <- length(rSums)
  cumSums <- cumsum(rSums)
  iter <- 1
  M <- tail(cumSums, 1)
  mRange <- 1:M
  maxIter <- KAPPA * M
  print(paste("There are", maxIter, "edge switches to perform"))
  while (iter < maxIter) {
    curInds  <- sample(mRange, 2)
    res1 <- gtools::binsearch(function(x) {cumSums[x]}, target = curInds[1], range = c(1, m0))
    row1 <- ifelse(res1$flag == "Lower Boundary", 1, tail(res1$where, 1))
    res2 <- gtools::binsearch(function(x) {cumSums[x]}, target = curInds[2], range = c(1, m0))
    row2 <- ifelse(res2$flag == "Lower Boundary", 1, tail(res2$where, 1))
    if (row1 == row2) {
      next
    }
    col1 <- min(which(cumsum(randMatrix[row1, ]) == curInds[1] - ifelse(row1 == 1, 0, cumSums[row1 - 1])))
    col2 <- min(which(cumsum(randMatrix[row2, ]) == curInds[2] - ifelse(row2 == 1, 0, cumSums[row2 - 1])))
    if (col1 == col2 || randMatrix[row1, col2] == 1 || randMatrix[row2, col1] == 1) {
      next
    }
    randMatrix[row1, col1] <- 0
    randMatrix[row2, col2] <- 0
    randMatrix[row1, col2] <- 1
    randMatrix[row2, col1] <- 1
    iter <- iter + 1
    if (iter %% 1000 == 0) {
      print(iter)
    }
  }
  binaryMatrix[goodRows, ] <- randMatrix
  if (!keepRownames) {
    rownames(binaryMatrix) <- paste0("R", 1:m)
  }
  if (!keepColnames) {
    colnames(binaryMatrix) <- paste0("C", 1:n)
  }
  binaryMatrix
}

randomizeInputFile = function(inputFile = "CPLEXInput_AutoImmune_EU_100000_15.tsv") {
  extension <- str_sub(inputFile, -4, -1)
  if (extension == ".csv") {
    inputTab <- read_csv(inputFile)
  } else {
    inputTab <- read_tsv(inputFile)
  }
  numSNPs  <- inputFile %>%
    str_extract_all("[0-9]+") %>%
    extract2(1) %>%
    tail(1) %>%
    as.integer()
  initMatrices <- prepareMatrices(inputTab, numSNPs = numSNPs, complement = FALSE, keepDuplicates = TRUE)
  phenotypes   <- initMatrices[[1]]
  stopifnot(!any(is.na(phenotypes)))
  if (sum(phenotypes) > prod(dim(phenotypes))/2) {
    phenotypesR   <- 1 - randomizeBinaryMatrix(1 - phenotypes,  keepColnames = TRUE)
  } else {  
    phenotypesR  <- randomizeBinaryMatrix(phenotypes, keepColnames = TRUE)
  }
  genotypes    <- initMatrices[[2]]
  genotypes[is.na(genotypes)] <- 0
  if (sum(genotypes) > prod(dim(genotypes))/2) {
    genotypesR   <- 1 - randomizeBinaryMatrix(1 - genotypes,  keepColnames = TRUE)
  } else {
    genotypesR   <- randomizeBinaryMatrix(genotypes,  keepColnames = TRUE)
  }
  outputFile <- str_replace(inputFile, extension, paste0("_Randomized", ".csv"))
  fullTab <- enframe(rownames(genotypesR), name = NULL, value = "ID") %>%
    bind_cols(as_tibble(genotypesR)) %>%
    bind_cols(as_tibble(phenotypesR))
  write_csv(fullTab, outputFile)
  fullTab
}
