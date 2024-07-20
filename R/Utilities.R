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
    repList <- lapply(1:N, function(x) { myOrder[lowerBounds[x]:upperBounds[x]] })
    output <- c(output, list(repList))
  }
  output
}

deleteRepeatedRows = function(inputMatrix) {
  uniqueRows <- extractUniqueRows(inputMatrix, repeats = FALSE)[[1]]
  uniqueRows
}

deleteRepeatedColumns = function(inputMatrix) {
  uniqueCols <- t(deleteRepeatedRows(t(inputMatrix)))
  uniqueCols
}

deleteZeroColumns = function(inputMatrix) {
  goodInds <- which(colSums(inputMatrix) > 0)
  inputMatrix %<>% 
    subset(select = goodInds)
  inputMatrix
}

### Removes all rows in matrix1 that also occur in matrix2; assumes, without checking, that matrix2 has no duplicate
removeOccurrences = function(matrix1, matrix2) {
  n2 <- nrow(matrix2)
  fullMatrix <- rbind(matrix2, matrix1)
  fullMatrix <- extractUniqueRows(fullMatrix, repeats = TRUE)
  repeats <- fullMatrix[[2]]
  goodGroups <- findIndicesInGroups(1:n2, repeats)
  outputMatrix <- fullMatrix[[1]][!goodGroups, , drop = FALSE]
  outputMatrix
}

findSubsetRows = function(inputMatrix, inputVector, skipColumns = NULL) {
  M <- nrow(inputMatrix)
  N <- ncol(inputMatrix)
  L <- length(inputVector)
  stopifnot(L == N)
  goodRows <- rep(TRUE, M)
  goodPos <- setdiff(which(!inputVector), skipColumns)
  for (ind in goodPos) {
    goodRows <- goodRows & (!inputMatrix[, ind])
  }
  goodRows
}

findSupersetRows = function(inputMatrix, inputVector, skipColumns = NULL) {
  M <- nrow(inputMatrix)
  N <- ncol(inputMatrix)
  L <- length(inputVector)
  stopifnot(L == N)
  goodRows <- rep(TRUE, M)
  goodPos <- setdiff(which(inputVector), skipColumns)
  for (ind in goodPos) {
    goodRows <- goodRows & inputMatrix[, ind]
  }
  goodRows
}

### Return the rows in matrix1 which are not supersets of those in matrix2
### If self = TRUE, assumes that matrix1 = matrix2 and returns the rows that are not PROPER supersets of other rows
### Note: assumes, without checking, that matrix2 does not have any comparable rows (in particular, any duplicates)!
pruneSupersets = function(matrix1, matrix2, self = FALSE) {
  if (!self) {
    matrix1 <- removeOccurrences(matrix1, matrix2)
  }
  M1 <- nrow(matrix1)
  M2 <- nrow(matrix2)
  stopifnot(ncol(matrix1) == ncol(matrix2))
  goodRows <- rep(TRUE, M1)
  mySize <- min(M1, M2)
  if (M1 <= M2) {
    # print("The first matrix is smaller or the same")
    for (ind in 1:M1) {
      eliminate <- findSubsetRows(matrix2, matrix1[ind,])
      if (self) {
        eliminate[ind] <- FALSE
      }
      goodRows[ind] <- !(any(eliminate))
    }
  }
  else {
    # print("The second matrix is smaller")
    for (ind in 1:M2) {
      eliminate <- findSupersetRows(matrix1, matrix2[ind,])
      goodRows <- goodRows & (!eliminate)
    }
  }
  outputMatrix <- matrix1[goodRows, , drop = FALSE]
  outputMatrix
}

extractMinimalRows = function(inputMatrix, keepDuplicates = FALSE) {
  if (!keepDuplicates) {
    inputMatrix <- deleteRepeatedRows(inputMatrix)
  }
  outputMatrix <- pruneSupersets(inputMatrix, inputMatrix, self = TRUE)
  outputMatrix
}

extractMinimalColumns = function(inputMatrix, keepDuplicates = FALSE) {
  outputMatrix <- t(extractMinimalRows(t(inputMatrix), keepDuplicates = keepDuplicates))
  outputMatrix
}

### Auxiliary function to create name pairs, separated by "I"
pasteI <- function(x, y) {
  output <- paste(y, x, sep = "I")
  output
}
