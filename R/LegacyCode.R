#' Extract the unique rows of a table; if repeats = TRUE, also returns the groups corresponding to each unique row
#' Partly based on the solution at r.789695.n4.nabble.com/Count-unique-rows-columns-in-a-matrix-td844731.html
#' @noRd
extractUniqueRowsOld = function(Table, repeats = TRUE) {
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

#' Return the rows in inputMatrix which are not proper supersets of other rows in the same matrix
#' Note: assumes, without checking, that inputMatrix does not have any comparable rows (in particular, any duplicates)!
#' @noRd
pruneSupersetsOld = function(inputMatrix) {
  M1 <- nrow(inputMatrix)
  goodRows <- rep(TRUE, M1)
  for (ind in 1:M1) {
    eliminate <- findSubsetRows(inputMatrix, inputMatrix[ind,])
    eliminate[ind] <- FALSE
    goodRows[ind]  <- !(any(eliminate))
  }
  outputMatrix <- inputMatrix[goodRows, , drop = FALSE]
  outputMatrix
}
