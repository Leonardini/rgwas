#' Extract the unique rows of a table; if repeats = TRUE, also returns the groups corresponding to each unique row
#' Partly based on the solution at r.789695.n4.nabble.com/Count-unique-rows-columns-in-a-matrix-td844731.html
#' @noRd
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

#' Delete any repeated rows of the input matrix
#' @noRd
deleteRepeatedRows = function(inputMatrix) {
  uniqueRows <- extractUniqueRows(inputMatrix, repeats = FALSE)[[1]]
  uniqueRows
}

#' Delete any all-zero columns of a logical input matrix
#' @noRd
deleteZeroColumns = function(inputMatrix) {
  goodInds <- which(colSums(inputMatrix) > 0)
  inputMatrix <- inputMatrix %>%
    subset(select = goodInds)
  inputMatrix
}

#' Find any rows of the input matrix which are subsets of the input vector; both are assumed to be logical
#' @noRd
findSubsetRows = function(inputMatrix, inputVector) {
  M <- nrow(inputMatrix)
  stopifnot(length(inputVector) == ncol(inputMatrix))
  goodRows <- rep(TRUE, M)
  for (ind in which(!inputVector)) {
    goodRows <- goodRows & (!inputMatrix[, ind])
  }
  goodRows
}

#' Return the rows in inputMatrix which are not proper supersets of other rows in the same matrix
#' Note: assumes, without checking, that inputMatrix does not have any comparable rows (in particular, any duplicates)!
#' @noRd
pruneSupersets = function(inputMatrix) {
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

#' Extract the minimal rows of a logical input matrix with respect to set inclusion
#' @noRd
extractMinimalRows = function(inputMatrix) {
  inputMatrix <- deleteRepeatedRows(inputMatrix)
  outputMatrix <- pruneSupersets(inputMatrix)
  outputMatrix
}

#' Extract the minimal columns of a logical input matrix with respect to set inclusion
#' @noRd
extractMinimalColumns = function(inputMatrix) {
  outputMatrix <- t(extractMinimalRows(t(inputMatrix)))
  outputMatrix
}

#' Parse a formula of type "CNF" or "DNF" into a standard matrix representation
#' @noRd
parseFormula = function(Formula, type) {
  outerString   <- ifelse(type == "CNF", " AND ", " OR ")
  innerString   <- ifelse(type == "DNF", " AND ", " OR ")
  splitFormula <- stringr::str_split(Formula, outerString)[[1]]
  fullySplitFormula <- purrr::map(splitFormula, function(x) {
    x %>%
      stringr::str_trim %>%
      stringr::str_remove("^\\(") %>%
      stringr::str_remove("\\)$") %>%
      stringr::str_split(innerString) %>%
      magrittr::extract2(1) })
  allUsedNames <- sort(unique(unlist(fullySplitFormula)))
  M <- length(allUsedNames)
  N <- length(splitFormula)
  parsedFormula <- matrix(0, M, N, dimnames = list(allUsedNames, 1:N))
  for (ind in 1:N) {
    parsedFormula[fullySplitFormula[[ind]], ind] <- 1
  }
  parsedFormula
}

#' Apply a formula in a standard matrix representation (of type "CNF" or "DNF") to a phenotype matrix (assumed to be logical)
#' @noRd
applyFormula = function(parsedFormula, phenotypes, type) {
  stopifnot(all(rownames(parsedFormula) %in% colnames(phenotypes)))
  phenotypes <- phenotypes %>%
    as.data.frame %>%
    dplyr::select(rownames(parsedFormula))
  outerFunction <- ifelse(type == "CNF", magrittr::and, magrittr::or)
  innerFunction <- ifelse(type == "DNF", magrittr::and, magrittr::or)
  allTerms <- map(1:ncol(parsedFormula), ~{purrr::reduce(phenotypes[, (parsedFormula[, .] == 1), drop = FALSE], innerFunction)})
  finalOutput   <- reduce(allTermds, outerFunction)
  if (class(finalOutput) != "matrix") { finalOutput <- matrix(finalOutput, ncol = 1) }
  finalOutput
}
