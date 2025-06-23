#' Extract the unique rows of a table; if repeats = TRUE, also returns the groups corresponding to each unique row
#' Based on the suggestion made by ChatGPT, replacing the more complicated original workflow
#' @noRd
extractUniqueRows = function(Table, repeats = TRUE) {
  Q = nrow(Table)
  if (Q <= 1) {
    output <- list(Table)
    if (repeats) {
      output <- c(output, list(Q))
    }
    return(output)
  }
  rown = rownames(Table)
  coln = colnames(Table)
  rownames(Table) = NULL
  colnames(Table) = paste0("C", seq_len(ncol(Table)))
  Table = dplyr::as_tibble(Table) %>%
    dplyr::group_by(!!!rlang::syms(colnames(Table)))
  groups = dplyr::group_rows(Table)
  Table = Table %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    as.matrix() %>%
    magrittr::set_rownames(rown[sapply(groups, dplyr::first)]) %>%
    magrittr::set_colnames(coln)
  output = list(Table)
  if (repeats) {
    output <- c(output, list(groups))
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

#' Return the rows in inputMatrix which are not proper supersets of other rows in the same matrix, following ChatGPT's suggestion
#' Note: assumes, without checking, that inputMatrix does not have any comparable rows (in particular, any duplicates)!
#' @noRd
pruneSupersets = function(inputMatrix) {
  res <- inputMatrix %*% t(!inputMatrix)
  # res[i,j] = n. entries where r_i = TRUE and r_j = FALSE; r_i is a superset of r_j iff res[i,j] > 0 and res[j,i] = 0
  supersets <- (res > 0 & t(res) == 0)
  goodRows <- which(rowSums(supersets) == 0)
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
  outerString   <- ifelse(type == "CNF", "AND", "OR")
  innerString   <- ifelse(type == "DNF", "AND", "OR")
  splitFormula <- stringr::str_split(Formula, outerString)[[1]]
  fullySplitFormula <- purrr::map(splitFormula, function(x) {
    x %>%
      stringr::str_trim() %>%
      stringr::str_remove("^\\(") %>%
      stringr::str_remove("\\)$") %>%
      stringr::str_split(innerString) %>%
      magrittr::extract2(1) %>%
      stringr::str_trim() })
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
  allTerms <- purrr::map(1:ncol(parsedFormula), ~{purrr::reduce(phenotypes[, (parsedFormula[, .] == 1), drop = FALSE], innerFunction)})
  finalOutput   <- purrr::reduce(allTerms, outerFunction)
  if (class(finalOutput) != "matrix") { finalOutput <- matrix(finalOutput, ncol = 1) }
  finalOutput
}
