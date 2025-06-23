test_that("extractUniqueRows works on duplicates", {
  X <- matrix(c(1,1,0,
                1,1,0,
                0,0,1), ncol = 3, byrow = TRUE)
  res <- extractUniqueRows(X)
  expect_equal(nrow(res[[1]]), 2)           # two unique rows
  expect_equal(length(res[[2]]), 2)         # grouping information returned
})

test_that("extractUniqueRows works on edge cases", {
  for (Q in 0:1) {
    X <- matrix(0, nrow = Q, ncol = 3, byrow = TRUE)
    res <- extractUniqueRows(X)
    expect_equal(nrow(res[[1]]), Q)           # correct number of unique rows
    expect_equal(length(res[[2]]), 1)         # correct grouping information
  }
})

test_that("deleteZeroColumns removes all-zero cols", {
  X <- matrix(c(1,0,0,
                0,0,1), ncol = 3, byrow = TRUE)
  out <- deleteZeroColumns(X)
  expect_equal(ncol(out), 2)
  expect_false(any(colSums(out) == 0))
})

test_that("findSubsetRows identifies logical subsets", {
  M <- matrix(c(TRUE, FALSE,
                TRUE, TRUE,
                FALSE, FALSE), ncol = 2, byrow = TRUE)
  v <- c(TRUE, FALSE)
  expect_equal(findSubsetRows(M, v), c(TRUE, FALSE, TRUE))
})

test_that("pruneSupersets removes logical supersets", {
  M <- matrix(c(TRUE, FALSE, TRUE,
                FALSE, TRUE, FALSE,
                FALSE, FALSE, TRUE), ncol = 3, byrow = TRUE)
  expect_equal(nrow(pruneSupersets(M)), 2)
})

test_that("parseFormula + applyFormula roundtrip for CNF/DNF", {
  pheno <- tibble::tibble(
    A = c(TRUE,  TRUE, FALSE, FALSE),
    B = c(TRUE, FALSE, TRUE,  FALSE),
    C = c(FALSE, TRUE, TRUE,  FALSE)
  )
  cnf <- "(A OR B) AND (C)"
  dnf <- "(A AND C) OR (B AND C)"
  for (f in c(cnf, dnf)) {
    type <- ifelse(identical(f, cnf), "CNF", "DNF")
    mat  <- parseFormula(f, type)
    vec  <- applyFormula(mat, pheno, type)
    expect_type(vec, "logical")
    expect_equal(nrow(vec), nrow(pheno))
  }
})

test_that("applyFormula typo fixed", {
  pheno <- tibble::tibble(A = c(TRUE), B = c(TRUE))
  mat <- matrix(c(1,1), nrow = 2, dimnames = list(c("A","B"), "1"))
  expect_no_error(applyFormula(mat, pheno, "CNF"))
})
