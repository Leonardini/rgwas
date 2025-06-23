test_that("mainDriver handles tiny synthetic dataset", {
  tf <- tempfile(fileext = "_1.csv")
  toy <- tibble::tibble(
    ID = sprintf("P%02d", 1:4),
    g0 = c(TRUE, FALSE, TRUE, FALSE),
    p0 = c(TRUE, TRUE, FALSE, FALSE)
  )
  readr::write_csv(toy, tf)
  res <- mainDriver(tf, extremeValue = log(1e-2),
                    K = 1, L = 1,
                    outputSummary = FALSE,
                    outputPhenotype = FALSE,
                    saveOutput = FALSE)
  expect_type(res, "list")
  expect_true(all(c("summary","phenotype") %in% names(res)))
})
