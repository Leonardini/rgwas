#' Run `reverseGWAS` on a dataset.
#'
#' Main driver for a multi-genotype multi-phenotype association with a single type, number of clauses and clause size.
#'
#' @param inputFile The absolute path to an input file. This can be a CSV or TSV file, either uncompressed or compressed with gz. NB: the last integer before the file extension must contain the number of genotype columns!
#' @param type The type of phenotype combinations to search over; the allowed values are "CNF", or conjunctive normal form (the default) and "DNF", or disjunctive normal form
#' @param objective The association function between genotype and complex phenotype to optimise; the allowed values are "agreement" (the default) and "covariance"
#' @param K,L The maximum number of clauses (K) and the maximum size of each clause (L) to use in searching for the combination phenotype. NB: high values may lead to overfitting to the training (e.g. discovery) data!
#' @param complement Integer between -2 and 1, or a value coercible to one. If complement =  1 or 2, each genotype's negated form is added in preprocessing. If complement = -1 or 2, each phenotype's negated form is added in preprocessing.
#' @param ext The extension of the uncompressed input file (CSV or TSV); by default, it is inferred automatically from the input filename.
#' @param index Offset value to be used in distinguishing between the iterations; the default is `1`.
#' @param extremeValue The target value of the association function (by default, the Fisher exact test p-value). The default is log(5e-8), i.e. genome-wide significance. NB: it must be in logarithmic space (i.e. non-positive)!
#' @param outputSummary Create a file with the summary statistics for the combination phenotype? The default is `TRUE`.
#' @param outputPhenotype Create a file with the best combination phenotype's value for each patient? The default is `FALSE`.
#' @param outputAssociations Create a file with the single best phenotype associated to each of the genotypes? The default is `FALSE`.
#' @param saveOutput Save the output as an RDS file? The default is `FALSE`.
#'
#' @returns A list containing a summary of the association statistics (one row per genotype) and the best combination phenotypes (one column per genotype).
#'
#' @examples # X <- mainDriver(inputFile = "TestInputN5000P5_3.csv", extremeValue = log(1e-3))
#' @export
mainDriver = function(inputFile, type = "CNF", objective = "agreement", K = 3, L = 3, complement = FALSE,
                      ext = stringr::str_sub(inputFile, -4), index = 1L, extremeValue = log(MAX_P),
                      outputSummary = TRUE, outputPhenotype = FALSE, outputAssociations = FALSE, saveOutput = FALSE) {
  miniFile <- inputFile
  if (stringr::str_sub(inputFile, -3) == ".gz") {
    ext <- stringr::str_sub(inputFile, -7, -4)
    miniFile <- stringr::str_sub(inputFile, 1, -4)
  }
  readFunction <- ifelse(ext == '.tsv', readr::read_tsv, readr::read_csv)
  inputTab <- readFunction(inputFile, col_types = readr::cols(ID = "c", .default = "l"))
  numSNPs <- stringr::str_extract(miniFile, paste0("([0-9]+)", ext)) %>%
    stringr::str_remove(ext) %>%
    readr::parse_integer()
  numPheno <- (ncol(inputTab) - 1 - numSNPs) * (1 + (complement %in% c(-1, 2)))
  numCombos <- choose(choose(numPheno, L), K)
  myP <- extremeValue %>% as.character()
  myF <- paste0("F", Inf)
  baseFilename <- stringr::str_remove(miniFile, ext)
  stage <- ifelse(OPTIMIZATION, "Optimization", "Triage")
  specString <- paste(c("", type, paste0("K", K), paste0("L", L), paste0("P", myP), myF, objective, stage, "Results"), collapse = "_")
  outFileSum <- stringr::str_replace(miniFile, ext, stringr::str_c(specString, '_Summary.csv'))
  if (!(file.exists(outFileSum))) {
    output <- createAndSolveILPs(inputTab, numSNPs = numSNPs, type = type, objective = objective, K = K, L = L, index = index,
                                 extremeValue = extremeValue, complement = complement, baseFilename = baseFilename, outputAssociations = FALSE, strata = TRUE)
    if (outputSummary && !is.null(output$summary)) {
      readr::write_csv(output$summary,      stringr::str_replace(miniFile, ext, stringr::str_c(specString,      '_Summary.csv')))
    }
    if (outputPhenotype && !is.null(output$phenotype))   {
      readr::write_csv(output$phenotype,    stringr::str_replace(miniFile, ext, stringr::str_c(specString,    '_Phenotype.csv')))
    }
    if (outputAssociations && !is.null(output$associations)) {
      readr::write_csv(output$associations, stringr::str_replace(miniFile, ext, stringr::str_c(specString, '_Associations.csv')))
    }
    if (saveOutput)    {
      saveRDS(object = output, file = stringr::str_replace(miniFile, ext, stringr::str_c(specString,              '.RDS')))
    }
    output
  } else {
    print(paste("Not executing the required computation because the output file", outFileSum, "already exists."))
    print("To force a recomputation, please delete the file and try again.")
  }
}

optimizingDriver = function(inputFile, startPValue = log(1e-10), factor = 2, type = "CNF", complement = -1, K = MY_K, L = MY_L, startIndex = 0) {
  stopifnot(startPValue <= 0)
  curPValue <- startPValue
  ind <- 0
  curOutput <- mainDriver(inputFile, extremeValue = curPValue, index = startIndex * THOUSAND + ind, complement = complement,
                  type = type, K = K, L = L, outputSummary = FALSE, outputPhenotype = FALSE, outputAssociations = FALSE)
  prevOutput <- curOutput
  while (curOutput$summary$status[1] %in% c("feasible", "optimal")) {
    prevOutput <- curOutput
    newPValue <- prevOutput$summary$LogFisherExactP[1]
    if (is.na(newPValue)) {
      print("Warning: The p-value was NA, so using the previous value times factor as a stopgap measure!")
      newPValue <- curPValue
    }
    curPValue <- newPValue - log(factor)
    ind <- ind + 1
    curOutput <- mainDriver(inputFile, extremeValue = curPValue, index = startIndex * THOUSAND + ind, complement = complement,
                    type = type, K = K, L = L, outputSummary = FALSE, outputPhenotype = FALSE, outputAssociations = FALSE)
    if (curOutput$summary$status[1] %in% c("timeout", "memoryout")) {
      print(paste("Warning: the optimization is terminated prematurely due to insufficient CPLEX resources (time or memory)"))
      break
    }
  }
  print(paste("The best feasible p-value to within a factor of", factor, "is", exp(prevOutput$summary$LogFisherExactP[1])))
  return(prevOutput)
}
