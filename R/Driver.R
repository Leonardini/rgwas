#' Main driver
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
