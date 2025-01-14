#' Extract solution from the output of CPLEX via Rcplex
#' @noRd
extractSolution = function(solution, prefix1 = PREFIX1, prefix2 = PREFIX2, optFactor = -1L, ignoreNegatives = !is.null(prefix2)) {
  usedVars = solution$xopt
  if (ignoreNegatives) {
    usedVars <- usedVars %>%
      magrittr::extract(. >= 0)
  }
  usedVars <- usedVars %>%
    magrittr::extract(!dplyr::near(., 0, tol = TOL)) %>%
    tibble::enframe() %>%
    dplyr::filter(stringr::str_starts(name, prefix1))
  if (!is.null(prefix2)) {
    usedVars <- usedVars %>%
      dplyr::mutate(sub = stringr::str_remove_all(name, paste0("^", prefix1)) %>%
                      stringr::str_split_fixed(prefix2, n = 2),
                    type = as.integer(sub[,1]),
                    iter = as.integer(sub[,2]))
  } else {
    usedVars <- usedVars %>%
      dplyr::mutate(type = stringr::str_remove_all(name, paste0("^", prefix1)) %>%
                      as.integer(),
                    iter = 1L)
  }
  if (nrow(usedVars) > 0) {
    usedVars <- as.matrix(slam::simple_triplet_matrix(i = usedVars$type, j = usedVars$iter, v = usedVars$value))
  } else {
    usedVars <- matrix(0, 0, 0)
  }
  optimum = solution$obj * optFactor
  statusCode = solution$status
  status = HOPELESS
  if (statusCode %in% c(101, 102)) {
    status = OPTIMAL
  }
  if (statusCode %in% c(103, 109, 110, 114, 115, 118, 119)) {
    status = INFEASIBLE
  }
  if (statusCode %in% c(104, 113, 120)) {
    status = FEASIBLE
  }
  if (statusCode %in% c(107, 108, 131, 132)) {
    status = TIMEOUT
  }
  if (statusCode %in% c(105, 106, 111, 112, 116, 117)) {
    status = MEMORYOUT
  }
  output <- list(usedVars = usedVars, optimum = optimum, status = status, time = NA, gap = NA)
  output
}

#' Postprocess a given solution of an ILP to put it into a convenient format: the phenotype and the clause
#' @noRd
postprocessSolution = function(origPhenotypes, solutionMatrix, type, complement = 0) {
  if (any(is.na(solutionMatrix)) || all(solutionMatrix == 0)) {
    return(list(solution = matrix(ifelse(type == "CNF", TRUE, FALSE), nrow(origPhenotypes), 1), clause = ""))
  }
  mode(solutionMatrix) = "logical"
  solutionMatrix <- solutionMatrix %>%
    deleteZeroColumns %>%
    extractMinimalColumns
  if ((complement < 0 || complement > 1) && type %in% c("CNF", "DNF")) { # negated phenotypes
    rown    <- rownames(solutionMatrix)
    singletons <- which(colSums(solutionMatrix) == 1)
    if (length(singletons) > 0) {
      singletonRows <- rown[which(solutionMatrix[, singletons, drop = FALSE], arr.ind = TRUE)[, "row"]]
      negSingletonRows <- stringr::str_c("NOT_", singletonRows) %>% stringr::str_remove_all("NOT_NOT_")
      negFound <- intersect(negSingletonRows, rown)
      solutionMatrix[negFound, ] <- 0
    }
    if (all(solutionMatrix == 0)) {
      return(list(solution = matrix(ifelse(type == "CNF", TRUE, FALSE), nrow(origPhenotypes), 1), clause = ""))
    }
  }
  usedInds <- which(solutionMatrix, arr.ind = TRUE)
  splitInds <- split(usedInds[,'row'], usedInds[,'col'])
  subPhenotypes <- sapply(splitInds,
                          function(x) {
                            apply(origPhenotypes[, rownames(solutionMatrix)[x], drop = FALSE], 1, ifelse(type == "CNF", any, all))
                          })
  curSolution <- apply(subPhenotypes, 1, ifelse(type == "CNF", all, any))
  clauses <- sapply(splitInds,
                    function(x) {paste(rownames(solutionMatrix)[x], collapse = ifelse(type == "CNF", " OR ", " AND "))})
  clause <- paste0("(", paste(clauses, collapse = paste(")", ifelse(type == "CNF", " AND ", " OR "), "(")), ")")
  output <- list(solution = curSolution, clause = clause)
  output
}

#' Compute a number of 2 x 2 statistics based on a tibble containing at a minimum the columns TP, FP, FN and TN
#' If pvals = TRUE, additionally compute the Fisher exact and chi-squared p-values (in log-space).
#' @noRd
computeStats = function(Tibble, pvals = TRUE) {
  Tibble <- Tibble %>%
    dplyr::mutate(Pos = TP + FN, Neg = TN + FP) %>%
    dplyr::mutate(PredPos = TP + FP, PredNeg = TN + FN) %>%
    dplyr::mutate(Correct = TP + TN, Incorrect = FP + FN) %>%
    dplyr::mutate(Total = TP + FP + TN + FN) %>%
    dplyr::mutate(TPR = TP/Pos, TNR = TN/Neg) %>%
    dplyr::mutate(PPV = TP/PredPos, NPV = TN/PredNeg) %>%
    dplyr::mutate(Accuracy = Correct/Total) %>%
    dplyr::mutate(OR = (TP / FP) * (TN / FN)) %>%
    dplyr::mutate(F1 = 1/((1/TPR + 1/PPV)/2)) %>%
    dplyr::mutate(Jaccard = TP/(TP + FP + FN)) %>%
    dplyr::mutate(MCC = (TP - Pos / Total * PredPos) / (sqrt(Pos / Total * Neg) * sqrt(PredPos / Total * PredNeg)))
  if (pvals) {
    Tibble <- Tibble %>% dplyr::mutate(
      LogChiSquaredP  = pchisq(Total * MCC^2, df = 1, lower.tail = FALSE, log.p = TRUE),
      LogFisherExactP = phyper(TP - 1, Pos, Neg, PredPos, lower.tail = FALSE, log.p = TRUE))
  }
  Tibble
}

#' Prepare the final output from a sequence of ILPs, one per SNP, in the form of 2 tables: summary and phenotype
#' @noRd
prepareOutput = function(output, complexPhenotypes, SNP, ID, type, numSegments, objective, singleRes, boundValues, ns, MH,
                         outputAssociations = FALSE) {
  myTab <- matrix(unlist(output), nrow = length(output), dimnames = list(names(output), names(output[[1]])), byrow = TRUE)
  myTab <- myTab %>%
    tibble::as_tibble() %>%
    # dplyr::mutate(simpleFormula = map_chr(formula, ~{simplifyFormula(., type = type)}), .after = formula) %>%
    dplyr::mutate(SNP = SNP) %>%
    dplyr::select(SNP, everything())
  myTab <- myTab %>%
    dplyr::mutate(numSegments = numSegments) %>%
    dplyr::mutate_at(c("TP", "FP", "TN", "FN"), as.integer) %>%
    computeStats()
  if (!all(is.na(MH))) {
    myTab <- myTab %>%
      dplyr::mutate(LogCMHExactP = MH)
  }
  myTab <- myTab %>%
    dplyr::mutate_at(c("time", "gap", "opt"), readr::parse_double)
  if (objective == "covariance") {
    myTab <- myTab %>%
      dplyr::mutate_at("opt", ~magrittr::divide_by(., ns))
  }
  myTab <- myTab %>%
    dplyr::bind_cols(tibble::tibble(bestSingle = singleRes$phenotype, bestSingleStat = singleRes$value, computedBounds = boundValues))
  if (!is.null(complexPhenotypes)) {
    complexPhenotypes <- complexPhenotypes %>%
      tibble::as_tibble() %>%
      dplyr::mutate_all(as.numeric) %>%
      dplyr::mutate(ID = ID) %>%
      dplyr::select(ID, everything())
  }
  fullOutput <- list(summary = myTab, phenotype = complexPhenotypes)
  if (outputAssociations) { fullOutput = c(fullOutput, list(associations = singleRes$allValues)) }
  fullOutput
}

simplifyFormula = function(initFormula, type = "CNF") {
  if (initFormula == "") {
    return(initFormula)
  } else {
    prepFormula = prepareFormula(initFormula, type = type)
    if (type == "CNF") {
      DNFformula = makeDNF(prepFormula)
    } else {
      DNFformula = prepFormula
    }
    if (ncol(DNFformula) > 1) {
      QCAInput = paste0(apply(DNFformula, 1, function(x) {paste(x, collapse = "*")}), collapse = " + ")
      DNFformula = admisc::simplify(QCAInput)
      DNFclauses = DNFformula %>%
        stringr::str_split("\\+") %>%
        magrittr::extract2(1) %>%
        stringr::str_trim() %>%
        stringr::str_split("\\*") %>%
        lapply(function(x) {
          x %>%
            stringr::str_trim() %>%
            stringr::str_replace_all("^~", "NOT_") %>%
            sort() })
      DNFform = paste0(sapply(DNFclauses, function(x) {paste0("(", paste(sort(x), collapse = " AND "), ")")}), collapse = " OR ")
    } else {
      DNFformula = matrix(stringr::str_replace_all(DNFformula, "^~", "NOT_"), ncol = ncol(DNFformula))
      DNFform = paste0(apply(DNFformula, 1, function(x) {paste0("(", paste(sort(x), collapse = " AND "), ")")}), collapse = " OR ")
    }
    return(DNFform)
  }
}

prepareFormula = function(initFormula, type = "CNF") {
  clauses = stringr::str_split(initFormula, ifelse(type == "CNF", "AND", "OR"))[[1]] %>%
    stringr::str_trim() %>%
    stringr::str_remove("^\\(") %>%
    stringr::str_remove("\\)$")
  fullClauses = clauses %>%
    stringr::str_split(ifelse(type == "CNF", "OR", "AND")) %>%
    lapply(function(x) {
      x %>%
        stringr::str_trim() %>%
        stringr::str_replace_all("^NOT_", "~") })
  if (type == "DNF") {
    L = length(fullClauses)
    Lens = sapply(fullClauses, length)
    maxL = max(Lens)
    Mat = matrix("", L, maxL)
    for (ind in 1:L) {
      Mat[ind, ] = fullClauses[[ind]][1]
      Mat[ind, 1:Lens[ind]] = fullClauses[[ind]]
    }
    fullClauses = Mat
  }
  fullClauses
}

makeDNF = function(preparedFormula) {
  allCombos = do.call(expand.grid, preparedFormula)
  uPheno    = unique(as.vector(stringr::str_remove_all(unlist(preparedFormula), "~")))
  badRows = c()
  for (pheno in uPheno) {
    curIndsPos = which(allCombos == pheno, arr.ind = TRUE)[, "row"]
    curIndsNeg = which(allCombos == paste0("~", pheno), arr.ind = TRUE)[, "row"]
    badRows = c(badRows, intersect(curIndsPos, curIndsNeg))
  }
  if (length(badRows) > 0) {
    allCombos = allCombos[-badRows, , drop = FALSE]
  }
  allCombos
}
