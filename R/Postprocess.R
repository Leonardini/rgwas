### Parsing function for the output from CPLEX
parseLines = function(resultLines, prefix1 = PREFIX1, prefix2 = PREFIX2, optFactor = -1L, ignoreNegatives = !is.null(prefix2)) {
  timeLines <- str_detect(resultLines, "^Solution time")
  stopifnot(sum(timeLines) == 1)
  timeLine <- resultLines[timeLines]
  myTime <- parse_double(str_match(timeLine, "([0-9]+[\\.]*[0-9]+)")[1,1])
  status <- INFEASIBLE
  failLines <- str_detect(resultLines, "^MIP - Integer infeasible") | str_detect(resultLines, "^MIP - Integer unbounded")
  if (any(failLines)) {
    print("The problem was found to be infeasible!")
    return(list(usedVars = NA, optimum = NA, status = status, time = myTime, gap = NA))
  }
  memoryLines <- str_detect(resultLines, "^MIP - Memory limit exceeded, no integer solution")
  if (any(memoryLines)) {
    print("The problem used up too much memory!")
    return(list(usedVars = NA, optimum = NA, status = MEMORYOUT, time = myTime, gap = NA))
  }
  feasibleLines <- str_detect(resultLines, "^MIP - Solution limit exceeded, integer feasible")
  if (any(feasibleLines)) {
    print("Found a feasible solution and stopped!")
    status <- FEASIBLE
    optLine <- resultLines[feasibleLines]
    optIndex <- which(feasibleLines)[1]
  }
  timeOutLines <- str_detect(resultLines, "^MIP - Time limit exceeded")
  if (any(timeOutLines)) {
    print("The problem was not solved within the time limit!")
    status <- TIMEOUT
    optLine <- resultLines[timeOutLines]
    optIndex <- which(timeOutLines)[1]
  }
  optLines <- str_detect(resultLines, "^MIP - Integer optimal")
  if (any(optLines)) {
    status <- OPTIMAL
    optLine <- resultLines[optLines]
  }
  if (status == FEASIBLE || status == TIMEOUT) {
    gapLine <- resultLines[optIndex + 1]
    myGap <- parse_double(str_sub(str_extract(gapLine, "[0-9.]*\\%"), end = -2))
  }
  stopifnot(length(optLine) == 1)
  optPos  <- str_locate(optLine, "=")
  optimum <- parse_double(str_sub(optLine, start = optPos[1,1] + 1)) * optFactor
  usedVars <- str_subset(resultLines, paste0("^", prefix1, "[0-9]+"))
  if (ignoreNegatives) {
    usedVars %<>% 
      str_subset("-0\\.", negate = TRUE)
  }
  if (!is.null(prefix2)) {
    myGroups <- str_match(usedVars, paste0(prefix1, "([0-9]+)", prefix2, "([0-9]+)","[ ]*","([0-9]+[\\.]*[0-9]+)"))
    usedVars <- tibble(type = parse_integer(myGroups[,2]), iter = parse_integer(myGroups[,3]), value = parse_double(myGroups[,4]))
  } else {
    myGroups <- str_match(usedVars, paste0(prefix1, "([0-9]+)", "[ ]*","([-]{0,1}[0-9]+[\\.]*[0-9]+)"))
    usedVars <- tibble(type = parse_integer(myGroups[,2]), iter = 1, value = parse_double(myGroups[,3]))
  }
  if (nrow(usedVars) > 0) {
    usedVars %<>%
      filter(!near(value, 0, tol = TOL)) %>%
      mutate(value = value/max(abs(value)))
    usedVars <- as.matrix(slam::simple_triplet_matrix(i = usedVars$type, j = usedVars$iter, v = usedVars$value))
  } else {
    usedVars <- matrix(0, 0, 0)
  }
  gap <- ifelse(status == OPTIMAL, 0, myGap)
  output <- list(usedVars = usedVars, optimum = optimum, status = status, time = myTime, gap = gap)
  output
}

### Solution extracting function for the output from CPLEX via Rcplex2
extractSolutionAPI = function(solution, prefix1 = PREFIX1, prefix2 = PREFIX2, optFactor = -1L, ignoreNegatives = !is.null(prefix2)) {
  usedVars = solution$xopt
  if (ignoreNegatives) {
    usedVars %<>%
      magrittr::extract(. >= 0)
  }
  usedVars %<>% 
    magrittr::extract(!near(., 0, tol = TOL)) %>%
    enframe() %>%
    filter(str_starts(name, prefix1))
  if (!is.null(prefix2)) {
    usedVars %<>%
      mutate(sub = str_remove_all(name, paste0("^", prefix1)) %>% str_split_fixed(prefix2, n = 2), type = as.integer(sub[,1]), iter = as.integer(sub[,2]))
  } else {
    usedVars %<>%
      mutate(type = str_remove_all(name, paste0("^", prefix1)) %>% as.integer(), iter = 1L)
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

### This function postprocesses a given solution of an ILP/MIQCP to put it into a convenient format: the phenotype and the clause
postprocessSolution = function(origPhenotypes, solutionMatrix, type, complement = 0) {
  if (type != "LCP") {
    if (any(is.na(solutionMatrix)) || all(solutionMatrix == 0)) {
      return(list(solution = matrix(ifelse(type == "CNF", TRUE, FALSE), nrow(origPhenotypes), 1), clause = ""))
    }
    mode(solutionMatrix) = "logical"
    solutionMatrix %<>%
      deleteZeroColumns %>%
      extractMinimalColumns(keepDuplicates = FALSE)
    if ((complement < 0 || complement > 1) && type %in% c("CNF", "DNF")) { ### negated phenotypes; look for ways to simplify!
      rown    <- rownames(solutionMatrix)
      singletons <- which(colSums(solutionMatrix) == 1)
      if (length(singletons) > 0) {
        singletonRows <- rown[which(solutionMatrix[, singletons, drop = FALSE], arr.ind = TRUE)[, "row"]]
        negSingletonRows <- str_c("NOT_", singletonRows) %>% str_remove_all("NOT_NOT_")
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
  } else {
    if (any(is.na(solutionMatrix)) || all(near(solutionMatrix, 0))) {
      return(list(solution = matrix(0, nrow(origPhenotypes), 1), clause = ""))
    }
    goodPos <- which(!(near(solutionMatrix, 0)))
    curSolution <- origPhenotypes[, rownames(solutionMatrix), drop = FALSE] %*% solutionMatrix
    solutionMatrix <- solutionMatrix[goodPos, 1, drop = FALSE]
    solOrder <- order(abs(solutionMatrix[,1]), decreasing = TRUE)
    clause <- paste0("(", solutionMatrix[solOrder, 1], ") * ", rownames(solutionMatrix)[solOrder], collapse = " + ")
  }
  output <- list(solution = curSolution, clause = clause)
  output
}

### This function computes a number of 2 x 2 statistics based on a tibble containing at a minimum the columns TP, FP, FN and TN
### If pvals = TRUE, additionally computes the Fisher exact and chi-squared p-values (in log-space).
computeStats = function(Tibble, pvals = TRUE) {
  Tibble %<>%
    mutate(Pos = TP + FN, Neg = TN + FP) %>%
    mutate(PredPos = TP + FP, PredNeg = TN + FN) %>%
    mutate(Correct = TP + TN, Incorrect = FP + FN) %>%
    mutate(Total = TP + FP + TN + FN) %>%
    mutate(TPR = TP/Pos, TNR = TN/Neg) %>%
    mutate(PPV = TP/PredPos, NPV = TN/PredNeg) %>%
    mutate(Accuracy = Correct/Total) %>%
    mutate(OR = (TP / FP) * (TN / FN)) %>%
    mutate(F1 = 1/((1/TPR + 1/PPV)/2)) %>%
    mutate(Jaccard = TP/(TP + FP + FN)) %>%
    mutate(MCC = (TP - Pos / Total * PredPos) / (sqrt(Pos / Total * Neg) * sqrt(PredPos / Total * PredNeg)))
  if (pvals) {
    Tibble %<>% mutate(
      LogChiSquaredP  = pchisq(Total * MCC^2, df = 1, lower.tail = FALSE, log.p = TRUE),
      LogFisherExactP = phyper(TP - 1, Pos, Neg, PredPos, lower.tail = FALSE, log.p = TRUE))
  }
  Tibble
}

### This function prepares the final output from a sequence of ILPs, one per SNP, in the form of 2 tables: summary and phenotype 
prepareOutput = function(output, complexPhenotypes, SNP, ID, type, numSegments, objective, singleRes, boundValues, ns, MH, 
                         outputAssociations = FALSE) {
  myTab <- matrix(unlist(output), nrow = length(output), dimnames = list(names(output), names(output[[1]])), byrow = TRUE)
  myTab %<>%
    as_tibble %>%
    # mutate(simpleFormula = map_chr(formula, ~{simplifyFormula(., type = type)}), .after = formula) %>%
    mutate(SNP = SNP) %>%
    select(SNP, everything())
  if (type != "LCP") {
    myTab %<>%
      mutate(numSegments = numSegments) %>%
      mutate_at(c("TP", "FP", "TN", "FN"), as.integer) %>%
      computeStats()
    if (!all(is.na(MH))) {
      myTab %<>%
        mutate(LogCMHExactP = MH)
    }
  }
  myTab %<>%
    mutate_at(c("time", "gap", "opt"), parse_double)
  if (objective == "covariance") {
    myTab %<>%
      mutate_at("opt", ~divide_by(., ns))
  }
  myTab %<>%
    bind_cols(tibble(bestSingle = singleRes$phenotype, bestSingleStat = singleRes$value, computedBounds = boundValues))
  if (type == "LCP") {
    myTab %<>% 
      bind_cols(tibble(bestPValue = singleRes$significance))
  }
  if (!is.null(complexPhenotypes)) {
    complexPhenotypes %<>%
      as_tibble %>%
      mutate_all(as.numeric) %>%
      mutate(ID = ID) %>%
      select(ID, everything())
  }
  fullOutput <- list(summary = myTab, phenotype = complexPhenotypes)
  if (outputAssociations) { fullOutput = c(fullOutput, list(associations = singleRes$allValues)) }
  fullOutput
}

simplifyFormula = function(initFormula, type = "CNF") {
  if (initFormula == "" || type == "LCP") {
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
        str_split("\\+") %>% 
        extract2(1) %>% 
        str_trim() %>% 
        str_split("\\*") %>% 
        lapply(function(x) { 
          x %>% 
            str_trim() %>% 
            str_replace_all("^~", "NOT_") %>% 
            sort() })
      DNFform = paste0(sapply(DNFclauses, function(x) {paste0("(", paste(sort(x), collapse = " AND "), ")")}), collapse = " OR ")
    } else {
      DNFformula = matrix(str_replace_all(DNFformula, "^~", "NOT_"), ncol = ncol(DNFformula))
      DNFform = paste0(apply(DNFformula, 1, function(x) {paste0("(", paste(sort(x), collapse = " AND "), ")")}), collapse = " OR ")
    }
    return(DNFform)
  }
}

prepareFormula = function(initFormula, type = "CNF") {
  clauses = str_split(initFormula, ifelse(type == "CNF", "AND", "OR"))[[1]] %>%
    str_trim() %>%
    str_remove("^\\(") %>%
    str_remove("\\)$")
  fullClauses = clauses %>%
    str_split(ifelse(type == "CNF", "OR", "AND")) %>%
    lapply(function(x) { 
      x %>% 
        str_trim() %>% 
        str_replace_all("^NOT_", "~") })
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
  clauseSizes = sapply(preparedFormula, length)
  allComboInds = combinat::hcube(clauseSizes)
  M = nrow(allComboInds)
  N = ncol(allComboInds)
  allCombos = matrix(NA, M, N)
  for (ind in 1:N) {
    curClause = preparedFormula[[ind]]
    curCol = allComboInds[, ind]
    allCombos[, ind] = curClause[curCol]
  }
  uPheno = unique(as.vector(str_remove_all(allCombos, "~")))
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
