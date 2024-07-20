library(combinat)
library(tidyverse)

MAX_TOP = 10

bruteForceOptimize = function(inputTab, numSNPs, type = "CNF", K = 3, L = 3, index = 0, complement = 2, 
                              extremeValue = MAX_P, outputAssociations = FALSE) {
  IDs <- inputTab$ID
  myMatrices <- prepareMatrices(inputTab, numSNPs, complement = ifelse(type == "DNF", 2, complement))
  phenotypes <- myMatrices$phenotypes
  genotypes  <- myMatrices$genotypes
  if (type == "DNF") {
    M <- ncol( genotypes)/2
    colnames( genotypes) <- colnames( genotypes)[c(M + (1:M), 1:M)]
    if (complement %in% c(-1, 0)) {  genotypes <-  genotypes[, -(1:M), drop = FALSE] }
    N <- ncol(phenotypes)/2
    colnames(phenotypes) <- colnames(phenotypes)[c(N + (1:N), 1:N)]
    if (complement %in% c(0, 1))  { phenotypes <- phenotypes[, -(1:N), drop = FALSE] }
  }
  # genoSums <- colSums(genotypes, na.rm = TRUE)
  # badGenos <- which(genoSums == 0 | genoSums == nrow(genotypes))
  # if (length(badGenos) > 0) {
  #   print(paste("Eliminating", length(badGenos), "monomorphic genotypes"))
  #   genotypes = genotypes[, -badGenos, drop = FALSE]
  #   if (ncol(genotypes) == 0) {
  #     print("Warning: there are no non-monomorphic genotypes in the input; returning NULLs!")
  #     output <- list(summary = NULL, phenotype = NULL, associations = NULL)
  #     return(output)
  #   }
  # }
  print(paste("There are", ncol(phenotypes), "phenotypes over", nrow(phenotypes), "patients"))
  print(paste("There are", ncol(genotypes), "genotypes over", nrow(genotypes), "patients"))
  print("Computing single associations")
  singleRes <- computeBestPValues(genotypes, phenotypes, cor = FALSE, fisher = TRUE)
  col1 <- colnames(phenotypes)[1]
  inFile <- paste0("BruteForceInput", index, ".tsv")
  testTab <- cbind(genotypes, phenotypes)
  fullTab <- as_tibble(testTab) %>%
    mutate_all(as.numeric) %>%
    bind_cols(ID = IDs) %>% 
    select(ID, everything())
  write_tsv(fullTab, inFile)
  fullResTab <- tibble()
  P <- 0
  for (myK in 1:K) {
    for (myL in 1:L) {
      print(paste("Processing K = ", myK, "and L = ", myL))
      cmdName <- paste0(ifelse(HPC, "/hpc/grid/wip_cmg_systems-immunology/chindl/", "/Users/lchindelevitch/"), 
                        paste0(ifelse(HPC, "", "Downloads/"), "ReverseGWAS/"), "rgwas/rgwas")
      arguments <- c("--input", inFile, "--pheno", col1, "-c", myK, "-v", myL, "--top-k", MAX_TOP, "-t", 1)
      result <- system2(cmdName, args = arguments, stdout = TRUE) %>%
        str_split_fixed(pattern = ",", n = 32) %>%
        magrittr::extract(, 1:10, drop = FALSE)
      coln <- result %>%
        magrittr::extract(1, , drop = TRUE)
      resTab <- result %>%
        magrittr::extract(-1, , drop = FALSE) %>%
        set_colnames(coln) %>%
        as_tibble() %>%
        select(-c(time, gap, opt, status)) %>%
        mutate_at(c("TP", "FP", "FN", "TN"), as.integer) %>%
        computeStats(logSpace = TRUE) %>%
        group_by(SNP) %>%
        arrange(LogFisherExactP, .by_group = TRUE) %>%
        slice(1) %>%
        ungroup %>%
        mutate(K = myK, L = myL) %>%
        select(c(K, L, everything()))
      if (nrow(resTab) > 0) {
        fullResTab %<>%
          bind_rows(resTab)
        P <- P + 1
      }
    }
  }
  extVs <- prepareExtremeValues(log(extremeValue), type, singleRes$value, SINGLE_BEST_RATIO, USE_TIGHTER_BOUND, logSpace = TRUE)
  fullResTab <- fullResTab %>%
    bind_cols(tibble(bestSingle = rep(singleRes$phenotype, P), bestSingleStat = rep(singleRes$value, P), computedBounds = NA)) %>%
    mutate(time = 0, gap = NA, opt = (2 * Accuracy - 1) * Total, status = ifelse(LogFisherExactP < extVs, FEASIBLE, INFEASIBLE)) %>%
    select(SNP, time, gap, opt, status, everything()) %>%
    group_by(SNP) %>%
    arrange(LogFisherExactP, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup %>%
    mutate(K = myK, L = myL)
  if (type == "DNF") {
    fullResTab <- fullResTab %>%
      mutate_at("formula", ~{str_replace_all(., "OR", "AND") %>% str_replace_all("\\)\\ AND\\ \\(", ") OR (")})
  }
  defValue <- ifelse(type == "CNF", TRUE, FALSE)
  complexPs <- matrix(defValue, nrow(phenotypes), ncol(genotypes), dimnames = list(rownames(phenotypes), colnames(genotypes)))
  for (ind in 1:nrow(fullResTab)) {
    curRow <- fullResTab %>% 
      slice(ind)
    if (curRow$status %in% c("feasible", "optimal")) {
      curSNP <- curRow$SNP
      parsedFormula <- parseFormula(curRow$formula, type = type)
      complexPs[, curSNP] <- applyFormula(parsedFormula, phenotypes, type)  
    }
  }
  output <- list(summary = fullResTab, phenotype = complexPs)
  file.remove(inFile)
  if (outputAssociations) { output = c(output, list(associations = singleRes$allValues %>% as_tibble(rownames = "association"))) }
  output
}

testBruteForce = function() {
  testN <- 1200000L
  initMat <- combinat::hcube(rep(2, 4)) - 1
  testMat <- initMat[rep(1:nrow(initMat), testN/nrow(initMat)), ] %>%
    set_colnames(paste0("P", 1:4)) %>%
    as_tibble() %>%
    mutate_all(as.logical) %>%
    mutate(SNP1 = P1 | P2, SNP2 = P1 & P2) %>%
    mutate(ID = as.character(1:testN)) %>%
    select(c(ID, SNP1, SNP2, everything()))
  Q <- bruteForceOptimize(testMat, numSNPs = 2, K = 2, L = 2, outputAssociations = TRUE)
  Q
}
