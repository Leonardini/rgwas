N_PHENOTYPES = c(30, 100, 300)
N_PATIENTS   = c(10000, 30000, 90000)
GOOD_SNPS    = rbind(SNP = c("rs140135_C", "rs1569414_T"), DIR = c(FALSE, TRUE))

makeBashScript = function(filename = "FullScript.sh") {
  n1 <- length(N_PHENOTYPES)
  n2 <- length(N_PATIENTS)
  n3 <- ncol(GOOD_SNPS)
  N <- n1 * n2 * n3
  Lines <- rep("", N + 1)
  Lines[1]  = "#!/bin/bash"
  pos <- 2
  for (ind1 in 1:n1) {
    curP <- N_PHENOTYPES[ind1]
    for (ind2 in 1:n2) {
      curN <- N_PATIENTS[ind2]
      for (ind3 in 1:n3) {
        curSNP <- GOOD_SNPS[, ind3]
        curS <- curSNP[1]
        curD <- curSNP[2]
        Lines[pos] <- paste("Rscript", "MetaDriver.R", "-p", curP, "-n", curN, "-s", curS, "-d", curD, rep(" &", (pos %% 2 == 0))) 
        pos <- pos + 1
      }
    }
  }
  fn <- file(filename)
  writeLines(Lines, fn)
  close(fn)
}

makeBashScriptForTimeouts = function(filename = "FullScriptForTimeouts.sh", sourceFile = "~ziemed/UKBB/CPLEXInput/SummarizedSweepResults.tsv", Dir = "~ziemed/UKBB/CPLEXInput/Non-RDSForDaniel/") {
  Tab <- read_tsv(sourceFile) %>% 
    filter(status == "timeout" & is.na(formula))
  N <- nrow(Tab)
  Lines <- rep("", N + 1)
  Lines[1]  = "#!/bin/bash"
  pos <- 2
  for (ind in 1:N) {
    curRow <- Tab %>%
      slice(ind)
    curP <- curRow$numPhenos
    curN <- curRow$numPatients
    curS <- curRow$SNP %>%
      str_remove("NOT")
    curD <- (str_sub(curRow$SNP, 1, 3) == "NOT")
    curK <- curRow$K
    curL <- curRow$L
    Lines[pos] <- paste("Rscript", "MetaDriver.R", "-p", curP, "-n", curN, "-s", curS, "-d", curD, "-k", curK, "-l", curL, "-i", pos - 1, rep(" &", (pos %% 3 != 1))) 
    pos <- pos + 1
  }
  fn <- file(filename)
  writeLines(Lines, fn)
  close(fn)
}