source("Driver.R")

postprocessSelected = function(inputFile = "UKBB_Diabetes_300kSNPs_feasible_discovery_all_chr_v3_june_full_sweep.csv", 
                               Dir = "/hpc/grid/scratch/chindl/InputFilesDiabetesNonNegated/") {
  resTab = read_csv(inputFile) %>%
    mutate(K = str_extract(file_name_short, "\\_K[0-9]+")) %>%
    mutate(L = str_extract(file_name_short, "\\_L[0-9]+")) %>%
    mutate_at(c("K", "L"), ~{ str_remove_all(., "\\_") %>% str_sub(2) %>% as.integer} ) %>%
    select(K, L, everything())
  print(paste("There are", nrow(resTab), "formulas to process!"))
  initDir = getwd()
  setwd(Dir)
  for (ind in 1:nrow(resTab)) {
    print(ind)
    curLine = resTab %>% slice(ind)
    curType = str_extract(curLine$file_name_short, "(C|D)NF")
    curName = curLine$file_name_short %>%
      str_remove("\\_K[0-9]+") %>%
      str_remove("\\_L[0-9]+") %>%
      str_remove("\\_(C|D)NF") %>%
      paste0("CPLEXInputUKBB_DiabetesComplications_HapMap300_", ., ".csv.gz")
    curFormula = curLine$formula
    curK = curLine$K
    curL = curLine$L
    curSNP = curLine$SNP
    BD = NA
    ext <- str_sub(curName, -7, -4)
    miniFile <- str_sub(curName, 1, -4)
    specString <- paste0(paste(c("", curType, "K", curK, "L", curL), collapse = "_"), ifelse(curSNP != "", paste0("_", curSNP), ""))
    fullFile = str_replace(miniFile, ext, str_c(specString, '_Phenotype.csv'))
    if (!file.exists(fullFile)) {
      valResults = generateValidationPhenotype(curName, curFormula, type = curType, K = curK, L = curL, bestDiscovery = BD, 
                                               outputSummary = FALSE, outputPhenotype = TRUE, outputAssociations = FALSE, SNP = curSNP)
    }
  }
  setwd(initDir)
  return(TRUE)
}

W = postprocessSelected()