inputFile <- "CPLEXInputs/CPLEXInput_UKBB_ICD10_SmallerSet_HapMap300_chr13_subset6b_341485_192.csv.gz"
MAX_COMBOS_BF <- 0
type = "CNF"
useSNPs = 5
VALIDATION = TRUE
extremeValue = ifelse(type == "LCP", MIN_CORREL, MAX_P)
miniFile <- inputFile
if (str_sub(inputFile, -3) == ".gz") {
  ext <- str_sub(inputFile, -7, -4)
  miniFile <- str_sub(inputFile, 1, -4)
}
readFunction <- ifelse(ext == '.tsv', read_tsv, read_csv)
inputTab <- readFunction(inputFile, col_types = cols(ID = "c", id = "c", .default = ifelse(type == "LCP", "d", "l")))
if ("id" %in% colnames(inputTab)) { inputTab <- inputTab %>% rename(ID = id) }
numSNPs <- str_extract(miniFile, paste0("([0-9]+)", ext)) %>% str_remove(ext) %>% parse_integer
baseFilename <- str_remove(miniFile, ext)
numPheno <- ncol(inputTab) - 1 - numSNPs
phenoCols <- inputTab[, -(1:(1 + numSNPs))]
cSums <- colSums(phenoCols, na.rm = TRUE)
myOrder <- order(cSums)
phenoCols <- phenoCols[, myOrder]
cSums = cSums[myOrder]
startPos <- floor(numPheno/2)
endPos   <- numPheno
fullOutputs <- vector('list')
for (ind in startPos:endPos) {
  print(ind)
  curInputTab <- bind_cols(inputTab[, 1:(1 + useSNPs)], phenoCols[, 1:ind])
  if (VALIDATION) {
    set.seed(MY_SEED)
    firstHalf = which(runif(nrow(curInputTab)) <= 0.5)
    curInputTab = curInputTab %>% slice(firstHalf)
  }
  curOutput <- createAndSolveILPs(curInputTab, numSNPs = useSNPs, type = type, objective = "agreement", K = 2, L = 2, index = 0, 
                             logSpace = TRUE, extremeValue = ifelse(type == "LCP", extremeValue, log(extremeValue)), 
                             complement = 0, baseFilename = NULL, outputAssociations = FALSE)
  print(paste("The average time required was", mean(curOutput$summary$time)))
  fullOutputs[[ind + 1 - startPos]] <- curOutput
}
