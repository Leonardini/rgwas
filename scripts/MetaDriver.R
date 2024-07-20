source("Driver.R")
MY_SEED = 123456789L
OBJECTIVE    = "correlation" ### "agreement"
TYPE         = "LCP" ### "CNF"
COMPLEMENT   = FALSE
KL_PAIRS     = rbind(K = c(1,2,3), L = c(1,2,3))

subsampleFile <- function(filename = "CPLEXInput_SelfReported_EU_100000_15.tsv", nPheno = NA, nPat = NA, Snp = NULL, Dir = TRUE) {
  set.seed(MY_SEED)
  numbers <- str_extract_all(filename, '[0-9]+')[[1]]
  outputFile  <- str_replace(filename, numbers[1], as.character(nPat))
  outputFile  <- str_replace(outputFile, numbers[2], paste0(nPheno, "-", ifelse(Dir, "NOT", ""), str_replace(Snp, "_", "-"), "_1"))
  if (!outputFile %in% list.files()) {
    ext <- str_sub(filename, -4)
    readFunction  <- ifelse(ext == '.tsv', read_tsv, read_csv)
    writeFunction <- ifelse(ext == '.tsv', write_tsv, write_csv)
    inputTab <- readFunction(filename, col_types = cols(ID = "c", .default = "l"))
    numSNPs <- str_extract(filename, paste0("([0-9]+)", ext)) %>%
      str_remove(ext) %>%
      parse_integer
    IDs  <- inputTab[, which(colnames(inputTab) == "ID")]
    SNPs <- inputTab[, which(colnames(inputTab) == Snp)]
    if (Dir) {
      SNPs[,1] <- !(SNPs[,1])
      colnames(SNPs) <- paste0("NOT", colnames(SNPs))
    }
    samplePat <- sample(1:nrow(inputTab), size = nPat, replace = FALSE)
    IDs  <- IDs [samplePat, ]
    SNPs <- SNPs[samplePat, ]
    phenotypes  <- inputTab[samplePat, -(1:(1 + numSNPs))]
    uniquePheno <- deleteRepeatedColumns(phenotypes)
    samplePheno <- sample(1:ncol(uniquePheno), size = nPheno, replace = FALSE)
    finalPheno  <- uniquePheno[, samplePheno, drop = FALSE] %>%
      as_tibble
    outputTab   <- bind_cols(IDs, SNPs, finalPheno)
    writeFunction(outputTab, outputFile)
  }
  outputFile
}

option_list = list(
  make_option(c("-p", "--phenotypes"), type = "integer", default = NULL, help = "number of phenotypes",      metavar = "integer"),
  make_option(c("-n", "--patients"),   type = "integer", default = NULL, help = "number of patients",        metavar = "integer"),
  make_option(c("-k", "--clauses"),    type = "integer", default = NULL, help = "number of clauses",         metavar = "integer"),
  make_option(c("-l", "--elements"),   type = "integer", default = NULL, help = "number of elements/clause", metavar = "integer"),
  make_option(c("-i", "--index"),      type = "integer", default = NULL, help = "helps make problem unique", metavar = "integer"),
  make_option(c("-s", "--SNP"),       type = "character",default = NULL, help = "RSID of SNP to process",    metavar = "string"),
  make_option(c("-d", "--direction"),  type = "logical", default = NULL, help = "if TRUE, negate the SNP",   metavar = "logical"))
optParser <- optparse::OptionParser(option_list = option_list)
opt       <- optparse::parse_args(optParser)
print("Here are the options that I parsed from your input:")
print(opt)
if (length(opt) > 1) { ### more than just the help parameter is available
  inputFile <- subsampleFile(nPheno = opt$phenotypes, nPat = opt$patients, Snp = opt$SNP, Dir = opt$direction)
  print(paste("The file I am now processing is", inputFile))
  result    <- mainDriverExtended(inputFile, type = TYPE, objective = OBJECTIVE, complement = COMPLEMENT, KLpairs = KL_PAIRS) ### Uncomment when running everything
  ### result <- mainDriver(inputFile, type = TYPE, objective = OBJECTIVE, complement = COMPLEMENT, K = opt$clauses, L = opt$elements, index = opt$i)
}
