source("Settings.R")
source("Driver.R")
source("ILP.R")

COMPLEMENT = FALSE  ### nothing    gets negated
# COMPLEMENT = TRUE ### genotypes  get  negated
# COMPLEMENT = 2    ### everything gets negated
# COMPLEMENT = -1   ### phenotypes get  negated
VALIDATION = TRUE
# VALIDATION = FALSE
SIMULATION = TRUE
# SIMULATION = FALSE

option_list = list(
  make_option(c("-f", "--filename"),   type = "character", default = NULL, help = "full path to input file",    metavar = "string"),
  make_option(c("-t", "--type"),       type = "character", default = NULL, help = "type (CNF/DNF/LCP)",         metavar = "string"),
  make_option(c("-k", "--clauses"),    type = "integer",   default = NULL, help = "number of clauses",          metavar = "integer"),
  make_option(c("-l", "--elements"),   type = "integer",   default = NULL, help = "number of elements/clause",  metavar = "integer"),
  make_option(c("-i", "--index"),      type = "integer",   default = NULL, help = "helps make problem unique",  metavar = "integer"),
  make_option(c("-o", "--objective"),  type = "character", default = NULL, help = "objective (agree/cov/corr)", metavar = "string"))
optParser <- optparse::OptionParser(option_list = option_list)
opt       <- optparse::parse_args(optParser)
print("Here are the options that I parsed from your input:")
print(opt)

if (length(opt) > 1) { ### more than just the help parameter is available
  print(paste("The file I am now processing is", opt$filename))
  ### result    <- mainDriverExtended(inputFile, type = TYPE, objective = OBJECTIVE, complement = COMPLEMENT, KLpairs = KL_PAIRS) ### Uncomment when running everything
  if (VALIDATION) {
    result <- validationDriver(opt$filename, type = opt$type, objective = opt$objective, complement = COMPLEMENT, index = opt$index,
                               Klist = opt$clauses, Llist = opt$elements, shuffle = (!SIMULATION), outputAssociations = (!SIMULATION))
  } else {
    ### result <- mainDriver(opt$filename, type = opt$type, objective = opt$objective, complement = COMPLEMENT, K = opt$clauses, L = opt$elements, index = opt$index)
    result <- preOptimizingDriver(opt$filename, type = opt$type, complement = COMPLEMENT, K = opt$clauses, L = opt$elements, index = opt$index)
  }
} else {
  source("Simulation.R")
  # setwd("../..")
  # Z = prepAllTriage(validation = TRUE)
}
