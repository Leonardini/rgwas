### Converts warnings into errors
options(warn = 0)

### List of rerquired packages
library(conflicted) ### used to resolve conflicts between the imported packages
library(admisc)     ### used to simplify logical combinations of variables
library(BiRewire)   ### used to simulate randomized versions of input matrices
library(combinat)   ### used to create all possible combinations of variables
library(glpkAPI)    ### used to define variables and constraints within an ILP
library(gtools)     ### used to conduct efficient binary search in simulations
library(magrittr)   ### used to simplify piping and allow short function names
library(Matrix)     ### used to compute cross-products and QR decompositions
library(optparse)   ### used to parse command-line options if run in batch mode
library(Rcplex2)    ### used to create optimization problems directly via an API
library(slam)       ### used to create sparse matrices in a convenient manner
library(tidyverse)  ### used to carry out tidy analysis of the input and output

conflicts_prefer(dplyr::filter)
conflicts_prefer(magrittr::set_names)

### Settings for driver functions in Driver.R
HPC = TRUE
OPTIMIZATION = TRUE ### TRUE means optimization mode; FALSE means triage mode
MY_SEED = 123456789L  
MY_TYPES = "DNF" 
MY_OBJECTIVES = "agreement" 
if (HPC) {
  KLIST = 1:3
  LLIST = 1:3
  MY_K = KLIST
  MY_L = LLIST
} else {
  MY_K = 3
  MY_L = 3
}
EXT = ".csv"
LOW_P = TRUE
MAX_P = ifelse(LOW_P, 1e-3, 5e-8)
SINGLE_BEST_RATIO = Inf   ### target ratio of the combined phenotype's summary statistic to that of the best singleton phenotype
USE_TIGHTER_BOUND = FALSE ### set to TRUE if the tighter (stronger) bound should be used, FALSE if the weaker one should be used
THOUSAND = 1000
TOP_DIR = Sys.getenv(x = "CPLEX_STUDIO_DIR201", unset = "", names = NA)
CPLEX_DIR = paste0(ifelse(TOP_DIR != "", TOP_DIR, paste0(ifelse(HPC, "/hpc/grid/wip_cmg_systems-immunology/", "/Applications/"), 
                                                         "CPLEX_Studio_Community201")), "/cplex/bin/x86-64_", ifelse(HPC, "linux", "osx"), "/")

### Settings for geometric functions in Geometry.R
WIDTH = 0.999 ### width of the field determining how much we can be off by from the boundary
M_STEP = 2L
ABOVE = 0L
BELOW = 1L
POSITIVE = 1L
NEGATIVE = -1L
SIGNS = c(positive = POSITIVE, negative = NEGATIVE)
TOL = 1e-5

### ILP and CPLEX solver settings in ILP.R
DELTA = 0.001 ### minimal amount by which to exceed the trivial bound
### internal names for ILP statuses
OPTIMAL    = "optimal"
INFEASIBLE = "infeasible"
TIMEOUT    = "timeout"
FEASIBLE   = "feasible"
MEMORYOUT  = "nomemory"
HOPELESS   = "hopeless"
### prefixes for variable names
PREFIX1 = "U"
PREFIX2 = "I"
### CPLEX solver settings, high-level
OPPORTUNISTIC = -1
AUTOMATIC     = 0
DETERMINISTIC = 1 # to be used when reproducibility is needed
PAR_MODE      = OPPORTUNISTIC
NUM_THREADS   = ifelse(HPC, 3L, 6L)
### CPLEX solver settings, low-level
TIME_LIMIT   = ifelse(OPTIMIZATION, 3600,  3600)
TREE_MEMORY  = ifelse(OPTIMIZATION, 20000, 7000)
PROBE_LEVEL  = ifelse(OPTIMIZATION, -1L, NA)
EMPHASIS     = ifelse(OPTIMIZATION,  3L, NA)
CUT_LEVEL    = ifelse(OPTIMIZATION,  2L, NA)
WRITE_LEVEL  = ifelse(OPTIMIZATION,  1L, NA)
REPAIR_TRIES = ifelse(OPTIMIZATION, -1L, NA)
VAR_SELECT   = ifelse(OPTIMIZATION,  4L, NA)
CUT_PASSES   = ifelse(OPTIMIZATION,  1L, NA)
HEURISTIC_FR = ifelse(OPTIMIZATION, -1L, NA)
NUM_SOL      = ifelse(OPTIMIZATION,  NA, 1L)
NUM_EMPHASIS = ifelse(OPTIMIZATION,  NA, 1L)
INTTOLERANCE = ifelse(OPTIMIZATION,  1e-7, 1e-7)

### Simulation settings in Simulation.R
MY_SEED  = 123456789L
PHENO_NAMES = paste0(c("RA", "CD", "UC", "IBD", "PsO", "AS", "Thyroid", "Celiac", "Lupus", "MS", "T1D"), ".pheno")