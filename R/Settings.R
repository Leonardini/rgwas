#' Settings for driver functions in Driver.R and Validation.R
OPTIMIZATION = TRUE #' TRUE means optimization mode; FALSE means triage mode
MY_SEED = 123456789L
MY_TYPES = "DNF"
MY_OBJECTIVES = "agreement"
KLIST = c(rep(1:3, 3),        1, 5)
LLIST = c(rep(1:3, each = 3), 5, 1)
MY_K = KLIST
MY_L = LLIST
EXT = ".csv"
LOW_P = TRUE
LOW_PVALUE  = 1e-3
HIGH_PVALUE = 5e-8
MAX_P = ifelse(LOW_P, LOW_PVALUE, HIGH_PVALUE)
THOUSAND = 1000

#' Settings for geometric functions in Geometry.R
#' Width of the field determining how much we can be off by from the boundary
WIDTH = 0.999
M_STEP = 2L
ABOVE = 0L
BELOW = 1L
POSITIVE = 1L
NEGATIVE = -1L
SIGNS = c(positive = POSITIVE, negative = NEGATIVE)
TOL = 1e-5

#' ILP and CPLEX solver settings in ILP.R
#' minimal amount by which to exceed the trivial bound
DELTA = 0.001
#' internal names for ILP statuses
OPTIMAL    = "optimal"
INFEASIBLE = "infeasible"
TIMEOUT    = "timeout"
FEASIBLE   = "feasible"
MEMORYOUT  = "nomemory"
HOPELESS   = "hopeless"
#' prefixes for variable names
PREFIX1 = "U"
PREFIX2 = "I"
#' CPLEX solver settings, high-level
OPPORTUNISTIC = -1
AUTOMATIC     = 0
DETERMINISTIC = 1 # to be used when reproducibility is needed
NUM_THREADS   = 6L
#' CPLEX solver settings, low-level
TIME_LIMIT   = 3600L
EMPHASIS     = 3L
PROBE_LEVEL  = -1L
VAR_SELECT   = 4L
CUT_LEVEL    = 2L
HEURISTIC_FR = -1L
