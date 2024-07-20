library(tidyverse)

phenoPath ="/Users/ziemed/UKBB/phenotypes/SelfReporting_Codes/"
tab = read_rds(paste0(phenoPath, "Field_20002_Self_Report_Matrix.RDS"))
               