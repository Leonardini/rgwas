library(tidyverse)

print("WRITING TOE SNPS...STARTS")

#toe = read_tsv("/Users/ziemed/UKBB/toe/toe_best_guess.txt")  %>%
toe = read_tsv("toe/toe_best_guess.txt") %>%
  filter(pvalue < 5e-8, chr==22) %>% 
  filter(str_detect(trait, "nflamm"))
allSnps = unique(c(toe$oldsnp, toe$newsnp))
write_lines(allSnps, "genotypes/derived/toeSNPs_Chr22.txt")

print("WRITING TOE SNPS...DONE")

