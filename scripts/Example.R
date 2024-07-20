library(tidyverse)
library(magrittr)
library(combinat)
RERUN = FALSE

if (RERUN) {
  bestP = Inf
  bestQ = c(A = Inf, B = Inf, C = Inf, D = Inf)
  for (A in 1:50) {
    for (B in 1:50) {
      print(c(A, B))
      for (C in 2:50) {
        if (A + B + C < 50) {
          D = 50 - (A + B + C)
          M1 = rbind(c(A, B), c(C, D))
          for (E in (1 - A):(B - 1)) {
            M2 = rbind(c(C - 1, D + 1), c(B - E, A + E))
            M12 = M1 + M2
            if (fisher.test(M12, alternative = "g")$p.value < 5e-8) {
              Z = array(NA, c(2,2,2))
              Z[, , 1] = M1
              Z[, , 2] = M2
              curP = mantelhaen.test(Z, alternative = "g", exact = TRUE)$p.value
              if (curP < bestP) {
                bestP = curP
                bestQ = c(A, B, C, D, E)
                # print("Found an improvement!")
              }
            }
          }
        }
      }
    }
  }
}

bestQ = c(27, 1, 21, 1, 0)
names(bestQ) = LETTERS[1:5]
bestQ %<>% as.list()
A <- bestQ$A
B <- bestQ$B
C <- bestQ$C
D <- bestQ$D
E <- bestQ$E
M1 = rbind(c(A, B), c(C, D))
M2 = rbind(c(C - 1, D + 1), c(B - E, A + E))
M12 = M1 + M2
M0 = rbind(colSums(M1), colSums(M2))
Z = array(NA, c(2,2,2))
Z[, , 1] = M1
Z[, , 2] = M2
PvalBestSingle = fisher.test(M0, alternative = "g")$p.value
print(paste("The best single phenotype-genotype association has p-value", signif(PvalBestSingle, 3)))
PvalCombined   = fisher.test(M12, alternative = "g")$p.value
print(paste("The composite phenotype-genotype association has p-value", signif(PvalCombined, 3)))
PvalAdjusted   = mantelhaen.test(Z, alternative = "g", exact = TRUE)$p.value
print(paste("The composite phenotype-genotype association stratified by the single phentoype has p-value", signif(PvalAdjusted, 3)))

Q = (combinat::hcube(rep(2, 3)) - 1) %>%
  set_colnames(c("Genotype", "SingleBestPhenotype", "CompositePhenotype")) %>%
  as_tibble() %>%
  mutate_all(as.logical)

Multiplicities = c(A + E, B - E, D, C, D + 1, C - 1, B, A)

extendedTab = Q[rep(1:8, Multiplicities), ]

### Sanity checks to make sure the numbers work out
Tab1 = table(extendedTab$Genotype, extendedTab$SingleBestPhenotype)
pvalBestSingle = fisher.test(Tab1, alternative = "g")$p.value
stopifnot(near(pvalBestSingle, PvalBestSingle))

Tab2 = table(extendedTab$Genotype, extendedTab$CompositePhenotype)
pvalCombined   = fisher.test(Tab2, alternative = "g")$p.value
stopifnot(near(pvalCombined, PvalCombined))

Tab3 = table(extendedTab$Genotype, extendedTab$CompositePhenotype, extendedTab$SingleBestPhenotype)
pvalAdjusted   = mantelhaen.test(Tab3, alternative = "g", exact = TRUE)$p.value
stopifnot(near(pvalAdjusted, PvalAdjusted))