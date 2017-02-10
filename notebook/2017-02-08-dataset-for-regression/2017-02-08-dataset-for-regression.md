This short workflow exports the data onto which various GLMs have been fitted: the read count ratio for 30 genes (as response variables) and 12 predictors across 579 samples.  The dataset (as a CSV file) was sent to Gabriel Hoffman for potential analysis.


```r
source("../../src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
```


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
names(e.vars) <- e.vars
E <- get.predictors()
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
S <- cbind(data.frame(sapply(Y[gene.ids], getElement, "S")), E[e.vars])
Q <- cbind(data.frame(sapply(Y[gene.ids], getElement, "Q")), E[e.vars])
N <- cbind(data.frame(sapply(Y[gene.ids], getElement, "N")), E[e.vars])
write.csv(S, file = "S-30-genes-12-predictors.csv")
write.csv(Q, file = "Q-30-genes-12-predictors.csv")
write.csv(N, file = "N-30-genes-12-predictors.csv")
```
