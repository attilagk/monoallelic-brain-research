## Preparations



Load functions


```r
source("~/projects/monoallelic-brain/src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
source("~/projects/monoallelic-brain/src/graphics.R")
```


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
names(e.vars) <- e.vars
```


```r
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
```


```r
set.seed(1976)
perms <- data.frame(cbind(seq_len(nrow(E)), replicate(n.perm <- 1e4, sample.int(nrow(E)))))
names(perms) <- c("U", paste0("P", seq_len(n.perm)))
sel.models <- c("wnlm.Q", "logi.S")
sel.vars <- e.vars[c(1, 3, 5, 8)]
system.time(Betas <- aggregate.CI.permut2(perms = perms, gene.ids = gene.ids, e.vars = e.vars,
                                          sel.vars = sel.vars, sel.models = sel.models,
                                          E = E, Y = Y[gene.ids], skip.CI = TRUE))
```

```
##      user    system   elapsed 
## 18748.892     5.076 18754.811
```

## Results

The expressions above obtain the null distribution of regression coefficients for the selected predictor(s) Age, Gender, Dx, Ancestry.1 under the selected model(s) wnlm.Q, logi.S based on 10<sup>4</sup> permutations.  Each panel within a plot shows a red vertical zero line of the null hypothesis $\beta_j = 0$, the null distribution as probability density (solid blue line) based on the permutations, and the corresponding p-value (dotted blue line and number on the bottom).

### Under wnlm.Q


```r
beta0densityplot("Age", mtype = "wnlm.Q")
```

<img src="figure/beta-age-null-wnlm-Q-1.png" title="plot of chunk beta-age-null-wnlm-Q" alt="plot of chunk beta-age-null-wnlm-Q" width="700px" />

<img src="figure/beta-gender-null-wnlm-Q-1.png" title="plot of chunk beta-gender-null-wnlm-Q" alt="plot of chunk beta-gender-null-wnlm-Q" width="700px" />

<img src="figure/beta-dx-scz-wnlm-Q-1.png" title="plot of chunk beta-dx-scz-wnlm-Q" alt="plot of chunk beta-dx-scz-wnlm-Q" width="700px" />

<img src="figure/beta-dx-aff-wnlm-Q-1.png" title="plot of chunk beta-dx-aff-wnlm-Q" alt="plot of chunk beta-dx-aff-wnlm-Q" width="700px" />

<img src="figure/beta-ancestry-1-null-wnlm-Q-1.png" title="plot of chunk beta-ancestry-1-null-wnlm-Q" alt="plot of chunk beta-ancestry-1-null-wnlm-Q" width="700px" />

### Under logi.S

<img src="figure/beta-age-null-logi-S-1.png" title="plot of chunk beta-age-null-logi-S" alt="plot of chunk beta-age-null-logi-S" width="700px" />

<img src="figure/beta-gender-null-logi-S-1.png" title="plot of chunk beta-gender-null-logi-S" alt="plot of chunk beta-gender-null-logi-S" width="700px" />

<img src="figure/beta-dx-scz-logi-S-1.png" title="plot of chunk beta-dx-scz-logi-S" alt="plot of chunk beta-dx-scz-logi-S" width="700px" />

<img src="figure/beta-dx-aff-logi-S-1.png" title="plot of chunk beta-dx-aff-logi-S" alt="plot of chunk beta-dx-aff-logi-S" width="700px" />

<img src="figure/beta-ancestry-1-null-logi-S-1.png" title="plot of chunk beta-ancestry-1-null-logi-S" alt="plot of chunk beta-ancestry-1-null-logi-S" width="700px" />

## Comparison to p-values from t-distribution

The calculations and the plot below compare the p-values from the permutations to those based on the theoretical t-distribution of the statistic.  Comparing these results to those from model checking reveals that the two approaches to p-value calculation agree (fall near the gray diagonal on the plots) as long as the model fits the data well.  Otherwise the t-distribution based approach tends to produce much lower p-values and therefore exaggerate the significance that $\beta_j\neq 0$.


```r
source("2016-10-03-permutation-test.R")
```


```r
M <- do.all.fits(Z = Y[gene.ids], G = E, preds = e.vars, sel.models = sel.models)
cf <- unlist(sapply(e.vars, function(e.v) predictor2coefs(M[[c(1, 1)]], e.v)))
both.p.val <-
    do.call(rbind,
            lapply(sel.models, function(mtype)
                   do.call(rbind,
                           lapply(cf, function(coef)
                                  do.call(rbind,
                                          lapply(gene.ids,
                                                 function(g) get.both.p.vals(mtype = mtype, gene = g, coef = coef, M = M, B = Betas)))))))
```

<img src="figure/p-val-tdist-vs-perm-1.png" title="plot of chunk p-val-tdist-vs-perm" alt="plot of chunk p-val-tdist-vs-perm" width="700px" />

Results are written to csv file


```r
write.csv(both.p.val, "../../results/regr-coefs.csv", row.names = FALSE)
```