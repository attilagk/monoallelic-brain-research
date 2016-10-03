### Preparations



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
perms <- data.frame(cbind(seq_len(nrow(E)), replicate(n.perm <- 1e3, sample.int(nrow(E)))))
names(perms) <- c("U", paste0("P", seq_len(n.perm)))
sel.models <- c("wnlm.Q")
sel.vars <- e.vars[c(1, 3, 5, 8)]
system.time(Betas <- aggregate.CI.permut2(perms = perms, gene.ids = gene.ids, e.vars = e.vars,
                                          sel.vars = sel.vars, sel.models = sel.models,
                                          E = E, Y = Y[gene.ids], skip.CI = TRUE))
```

```
##    user  system elapsed 
## 587.408   0.040 587.412
```

The expressions above obtain the null distribution of regression coefficients for the selected predictor(s) Age, Gender, Dx, Ancestry.1 under the selected model(s) wnlm.Q based on 1000 permutations.


```r
beta0densityplot("Age")
```

<img src="figure/beta-age-null-wnlm-Q-1.png" title="plot of chunk beta-age-null-wnlm-Q" alt="plot of chunk beta-age-null-wnlm-Q" width="700px" />


```r
beta0densityplot("GenderMale")
```

<img src="figure/beta-gender-null-wnlm-Q-1.png" title="plot of chunk beta-gender-null-wnlm-Q" alt="plot of chunk beta-gender-null-wnlm-Q" width="700px" />


```r
beta0densityplot("DxSCZ")
```

<img src="figure/beta-dx-scz-wnlm-Q-1.png" title="plot of chunk beta-dx-scz-wnlm-Q" alt="plot of chunk beta-dx-scz-wnlm-Q" width="700px" />


```r
beta0densityplot("DxAFF")
```

<img src="figure/beta-dx-aff-wnlm-Q-1.png" title="plot of chunk beta-dx-aff-wnlm-Q" alt="plot of chunk beta-dx-aff-wnlm-Q" width="700px" />


```r
beta0densityplot("Ancestry.1")
```

<img src="figure/beta-ancestry-1-null-wnlm-Q-1.png" title="plot of chunk beta-ancestry-1-null-wnlm-Q" alt="plot of chunk beta-ancestry-1-null-wnlm-Q" width="700px" />
