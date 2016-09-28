


```r
source("~/projects/monoallelic-brain/src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
source("~/projects/monoallelic-brain/src/graphics.R")
source("../2016-08-21-likelihood-surface/2016-08-21-likelihood-surface.R")
```

## Near collinearity: RIN and RIN2


```r
# explanatory variables (a.k.a. predictors)
e.vars2 <- c("Age",
               "Institution",
               "Gender",
               "PMI",
               "Dx",
               "RIN",
               "RIN2",
               "RNA_batch",
               "Ancestry.1", "Ancestry.2", "Ancestry.3", "Ancestry.4", "Ancestry.5" )
e.vars2[7]
```

```
## [1] "RIN2"
```

```r
e.vars <- e.vars2[-7]
```


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
```


```r
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
```


```r
# exclude unweighed aggregates UA.8 and UA from fitting
to.fit.ids <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
sel.models <- c("logi.S", "wnlm.Q", "wnlm.R"); names(sel.models) <- sel.models
M <- do.all.fits(Y[to.fit.ids], G = E, preds = e.vars, sel.models = sel.models)
M2 <- do.all.fits(Y[to.fit.ids], G = E, preds = e.vars2, sel.models = sel.models)
```

### Terms

#### PEG3: RIN and RIN2


```r
par(mfrow = c(3, 4))
termplot(M2$wnlm.Q$PEG3, terms = c(1:5, 8:12, 6:7))
```

<img src="figure/terms-PEG3-RIN-RIN2-1.png" title="plot of chunk terms-PEG3-RIN-RIN2" alt="plot of chunk terms-PEG3-RIN-RIN2" width="700px" />

#### PEG3: only RIN

<img src="figure/terms-PEG3-RIN-1.png" title="plot of chunk terms-PEG3-RIN" alt="plot of chunk terms-PEG3-RIN" width="700px" />

#### KCNK9: RIN and RIN2

<img src="figure/terms-KCNK9-RIN-RIN2-1.png" title="plot of chunk terms-KCNK9-RIN-RIN2" alt="plot of chunk terms-KCNK9-RIN-RIN2" width="700px" />

#### KCNK9: only RIN

<img src="figure/terms-KCNK9-RIN-1.png" title="plot of chunk terms-KCNK9-RIN" alt="plot of chunk terms-KCNK9-RIN" width="700px" />

### Estimates and CIs


```r
Betas <- lapply(M, function(m) { x <- get.estimate.CI(m); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
Betas2 <- lapply(M2, function(m) { x <- get.estimate.CI(m); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
```

<img src="figure/betas-RIN-RIN2-1.png" title="plot of chunk betas-RIN-RIN2" alt="plot of chunk betas-RIN-RIN2" width="700px" />

<img src="figure/betas-RIN-1.png" title="plot of chunk betas-RIN" alt="plot of chunk betas-RIN" width="700px" />

### Likelihood surface


```
## Error in m$model[["Y"]][, 1]: incorrect number of dimensions
```



<img src="figure/ll-surf-PEG3-RIN-RIN2-1.png" title="plot of chunk ll-surf-PEG3-RIN-RIN2" alt="plot of chunk ll-surf-PEG3-RIN-RIN2" width="700px" />

<img src="figure/ll-surf-PEG3-RIN-1.png" title="plot of chunk ll-surf-PEG3-RIN" alt="plot of chunk ll-surf-PEG3-RIN" width="700px" />
