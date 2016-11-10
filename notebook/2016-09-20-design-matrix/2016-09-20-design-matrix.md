## Motivation

Previous analysis used a design matrix $X$ with the following shortcomings

1. the baseline level of certain factors in $X$ was automatically set; in particular for Dx the baseline level was AFF instead of Control
1. $X$ contained a pair of nearly collinear predictors RIN and RIN2 (squared)

These shortcomings are corrected here and the corresponding changes in results are studied.  The results show that the corrections afford minor improvements:

1. SCZ is now directly contrasted to Control, revealing significant effect for some genes
1. with the removal of RIN2 the term associated with RIN is no more much larger than those with other predictors

## Results

### Preliminaries


```
## Stackoverflow is a great place to get help:
## http://stackoverflow.com/tags/ggplot2.
```

```
## Loading required package: RColorBrewer
```

```
## 
## Attaching package: 'latticeExtra'
```

```
## The following object is masked from 'package:ggplot2':
## 
##     layer
```


```r
source("../../src/import-data.R")
source("../../src/fit-glms.R")
source("../../src/likelihood-surface.R")
source("../../src/graphics.R")
```


```r
# explanatory variables (a.k.a. predictors)
e.vars2 <- c("Age",
               "Institution",
               "Gender",
               "PMI",
               "Dx",
               "RIN",
               "RIN2", # to be removed
               "RNA_batch",
               "Ancestry.1", "Ancestry.2", "Ancestry.3", "Ancestry.4", "Ancestry.5" )
e.vars <- e.vars2[-7] # removing RIN2
```




```r
E <- get.predictors() # default arguments
E.al <- get.predictors(adj.levels = FALSE) # automatic (unadjusted) levels
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
to.fit.ids <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE) # exclude unweighed aggregates UA.8 and UA from fitting
sel.models <- c("logi.S", "wnlm.Q", "wnlm.R", "unlm.Q"); names(sel.models) <- sel.models
M <- do.all.fits(Y[to.fit.ids], G = E, preds = e.vars, sel.models = sel.models)
M.al <- do.all.fits(Y[to.fit.ids], G = E.al, preds = e.vars, sel.models = sel.models)
M2 <- do.all.fits(Y[to.fit.ids], G = E, preds = e.vars2, sel.models = sel.models)
```

### Estimates and CIs


```r
Betas <- lapply(M, function(m) { x <- get.estimate.CI(m); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
Betas.al <- lapply(M.al, function(m) { x <- get.estimate.CI(m); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
Betas2 <- lapply(M2, function(m) { x <- get.estimate.CI(m); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
```

<img src="figure/betas-RIN-1.png" title="plot of chunk betas-RIN" alt="plot of chunk betas-RIN" width="700px" />

<img src="figure/betas-RIN-RIN-autolev-1.png" title="plot of chunk betas-RIN-RIN-autolev" alt="plot of chunk betas-RIN-RIN-autolev" width="700px" />

<img src="figure/betas-RIN-RIN2-1.png" title="plot of chunk betas-RIN-RIN2" alt="plot of chunk betas-RIN-RIN2" width="700px" />

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

### Likelihood surface





<img src="figure/ll-surf-wnlm-Q-PEG3-RIN-RIN2-1.png" title="plot of chunk ll-surf-wnlm-Q-PEG3-RIN-RIN2" alt="plot of chunk ll-surf-wnlm-Q-PEG3-RIN-RIN2" width="700px" />

<img src="figure/ll-surf-wnlm-Q-PEG3-RIN-1.png" title="plot of chunk ll-surf-wnlm-Q-PEG3-RIN" alt="plot of chunk ll-surf-wnlm-Q-PEG3-RIN" width="700px" />
