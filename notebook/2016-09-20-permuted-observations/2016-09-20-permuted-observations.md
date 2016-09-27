## Introduction

There are several limitations of the (generalized) linear models (e.g. wnlm.R, logi.S) of read count ratio $S$ regressed on predictors such as Age, Institution, etc.  These limitations might lead to bias in the estimation of regression coefficients $\beta_j$, on which some of the main conclusions from this study are based.  In particular, a nonzero $\beta_j$ associated with some predictor indicates either significant effect, or bias.  These are only the basic cases; typically we expect to deal with some combination of the two.

The analysis herein aims to distinguish between the two basic cases.  This done by randomly re-permuting observations for each predictor (in other words individuals are randomly relabeled for each predictor).  Any significant result associated with the permuted predictor must either come from bias or occur "by chance" with probability 1 less the chosen confidence level for $\beta_j$.

## Calculations

### Preparations



Load functions


```r
source("~/projects/monoallelic-brain/src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
source("~/projects/monoallelic-brain/src/graphics.R")
source("2016-09-20-permuted-observations.R")
```


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
```

Get data: observations on predictors (explanatory variables) and on the higher and lower read count from selected genes (more details in a previous post)

```r
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
```

### Permutations


```r
set.seed(1976)
perm.obs <- sample.int(nrow(E))
names(e.vars) <- e.vars
EP <- lapply(e.vars, function(v) { E1 <- E; E1[[v]] <- E[[v]][perm.obs]; return(E1) })
EP$Unpermuted <- E
```


```r
set.seed(1976)
perms <- data.frame(cbind(seq_len(nrow(E)), replicate(20, sample.int(nrow(E)))))
names(perms) <- c("U", paste0("P", seq_len(length(perms) - 1)))
```

### Model fitting and CI calculation


```r
# exclude unweighed aggregates UA.8 and UA from fitting
to.fit.ids <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
sel.models <- c("logi.S", "wnlm.R"); names(sel.models) <- sel.models
M <- lapply(EP, function(e) do.all.fits(Y[to.fit.ids], G = e, preds = e.vars, sel.models = sel.models))
```


```r
Betas.Unpermuted <- lapply(M$Unpermuted,
                           function(m) { x <- get.estimate.CI(m); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
Betas.Unpermuted.95 <- lapply(M$Unpermuted,
                           function(m) { x <- get.estimate.CI(m, conf.lev = 0.95); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
```


```r
Betas.Permuted <- lapply(sel.models, get.estimate.CI.permut, M = M, e.v = e.vars, conf.lev = 0.99)
Betas.Permuted.95 <- lapply(sel.models, get.estimate.CI.permut, M = M, e.v = e.vars, conf.lev = 0.95)
```


```r
Betas <- aggregate.CI.permut2(perms = perms, gene.ids = gene.ids, e.vars = e.vars,
                              sel.models = sel.models, E = E, Y = Y[gene.ids], conf.lev=0.99)
```

## Results

Both the *unpermuted* and *permuted* data was fit by two models: logi.S and wnlm.R.  These four combinations are shown in the four plots below at confidence level of $99\%$  The following observations and conclusions may be noted:

Under wnlm.R there are 43 significant coefficients before and 27 after permutation.  In contrast, under logi.S there are 165 and 180 (before and after).  But taking only those coefficients that are significant under *both* models yields 42 and 24 (before and after, respectively).

Under logi.S the following pattern is observed.  For some genes (e.g. GRB10, ZDBF2) permutation tends to abolish significance---if observed---, whereas for other genes (e.g. SNRPN, SNURF) there are many significant coefficients after permutation.  This suggests systematic differences between genes: better fit of logi.S to the data for the former set of genes and poorer for the latter gene set.  The poor fit likely results in biased coefficient estimates, which explains the many significant coefficients after permutation.


```r
my.segplot(data = Betas.Unpermuted$logi.S, xlim = my.xlim, main = "Unpermuted under logi.S")
```

<img src="figure/unpermuted-logi-S-1.png" title="plot of chunk unpermuted-logi-S" alt="plot of chunk unpermuted-logi-S" width="700px" />

<img src="figure/permuted-logi-S-1.png" title="plot of chunk permuted-logi-S" alt="plot of chunk permuted-logi-S" width="700px" />

<img src="figure/unpermuted-wnlm-R-1.png" title="plot of chunk unpermuted-wnlm-R" alt="plot of chunk unpermuted-wnlm-R" width="700px" />

<img src="figure/permuted-wnlm-R-1.png" title="plot of chunk permuted-wnlm-R" alt="plot of chunk permuted-wnlm-R" width="700px" />

### Repeated permutation

#### Age


```r
my.segplot2("Age", "logi.S", -2) # fit has not converged for gene ranked 2 (TMEM261P1)
```

<img src="figure/permuted-age-logi-S-1.png" title="plot of chunk permuted-age-logi-S" alt="plot of chunk permuted-age-logi-S" width="700px" />

<img src="figure/permuted-age-wnlm-R-1.png" title="plot of chunk permuted-age-wnlm-R" alt="plot of chunk permuted-age-wnlm-R" width="700px" />

#### Gender

<img src="figure/permuted-gender-logi-S-1.png" title="plot of chunk permuted-gender-logi-S" alt="plot of chunk permuted-gender-logi-S" width="700px" />

<img src="figure/permuted-gender-wnlm-R-1.png" title="plot of chunk permuted-gender-wnlm-R" alt="plot of chunk permuted-gender-wnlm-R" width="700px" />

#### Dx

<img src="figure/permuted-dx-control-logi-S-1.png" title="plot of chunk permuted-dx-control-logi-S" alt="plot of chunk permuted-dx-control-logi-S" width="700px" />

<img src="figure/permuted-dx-control-wnlm-R-1.png" title="plot of chunk permuted-dx-control-wnlm-R" alt="plot of chunk permuted-dx-control-wnlm-R" width="700px" />

#### Ancestry

<img src="figure/permuted-ancestry-1-logi-S-1.png" title="plot of chunk permuted-ancestry-1-logi-S" alt="plot of chunk permuted-ancestry-1-logi-S" width="700px" />

<img src="figure/permuted-ancestry-1-wnlm-R-1.png" title="plot of chunk permuted-ancestry-1-wnlm-R" alt="plot of chunk permuted-ancestry-1-wnlm-R" width="700px" />

#### Institution

<img src="figure/permuted-institution-penn-logi-S-1.png" title="plot of chunk permuted-institution-penn-logi-S" alt="plot of chunk permuted-institution-penn-logi-S" width="700px" />

<img src="figure/permuted-institution-penn-wnlm-R-1.png" title="plot of chunk permuted-institution-penn-wnlm-R" alt="plot of chunk permuted-institution-penn-wnlm-R" width="700px" />

#### PMI

<img src="figure/permuted-pmi-logi-S-1.png" title="plot of chunk permuted-pmi-logi-S" alt="plot of chunk permuted-pmi-logi-S" width="700px" />

<img src="figure/permuted-pmi-wnlm-R-1.png" title="plot of chunk permuted-pmi-wnlm-R" alt="plot of chunk permuted-pmi-wnlm-R" width="700px" />



## Conclusions

1. As discussed before, wnlm.R appears less sensitive but less biased than logi.S.
1. The bias appears greatly vary across genes; model checking will help filter out genes with much bias due to poor fit
1. It is beneficial to use a consensus approach by defining significance that is fulfilled both under wnlm.R and logi.S.
