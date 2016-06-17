All functions have been reimplemented partly because the new .csv files created from html tables are organized differently than the ones from Ifat, and partly to make code simpler, more usable and maintainable.

### New importer functions

Using the new implementation...

```r
source("~/projects/monoallelic-brain/src/import-data.R")
```


```r
# default arguments given explictely to both function calls
E <- get.predictors()
Y <- get.readcounts(gene.ids = gene.ids)
```

```
## Warning in max(y, na.rm = TRUE): no non-missing arguments to max; returning
## -Inf
```

Get data and fitted models obtained with my previous implementation (which was shown to give results consistent with Ifat's)...

```r
source("../2016-04-22-glm-for-s-statistic/2016-04-22-glm-for-s-statistic.R")
source("../2016-04-22-glm-for-s-statistic/2016-04-22-glm-for-s-statistic-run.R")
```

The old (left arguments) and new (right arguments) implementation agree perfectly on $N_g$ and $\bar{R}_{i;k\mathcal{G}_j}$.  For instance:

```r
c(identical(d$N_MEST, Y$MEST$N),
identical(d$R_MEST, Y$MEST$R),
identical(d$R_avg8, Y$UA.8$R))
```

```
## [1] TRUE TRUE TRUE
```

But there are slight differences regarding $S_g$ because the new implementation calculates it afresh from $L_g$ and $H_g$ whereas the old implementation imported rounded numbers from Ifat's `pop_skew_3June15.txt` file loosing some precision.

```r
c(identical(d$S_MEST, Y$MEST$S), all.equal(d$S_MEST, Y$MEST$S))
```

```
## [1] FALSE  TRUE
```
Moreover:

```r
c(identical(d$S_avg8, Y$WA.8$S), all.equal(d$S_avg8, Y$WA.8$S))
```

```
## [1] "FALSE"                                
## [2] "Mean relative difference: 0.001112457"
```

```r
c(identical(d$S_avg8, Y$UA.8$S), all.equal(d$S_avg8, Y$UA.8$S))
```

```
## [1] "FALSE"                                
## [2] "Mean relative difference: 0.009172348"
```
which shows that the `S_avg8` statistic corresponds to the `WA.8` weighted average $\bar{S}_{i}$ more closely than to the unweighted `UA.8`.  Analyzing the code of the old implementation in `../2016-04-22-glm-for-s-statistic/2016-04-22-glm-for-s-statistic.R` confirms this.

Note that the names of predictors have been simplified in the new implementation:

```r
str(E[ , 1:15])
```

```
## 'data.frame':	579 obs. of  15 variables:
##  $ Institution  : Factor w/ 3 levels "MSSM","Penn",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ Gender       : Factor w/ 2 levels "Female","Male": 2 2 1 1 2 1 1 2 1 2 ...
##  $ Age          : int  42 58 28 36 52 78 49 62 60 51 ...
##  $ PMI          : num  22.3 19.5 22.8 17.3 22.2 16 20 20.8 24 21.3 ...
##  $ Dx           : Factor w/ 3 levels "AFF","Control",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ RIN          : num  6.9 7 6.9 6.9 7.6 7.4 8 8.4 7.7 7.9 ...
##  $ RIN2         : num  47.6 49 47.6 47.6 57.8 ...
##  $ RNA_lib_batch: Factor w/ 9 levels "0","A","B","C",..: 6 6 5 6 2 7 8 7 6 3 ...
##  $ Ancestry_EV.1: num  0.0214 0.0213 0.0202 0.0213 0.0213 ...
##  $ Ancestry_EV.2: num  0.00459 0.03477 -0.00671 0.02264 -0.00656 ...
##  $ Ancestry_EV.3: num  -0.003252 0.002797 0.000894 0.003056 0.006738 ...
##  $ Ancestry_EV.4: num  0.0381 -0.0202 0.0461 -0.0182 0.041 ...
##  $ Ancestry_EV.5: num  0.000824 -0.00345 -0.005654 -0.007156 -0.013422 ...
##  $ SV1          : num  0.02559 -0.00105 -0.00124 0.02031 -0.02787 ...
##  $ SV2          : num  -0.00651 -0.04842 -0.00523 -0.01151 -0.00416 ...
```

### New functions for fitting


```r
source("~/projects/monoallelic-brain/src/fit-glms.R")
```
Compare estimated coefficients using the new and my old implementation under normal linear model fitted to unweighted average (`UA`) $R$ as response to find perfect agreement:

```r
old <- coef(m$avg8$nlm.R)
new <- coef(do.fit(Y$UA.8$R, X = E, e.v = e.vars, family = gaussian))
all.equal(old, new, tolerance = 0, check.attributes = FALSE)
```

```
## [1] TRUE
```
Using weighted (`WA`) average $S$ as response agrees reasonably but not perfectly since $S$ shows differences due to rounding between old and new importers (as discussed above):

```r
old <- coef(m$avg8$nlm.S)
new <- coef(do.fit(Y$WA.8$S, X = E, e.v = e.vars, family = gaussian))
all.equal(old, new, check.attributes = FALSE)
```

```
## [1] "Mean relative difference: 0.02600831"
```
In case of unweighted (`UA`) average $S$ the agreement is poor as expected based on the results concerning import:

```r
old <- coef(m$avg8$nlm.S)
new <- coef(do.fit(Y$UA.8$S, X = E, e.v = e.vars, family = gaussian))
all.equal(old, new, check.attributes = FALSE)
```

```
## [1] "Mean relative difference: 0.321898"
```
Again with weighted (`WA`) average $S$ as response the agreement is also reasonable under the logistic model:

```r
logi.S <- function(g)
    do.fit(y = cbind(Y[[g]]$H, Y[[g]]$L), X = E, e.v = e.vars, family = binomial)
old <- coef(m$avg8$logi.S)
new <- coef(logi.S("WA.8"))
all.equal(old, new, check.attributes = FALSE)
```

```
## [1] "Mean relative difference: 0.06142334"
```
Using weighted average (`WA`) but subjecting $S$ to an affine transformation $T$ such that $T(S)$ is supported on the interval $[0,1]$ as opposed to $S$'s support on $[1/2,1]$ results in less reasonable agreement between implementations because the "rounding differences" are further amplified:

```r
affine.transform.S <- function(y) {
    H2 <- as.integer((y$S * 2 - 1) * y$N)
    C <- cbind(H2[], y$N - H2[])
    C[ C < 0 & ! is.na(C) ] <- 0
    return(C)
}
logi2.S <- function(g)
    do.fit(y = affine.transform.S(Y[[g]]), X = E, e.v = e.vars, family = binomial)
old <- coef(m$avg8$logi.S)
new <- coef(logi2.S("WA.8"))
all.equal(old, new, check.attributes = FALSE)
```

```
## [1] "Mean relative difference: 0.1552438"
```
