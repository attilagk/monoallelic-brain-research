
```
## Loading required package: ggplot2
```

```
## Loading required package: foreach
```

```
## foreach: simple, scalable parallel programming from Revolution Analytics
## Use Revolution R for scalability, fault tolerance and more.
## http://www.revolutionanalytics.com
```

```
## Loading required package: Biobase
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, cbind, colnames,
##     do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, lengths, Map, mapply,
##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
##     Position, rank, rbind, Reduce, rownames, sapply, setdiff,
##     sort, table, tapply, union, unique, unsplit, which, which.max,
##     which.min
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: iterators
```

Import data


```r
# selected set of genes
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
# predictors
names(e.vars) <- e.vars
E <- get.predictors()[e.vars]
# response: Q, quasi-log transformed read count ratio
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
Q <- data.frame(sapply(Y[gene.ids], getElement, "Q"), check.names = FALSE)
# response: S, untransformed
S <- data.frame(sapply(Y[gene.ids], getElement, "S"), check.names = FALSE)
```

Fixed effect models


```r
M <- do.all.fits(Z = Y[gene.ids], G = E)
```

Define model formula:


```r
(fchar <- as.character(formula(M$unlm.Q$MEST)))
```

```
## [1] "~"                                                                                                                       
## [2] "Y"                                                                                                                       
## [3] "Age + Institution + Gender + PMI + Dx + RIN + RNA_batch + Ancestry.1 + Ancestry.2 + Ancestry.3 + Ancestry.4 + Ancestry.5"
```

```r
fm <- list()
fm$fixed <- as.formula(paste(fchar[1], fchar[3]))
fm$mixed.1 <- ~ Age + (1|Institution) + (1|Gender) + PMI + (1|Dx) + RIN +(1|RNA_batch) + Ancestry.1 + Ancestry.2 + Ancestry.3 + Ancestry.4 + Ancestry.5
fm$mixed.2 <- ~ Age + (1|Institution) + Gender + PMI + Dx + RIN +(1|RNA_batch) + Ancestry.1 + Ancestry.2 + Ancestry.3 + Ancestry.4 + Ancestry.5
```


```r
M$unlm.Q.varPart <- fitVarPartModel(t(Q), fm$fixed, E)
M$mixed.1 <- fitVarPartModel(t(Q), fm$mixed.1, E)
```

```
## Projected memory usage: > 3.2 Mb 
## Projected run time: ~ 0.02 min
```

```
## Loading required package: Matrix
```

```r
#M$mixed.2 <- fitVarPartModel(t(Q), fm$mixed.2, E)
```

Choose a gene randomly


```r
set.seed(1968)
(gene <- sample(gene.ids, size = 1))
```

```
##   RP11-909M7.3 
## "RP11-909M7.3"
```

Comparison


```r
l <- lapply(M[c("unlm.Q", "unlm.Q.varPart")], function(l.m) coef(l.m[[gene]]))
all.equal(l[[1]], l[[2]])
```

```
## [1] TRUE
```


```r
xyplot(y <- coef(M$unlm.Q.varPart[[gene]]) ~ coef(M$unlm.Q.varPart[[gene]]),
                   xlab = "fixed effects (unlm.Q)", ylab = "fixed effects (calling fitVarPartModel",
                   main = "Two implementations of fixed")
```

<img src="figure/unnamed-chunk-8-1.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" width="700px" />

```r
y <- unlist(coef(M$mixed.1[[gene]])[[1]][1, , drop = TRUE])
xyplot(y ~ coef(M$unlm.Q[[gene]])[names(y)], xlab = "fixed effects (unlm.Q)", ylab = "mixed effects (based on unlm.Q)",
       main = "Mixed vs fixed")
```

<img src="figure/unnamed-chunk-8-2.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" width="700px" />

```r
xyplot(y ~ coef(M$wnlm.Q[[gene]])[names(y)], xlab = "fixed effects (wnlm.Q)", ylab = "mixed effects (based on unlm.Q)",
       main = "Mixed vs weighted fixed")
```

<img src="figure/unnamed-chunk-8-3.png" title="plot of chunk unnamed-chunk-8" alt="plot of chunk unnamed-chunk-8" width="700px" />
