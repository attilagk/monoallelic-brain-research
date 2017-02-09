Load packages...



Import data, define model formula:


```r
# selected set of genes
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
# predictors
names(e.vars) <- e.vars
E <- get.predictors()[e.vars]
# response: Q, quasi-log transformed read count ratio
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
Q <- data.frame(sapply(Y[gene.ids], getElement, "Q"))
# response: R, rank transformed read count ratio
R <- data.frame(sapply(Y[gene.ids], getElement, "R"))
# response: S, untransformed
S <- data.frame(sapply(Y[gene.ids], getElement, "S"))
# model formula
form <- ~ Age + (1|Institution) + (1|Gender) + PMI + (1|Dx) + RIN +(1|RNA_batch) + Ancestry.1 + Ancestry.2 + Ancestry.3 + Ancestry.4 + Ancestry.5
```

Fitting the linear mixed model to various transformations:


```r
vp.Q <- fitExtractVarPartModel( t(Q), form, E )
```

```
## Projected run time: ~ 0.03 min
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

```r
#vp.R <- fitExtractVarPartModel( t(R), form, E )
vp.S <- fitExtractVarPartModel( t(S), form, E )
```

```
## Projected run time: ~ 0.02 min
```

```
## Warning: Some predictor variables are on very different scales: consider
## rescaling

## Warning: Some predictor variables are on very different scales: consider
## rescaling
```

## Results


```r
main <- "Q: quasi-log transformation"
plotVarPart(vp.Q, main = main)
```

<img src="figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" width="700px" />

```r
plotPercentBars(vp.Q)
```

<img src="figure/unnamed-chunk-4-2.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" width="700px" />




```r
main <- "S: untransformed"
plotVarPart(vp.S, main = main)
```

<img src="figure/unnamed-chunk-6-1.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" width="700px" />

```r
plotPercentBars(vp.S)
```

<img src="figure/unnamed-chunk-6-2.png" title="plot of chunk unnamed-chunk-6" alt="plot of chunk unnamed-chunk-6" width="700px" />
