Hello, World!




```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
```

Prepare data frame including all 30 selected genes


```r
dat <- merge.data(gene.ids = gene.ids)
```


```r
fm.1 <- Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) + (scale(Age) | Gene)
fm.3 <- Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) + (1 | Gene:Institution) + (scale(Age) + scale(RIN) + scale(PMI) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + Ancestry.1
```


```r
M1 <- list()
M1$unlm.Q <- lmer(fm.1, data = dat)
M1$wnlm.Q <- lmer(fm.1, data = dat, weights = N)
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model is nearly unidentifiable: very large eigenvalue
##  - Rescale variables?
```

```r
M1$unlm.S <- lmer(reformulate(as.character(fm.1)[3], response = "S"), data = dat)
M1$logi.S <- glmer(reformulate(as.character(fm.1)[3], response = "H.N"), data = dat, family = binomial)
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model failed to converge with max|grad| = 0.00436075 (tol = 0.001, component 1)

## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model is nearly unidentifiable: very large eigenvalue
##  - Rescale variables?
```


```r
get.diagnostic.data <- function(lM, sel.col = 2:6) {
    helper <- function(x)
        cbind(data.frame(Residual = residuals(m <- lM[[x]], type = "deviance"), Fitted.value = predict(m), Family = x),
              model.frame(m)[sel.col])
    l <- lapply(names(lM), helper)
    do.call(rbind, l)
}
```


```r
diag.M1 <- get.diagnostic.data(M1)
qqmath(~ Residual | Family, data = diag.M1) 
```

<img src="figure/unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" width="700px" />
