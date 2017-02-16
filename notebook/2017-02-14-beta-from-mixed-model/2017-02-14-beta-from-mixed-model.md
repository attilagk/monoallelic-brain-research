## Motivation

For certain genes $g$ allelic bias was found to depend on biological predictors (Age, Dx, Ancestry) in our previous analysis.  That analysis fitted various **fixed effects** multiple regression models (GLM family) using as response variable the read count ratio $S_{ig}$ or its quasi-log transformed version $Q_{ig}$.  After checking model fit, we primarily based our results to a weighted, normal linear model using $Q_{ig}$, which we called *wnlm.Q* (as opposed to its unweighted variation *unlm.Q*).  Significant dependence of the allelic bias of gene $g$ on (some level of) the predictor $p$ was assessed in terms of the estimated regression coefficient $\hat{\beta}_{pg}$.

But decomposition of variance using `variancePartition` showed that a predictor that significantly affected allelic bias for a certain gene did not necessarily explain substantial fraction of variance.  This was the case in spite of using the same response and predictors in the **mixed effects** model for `variancePartition` as for fixed effects models like wnlm.Q.

What is the reason for this discrepancy?

1. $\hat{\beta}_{pg}$ under the fixed effects model largely deviates from that under the mixed effects model because of
    * differences between the two model frameworks
    * programming bug
1. $\hat{\beta}_{pg}$ under the fixed effects agrees with that under the mixed effects model
    * in this case the discrepancy arises from the somewhat different questions that the two analyses address

The analysis below appears to support the second case rather than the first.

## Calculations


```r
library(variancePartition)
library(doParallel)
library(lattice)
lattice.options(default.args = list(as.table = TRUE))
lattice.options(default.theme = "standard.theme")
opts_chunk$set(dpi = 144)
opts_chunk$set(out.width = "700px")
opts_chunk$set(dev = c("png", "pdf"))
source("../../src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
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

Fit various fixed effect models to data for 30 selected genes


```r
M <- do.all.fits(Z = Y[gene.ids], G = E)
```

Define model formula for fitting using the  `fitVarPartModel`


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

Fit fixed and mixed effects models using `variancePartition` 


```r
M$unlm.Q.varPart <- fitVarPartModel(t(Q), fm$fixed, E)
M$mixed.1 <- fitVarPartModel(t(Q), fm$mixed.1, E)
```

```
## Projected memory usage: > 3.2 Mb 
## Projected run time: ~ 0.03 min
```

```
## Loading required package: Matrix
```

```r
#M$mixed.2 <- fitVarPartModel(t(Q), fm$mixed.2, E)
```

Extract $\hat{\beta}_{pg}$ for each gene $g$ under each model:


```r
coef.names <- names(coef(M$unlm.Q[[1]]))
cf <- data.frame(lapply(M[- length(M)],
                        function(l.m)
                            stack(lapply(l.m,
                                         function(m)
                                             coef(m)[ coef.names ]))$values))
cf$mixed.1 <-
    stack(lapply(M$mixed.1,
                 function(m)
                     unlist(coef(m)[[1]][1, , drop = TRUE])[ coef.names ]))$values
#cf$coefficient <- rep(coef.names, length(gene.ids))
cf$coefficient <- factor(rep(coef.names, length(gene.ids)), ordered = TRUE, levels = coef.names)
cf$gene <- factor(rep(gene.ids, each = length(coef.names)), ordered = TRUE, levels = gene.ids)
```

## Results: estimated coefficients

To test various implementations of the same fixed effects model:

<img src="figure/fixed-fixed-1.png" title="plot of chunk fixed-fixed" alt="plot of chunk fixed-fixed" width="700px" /><img src="figure/fixed-fixed-2.png" title="plot of chunk fixed-fixed" alt="plot of chunk fixed-fixed" width="700px" />

Consistent with the above graphs, the coefficients from one implementation are precisely equal to those from the other implementation:


```r
with(cf, all.equal(unlm.Q, unlm.Q.varPart))
```

```
## [1] TRUE
```

Compare the mixed effects model to the fixed effects model unml.Q; note that both are unweighted

<img src="figure/mixed-1-fixed-1.png" title="plot of chunk mixed-1-fixed" alt="plot of chunk mixed-1-fixed" width="700px" /><img src="figure/mixed-1-fixed-2.png" title="plot of chunk mixed-1-fixed" alt="plot of chunk mixed-1-fixed" width="700px" />

Compare the *unweighted* mixed effects model to *weighted* the fixed effects model wnml.Q

<img src="figure/mixed-1-fixed-weighted-1.png" title="plot of chunk mixed-1-fixed-weighted" alt="plot of chunk mixed-1-fixed-weighted" width="700px" /><img src="figure/mixed-1-fixed-weighted-2.png" title="plot of chunk mixed-1-fixed-weighted" alt="plot of chunk mixed-1-fixed-weighted" width="700px" />
