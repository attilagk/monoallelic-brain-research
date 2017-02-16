## Motivation

For certain genes $g$ allelic bias was found to depend on biological predictors (Age, Dx, Ancestry) in our previous analysis.  That analysis fitted various **fixed effects** multiple regression models (GLM family) using as response variable the read count ratio $S_{ig}$ or its quasi-log transformed version $Q_{ig}$.  After checking model fit, we primarily based our results to a weighted, normal linear model using $Q_{ig}$, which we called *wnlm.Q* (as opposed to its unweighted variation *unlm.Q*).  Significant dependence of the allelic bias of gene $g$ on (some level of) the predictor $p$ was assessed in terms of the estimated regression coefficient $\hat{\beta}_{pg}$.

But decomposition of variance using `variancePartition` showed that a predictor that significantly affected allelic bias for a certain gene did not necessarily explain substantial fraction of variance.  This was the case in spite of using the same response and predictors in the **mixed effects** model for `variancePartition` as for fixed effects models like wnlm.Q.  A notable difference is that the mixed effects model was based on the unweighted unlm.Q rather than the weighted wnlm.Q.

What is the reason for this discrepancy?  This question can be approached based on two cases giving altogether four explanations:

1. $\hat{\beta}_{pg}$ under the fixed effects model largely deviates from that under the mixed effects model because of
    * treating certain predictors' effect random instead of fixed impacts $\hat{\beta}_{pg}$ for other predictors
    * weighting, which was used for the fixed effects model wnlm.Q but not for the mixed effects model
    * programming bug
1. $\hat{\beta}_{pg}$ under the fixed effects agrees with that under the mixed effects model
    * in this case the discrepancy arises from the somewhat different questions that the two analyses address

The analysis below excludes programming bug and appears to suggest that the discrepancy is due to a mixture of the remaining three explanations.

## Calculations


```r
library(variancePartition)
library(doParallel)
library(lattice)
library(latticeExtra)
lattice.options(default.args = list(as.table = TRUE))
lattice.options(default.theme = "standard.theme")
opts_chunk$set(dpi = 144)
opts_chunk$set(out.width = "700px")
opts_chunk$set(dev = c("png", "pdf"))
source("../../src/import-data.R")
source("../../src/graphics.R")
source("../../src/fit-glms.R")
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
# total read count, used as weights
N <- data.frame(sapply(Y[gene.ids], getElement, "N"), check.names = FALSE)
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
fm$fixed.1 <- as.formula(paste(fchar[1], fchar[3]))
fm$mixed.1 <- ~ Age + (1|Institution) + (1|Gender) + PMI + (1|Dx) + RIN +(1|RNA_batch) + Ancestry.1 + Ancestry.2 + Ancestry.3 + Ancestry.4 + Ancestry.5
fm$mixed.2 <- ~ Age + (1|Institution) + Gender + PMI + Dx + RIN +(1|RNA_batch) + Ancestry.1 + Ancestry.2 + Ancestry.3 + Ancestry.4 + Ancestry.5
```

Fit fixed and mixed effects models using `variancePartition` 


```r
M$fixed.1 <- fitVarPartModel(t(Q), fm$fixed.1, E)
M$mixed.1 <- fitVarPartModel(t(Q), fm$mixed.1, E)
```

```
## Projected memory usage: > 3.2 Mb 
## Projected run time: ~ 0.03 min
```

```r
# the next expressions give errors
M$mixed.1.w <- fitVarPartModel(t(Q), fm$mixed.1, E, weightsMatrix = t(N))
```

```
## Projected memory usage: > 3.2 Mb 
## Projected run time: ~ 0.08 min
```

```
## Error in {: task 1 failed - "Invalid grouping factor specification, Institution"
```

```r
M$mixed.2 <- fitVarPartModel(t(Q), fm$mixed.2, E)
```

```
## Projected memory usage: > 3 Mb 
## Projected run time: ~ 0.01 min
```

```
## Error in checkModelStatus(fitInit, showWarnings = showWarnings, colinearityCutoff): Categorical variables modeled as fixed effect: Gender 
## The results will not behave as expected and may be very wrong!!
```

Weights for fitting are expected from `voom`; using total read counts as weights failed with an error for some reason.  Also, there's an error, when some factors are treated to convey fixed while the others as random effects.


```r
# extract variance partitions from model objects
vp <- lapply(M[c("fixed.1", "mixed.1")], function(m) extractVarPart(m))
# reshape into long format
vpl <- data.frame(lapply(vp, function(v) stack(v[ , c(e.vars, "Residuals")])$values))
vpl$predictor <- factor(rep(c(e.vars, "Residuals"), each = length(gene.ids)), ordered = TRUE, levels = c(e.vars, "Residuals"))
vpl$gene <- factor(rep(gene.ids, length(c(e.vars, "Residuals"))), ordered = TRUE, levels = gene.ids)
vpl <- cbind(vpl[c("predictor", "gene")], vpl[- grep("predictor|gene", names(vpl))])
vpll <- reshape(vpl, varying = c("fixed.1", "mixed.1"), v.names = "relative.variance",
                timevar = "model.type", times = c("fixed.1", " mixed.1"), direction = "long")
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
cf <- cbind(cf[c("coefficient", "gene")], cf[- grep("coefficient|gene", names(cf))])
```

Calculate 99% confidence intervals for each $\beta_{pg}$ (also extracts $\hat{\beta}_{pg}$):


```r
beta.hat.CI <- lapply(M[c("unlm.Q", "wnlm.Q")], get.estimate.CI)
```

Finally, load objects for plotting:


```r
source("2017-02-14-beta-from-mixed-model.R")
```

## Results

### Coefficients $\hat{\beta}_{pg}$ compared to variance partitions

Compare the estimated regression coefficients $\{\hat{\beta}_{pg}\}$ to the corresponding (estimated) variance partitions.  Note that the 99% confidence interval for a $\beta_{pg}$ relates to how significantly the response depends on predictor $p$ for gene $g$ in case of continuous predictors (covariates).  For discrete predictors $p$ is a level of the predictor compared to some other level chosen as baseline.

Begin with the "basic" fixed effects model *unlm.Q*, also called *fixed.1* in this document:


```r
my.segplot(beta.hat.CI$unlm.Q, main = expression(paste("99 % CI for ", beta, " under unlm.Q")))
```

<img src="figure/unlm-Q-CI-1.png" title="plot of chunk unlm-Q-CI" alt="plot of chunk unlm-Q-CI" width="700px" />

Now plot variance partitions for *fixed.1* and contrast this with *mixed.1* the mixed effects model derived from unlm.Q by turning all discrete predictors from fixed to random effects:


```r
tp <- trellis.par.get()
my.col <- c(rainbow(length(e.vars))[ c(seq(1, length(e.vars), by = 3), seq(2, length(e.vars), by = 3), seq(3, length(e.vars), by = 3)) ], "gray")
trellis.par.set(superpose.polygon = list(alpha = 1, col = my.col, border = "black", lty = 1, lwd = 0.5))
barchart(gene ~ relative.variance | model.type, data = vpll, groups = predictor, stack = TRUE, auto.key = list(columns = 3), xlim = c(0, 0.6), layout = c(2, 1))
```

<img src="figure/fixed-varpart-1.png" title="plot of chunk fixed-varpart" alt="plot of chunk fixed-varpart" width="700px" />

The next plot illustrates more directly the differences between the fixed and mixed model with respect to variance partitions:


```r
my.main <- "Variance partitions: mixed vs fixed effects model"
my.xlab <- "fixed effects: unlm.Q (fixed.1)"
my.ylab <- "mixed effects: mixed.1"
my.plot(x = "fixed.1", y = "mixed.1", group.by = "predictor", lbl.type = "gene", dt = vpl, main = my.main, xlab = my.xlab, ylab = my.ylab)
```

<img src="figure/fixed-mixed-varpart-1.png" title="plot of chunk fixed-mixed-varpart" alt="plot of chunk fixed-mixed-varpart" width="700px" />

The same results grouped by genes:


```r
update(my.plot(x = "fixed.1", y = "mixed.1", group.by = "gene", lbl.type = "predictor", dt = vpl, subset = vpl$predictor != "Residuals", main = my.main, xlab = my.xlab, ylab = my.ylab), scales = list(relation = "free"))
```

<img src="figure/fixed-mixed-varpart-genes-1.png" title="plot of chunk fixed-mixed-varpart-genes" alt="plot of chunk fixed-mixed-varpart-genes" width="700px" />

Next, return to estimates and confidence intervals for $\beta_{pg}$, but this time under our preferred (based on fit quality) fixed effects model *wnlm.Q*.  The same results are reported in the manuscript.

<img src="figure/wnlm-Q-CI-1.png" title="plot of chunk wnlm-Q-CI" alt="plot of chunk wnlm-Q-CI" width="700px" />

### Fixed vs mixed effects models in terms of $\hat{\beta}$ 

The next set of plots present estimated coefficients $\{\hat{\beta}_{pg}\}$ under fixed effects models fitted with my own scripts and compare those to the corresponding $\{\hat{\beta}_{pg}\}$ under fixed or mixed effects models fitted with  the `variancePartition` package.

#### Fixed effects w/o `variancePartition`

The first tests does not involve any random effects.  Rather it checks various implementations of the same fixed effects model (unlm.Q):

<img src="figure/fixed-fixed-1.png" title="plot of chunk fixed-fixed" alt="plot of chunk fixed-fixed" width="700px" />
<img src="figure/fixed-fixed-coef-1.png" title="plot of chunk fixed-fixed-coef" alt="plot of chunk fixed-fixed-coef" width="700px" />

Consistent with the above graphs, the coefficients from one implementation are precisely equal to those from the other implementation:


```r
with(cf, all.equal(unlm.Q, fixed.1))
```

```
## [1] TRUE
```

#### Mixed vs fixed effects

Compare the mixed effects model *mixed.1* (based on unlm.Q) to the fixed effects model unml.Q; note that both mixed.1 and unlm.Q are unweighted

<img src="figure/mixed-1-fixed-1.png" title="plot of chunk mixed-1-fixed" alt="plot of chunk mixed-1-fixed" width="700px" />

<img src="figure/mixed-1-fixed-coef-1.png" title="plot of chunk mixed-1-fixed-coef" alt="plot of chunk mixed-1-fixed-coef" width="700px" />

#### wlnm.Q: weighted fixed effects model

Compare the (unweighted) mixed effects model mixed.1 to the *weighted* fixed effects model wnml.Q.  Weighting increases the deviation between the two sets of $\hat{\beta}_{pg}$.

<img src="figure/mixed-1-fixed-weighted-1.png" title="plot of chunk mixed-1-fixed-weighted" alt="plot of chunk mixed-1-fixed-weighted" width="700px" />
<img src="figure/mixed-1-fixed-weighted-coef-1.png" title="plot of chunk mixed-1-fixed-weighted-coef" alt="plot of chunk mixed-1-fixed-weighted-coef" width="700px" />

#### Weighted vs unweighted fixed effects

To isolate the impact of weighting:

<img src="figure/fixed-unweighted-weighted-1.png" title="plot of chunk fixed-unweighted-weighted" alt="plot of chunk fixed-unweighted-weighted" width="700px" />
<img src="figure/fixed-unweighted-weighted-coef-1.png" title="plot of chunk fixed-unweighted-weighted-coef" alt="plot of chunk fixed-unweighted-weighted-coef" width="700px" />

