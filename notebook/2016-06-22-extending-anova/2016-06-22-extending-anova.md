## Introduction

ANOVA and effects were used before...TODO

The least square estimates are given by $\hat{\beta} = (X^T X)^{-1} X^T y$.  Permuting the coefficients is equivalent to permuting the bases of the linear transformation $(X^T X)^{-1} X^T$ or, equivalently, permuting the *rows* of the corresponding matrix or, equivalently again, permuting the *columns* of the design matrix.  Therefore, when explanatory variables are reordered and the columns of the design matrix are permuted accordingly then the least square estimates still remain the same up to their relative order.

Take two permutations of explanatory variables:

1. **forward** has been used in all my and Ifat's previous analysis
2. the **reverse** of the above


## Preparation

Relevant scripts


Import data; note that the set of **selected genes have been updated** based on later analysis

```r
E <- get.predictors() # default arguments
# updated gene set
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
#nobs <- as.data.frame(lapply(list(unfiltered=Y, filtered=Y.f), function(y) sapply(y, function(x) sum(! is.na(x[[1]])))))
```

Fit `wnlm.Q`, `wnlm.R` and `logi.S` using both the forward and the reverse order permutation.

```r
# e.vars defined in fit-glms.R
e.v <- list(forward = e.vars, reverse = rev(e.vars))
# exclude unweighed aggregates UA.8 and UA from fitting
to.fit.ids <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
to.fit.ids <- grep("^WA.8$", to.fit.ids, value = TRUE, invert = TRUE)
M <- lapply(list(unlm.Q = "unlm.Q", wnlm.Q = "wnlm.Q", unlm.R = "unlm.R", wnlm.R = "wnlm.R", logi.S = "logi.S", logi2.S = "logi2.S"),
            function(m) lapply(e.v,
                               function(v) do.all.fits(Y[to.fit.ids], preds = v, sel.models = m)[[1]]))
# a list of (sub)lists mirrioring the structure of M (sublists: forward or reverse)
f.ids <- lapply(M, lapply, function(m) ! sapply(m, is.null))
# the fit for TMEM261P1 has not converged under logi.S
f.ids$logi.S$forward["TMEM261P1"] <- FALSE
f.ids$logi.S$reverse["TMEM261P1"] <- FALSE
```

Consistent with the theory above, wherever the fit converged (with the exception of TMEM261P1), the order has no impact on the regression coefficients.  The next `R` expression compares coefficient estimates between forward and reverse under `logi.S` for every gene or aggregate and reports any difference:

```r
grep("TRUE", sapply(to.fit.ids, function(g) all.equal(coef(M$logi.S$forward[[g]]), coef(M$logi.S$reverse[[g]])[ names(coef(M$logi.S$forward[[g]])) ])), invert = TRUE, value = TRUE)
```

```
## named character(0)
```

## Components of variation

ANOVA allows *direct* comparison of predictors (explanatory variables) by partitioning all explained/systematic variation in the response ($S$ or $R$) among them.  Similarly, the directly comparable effect of each of the $p$ regression coefficients on the response (an $n$-vector) can be obtained from QR decomposition of the response into a set of $n$ orthogonal vectors, in which a subset of $p$ vectors corresponds to the coefficients. 

In the above sentences "directly comparable" means that the effects are on a uniform scale for all predictors/coefficients for a given gene as that scale is closely related to the one on which the response varies.  But the total variation of the response $S$ itself shows variation across genes, which leaves two options in gene-gene comparisons: compare components of variation with or without correction for across gene variation of the total variance.  The first one will be referred below as **genes on uniform scale** and the second as **genes on relative scale**.

### Comparing predictors with genes on uniform scale


```r
A.long <- lapply(lapply(M, lapply, l.anova), reshape.2, type = "anova")
```



Under `logi.S`:

```
## $anova.logi.S
```

<img src="figure/anova-fw-rv-logi-S-1.png" title="plot of chunk anova-fw-rv-logi-S" alt="plot of chunk anova-fw-rv-logi-S" width="700px" />

The same tendencies emerge under `wnlm.Q`:
<img src="figure/anova-fw-rv-wnlm-Q-1.png" title="plot of chunk anova-fw-rv-wnlm-Q" alt="plot of chunk anova-fw-rv-wnlm-Q" width="700px" />

Similarly under `wnlm.R`:
<img src="figure/anova-fw-rv-wnlm-R-1.png" title="plot of chunk anova-fw-rv-wnlm-R" alt="plot of chunk anova-fw-rv-wnlm-R" width="700px" />


```r
Ef.long <- lapply(M, function(m) { x <- l.l.effects(m); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
```

Under `logi.S`:
<img src="figure/effects-fw-rv-logi-S-1.png" title="plot of chunk effects-fw-rv-logi-S" alt="plot of chunk effects-fw-rv-logi-S" width="700px" />

Again, similar tendencies are observed under `wnlm.Q`:
<img src="figure/effects-fw-rv-wnlm-Q-1.png" title="plot of chunk effects-fw-rv-wnlm-Q" alt="plot of chunk effects-fw-rv-wnlm-Q" width="700px" />

Similarly under `wnlm.R`:
<img src="figure/effects-fw-rv-wnlm-R-1.png" title="plot of chunk effects-fw-rv-wnlm-R" alt="plot of chunk effects-fw-rv-wnlm-R" width="700px" />

### Figure for manuscript

<img src="figure/anova-effects-fw-rv-logi-S-1.png" title="plot of chunk anova-effects-fw-rv-logi-S" alt="plot of chunk anova-effects-fw-rv-logi-S" width="700" />

### Comparison with genes on relative scale

<img src="figure/anova-fw-rv-logi-S-trellis-1.png" title="plot of chunk anova-fw-rv-logi-S-trellis" alt="plot of chunk anova-fw-rv-logi-S-trellis" width="700px" />

<img src="figure/effects-fw-rv-logi-S-trellis-1.png" title="plot of chunk effects-fw-rv-logi-S-trellis" alt="plot of chunk effects-fw-rv-logi-S-trellis" width="700px" />

<img src="figure/anova-fw-rv-wnlm-Q-trellis-1.png" title="plot of chunk anova-fw-rv-wnlm-Q-trellis" alt="plot of chunk anova-fw-rv-wnlm-Q-trellis" width="700px" />

<img src="figure/effects-fw-rv-wnlm-Q-trellis-1.png" title="plot of chunk effects-fw-rv-wnlm-Q-trellis" alt="plot of chunk effects-fw-rv-wnlm-Q-trellis" width="700px" />

<img src="figure/anova-fw-rv-wnlm-R-trellis-1.png" title="plot of chunk anova-fw-rv-wnlm-R-trellis" alt="plot of chunk anova-fw-rv-wnlm-R-trellis" width="700px" />

<img src="figure/effects-fw-rv-wnlm-R-trellis-1.png" title="plot of chunk effects-fw-rv-wnlm-R-trellis" alt="plot of chunk effects-fw-rv-wnlm-R-trellis" width="700px" />

### Another view, genes on uniform scale

This display is meant to be comparable to a similarly structured trellis display of estimated regression coefficients.

#### Under logi.S

<img src="figure/effects-fw-rv-logi-S-trellis-coef-cond-1.png" title="plot of chunk effects-fw-rv-logi-S-trellis-coef-cond" alt="plot of chunk effects-fw-rv-logi-S-trellis-coef-cond" width="700px" />

#### Under logi2.S

<img src="figure/effects-fw-rv-logi2-S-trellis-coef-cond-1.png" title="plot of chunk effects-fw-rv-logi2-S-trellis-coef-cond" alt="plot of chunk effects-fw-rv-logi2-S-trellis-coef-cond" width="700px" />

#### Under wnlm.Q

<img src="figure/effects-fw-rv-wnlm-Q-trellis-coef-cond-1.png" title="plot of chunk effects-fw-rv-wnlm-Q-trellis-coef-cond" alt="plot of chunk effects-fw-rv-wnlm-Q-trellis-coef-cond" width="700px" />

#### Under unlm.Q

<img src="figure/effects-fw-rv-unlm-Q-trellis-coef-cond-1.png" title="plot of chunk effects-fw-rv-unlm-Q-trellis-coef-cond" alt="plot of chunk effects-fw-rv-unlm-Q-trellis-coef-cond" width="700px" />

#### Under wnlm.R

<img src="figure/effects-fw-rv-wnlm-R-trellis-coef-cond-1.png" title="plot of chunk effects-fw-rv-wnlm-R-trellis-coef-cond" alt="plot of chunk effects-fw-rv-wnlm-R-trellis-coef-cond" width="700px" />

#### Under unlm.R

<img src="figure/effects-fw-rv-unlm-R-trellis-coef-cond-1.png" title="plot of chunk effects-fw-rv-unlm-R-trellis-coef-cond" alt="plot of chunk effects-fw-rv-unlm-R-trellis-coef-cond" width="700px" />

## Estimate and CI for regression coefficients

This and the next section do not deal with analysis of variance specifically but rather with the estimated regression coefficient $\beta_{jg}$ for each column $X_j$ of the design matrix and for each gene $g$.  But because of computational convenience these are presented here.


```r
Betas <- lapply(M, function(m) { x <- get.estimate.CI(m$forward); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
```


```r
my.segplot(data = Betas$logi.S, xlim = my.xlim)
```

<img src="figure/reg-coef-logi-S-1.png" title="plot of chunk reg-coef-logi-S" alt="plot of chunk reg-coef-logi-S" width="700px" />

<img src="figure/reg-coef-logi2-S-1.png" title="plot of chunk reg-coef-logi2-S" alt="plot of chunk reg-coef-logi2-S" width="700px" />

<img src="figure/reg-coef-wnlm-Q-1.png" title="plot of chunk reg-coef-wnlm-Q" alt="plot of chunk reg-coef-wnlm-Q" width="700px" />

<img src="figure/reg-coef-unlm-Q-1.png" title="plot of chunk reg-coef-unlm-Q" alt="plot of chunk reg-coef-unlm-Q" width="700px" />

<img src="figure/reg-coef-wnlm-R-1.png" title="plot of chunk reg-coef-wnlm-R" alt="plot of chunk reg-coef-wnlm-R" width="700px" />

<img src="figure/reg-coef-unlm-R-1.png" title="plot of chunk reg-coef-unlm-R" alt="plot of chunk reg-coef-unlm-R" width="700px" />

## Comparison of models

The pairwise model comparisons in terms of $\hat{\beta}$ under either of two selected models are meant to assess sensitivity of the results to various aspects of modeling:

* logistic vs normal linear model
* scaling of the logit link function
* data transformations for the normal linear model: $Q$ vs $R$ statistic
* weighting for the normal linear model

Each panel in the plot shows the theoretical zero line $\beta_{jg}=0$ under each of the two models (horizontal and vertical dashed lines).  The plotting symbols are color coded according to gene rank (rainbow, red to violet); the last "rank" #31 (violet) corresponds to the weighted average `WA` of read count ratio over genes.  The plotting symbols also display the rank with numbers, see the key on the top.  Genes acceptably fitted by both models are distinguished with a diamond symbol and **bold font** from those that could be fitted only by wnlm.Q.



Filtering out results for logi.S if the fit is bad for a given gene, using decisions based on model checking analysis


```r
logi.S.OK <- read.csv("../../results/model-checking.csv", row.names = "gene")["logi.S.fit.OK"]
# filtered and unfiltered long format Betas
Betas.l.f <- Betas.l <- do.call(cbind, Betas)
# perform filtering by replacing data with NAs
Betas.l.f[Betas.l.f$logi.S.Gene %in% c(rownames(logi.S.OK)[! logi.S.OK$logi.S.fit.OK], "WA"), grep("logi\\.S\\.[ELU]", names(Betas.l.f))] <- NA
```

### Logistic vs normal linear model

The plot shows overall agreement between logi.S and wnlm.Q.  Genes of disagreement tend to be those for which logi.S fits poorly to the data.

Without filtering genes with poor fit by logi.S


```r
mtype.compare.plot(mtypeA = "logi.S", mtypeB = "wnlm.Q", dt = Betas.l, do.key = TRUE)
```

<img src="figure/logi-S-wnlm-Q-compare-1.png" title="plot of chunk logi-S-wnlm-Q-compare" alt="plot of chunk logi-S-wnlm-Q-compare" width="700px" />

With filtering


```r
mtype.compare.plot(mtypeA = "logi.S", mtypeB = "wnlm.Q", dt = Betas.l.f, do.key = TRUE)
```

<img src="figure/logi-S-filtered-wnlm-Q-compare-1.png" title="plot of chunk logi-S-filtered-wnlm-Q-compare" alt="plot of chunk logi-S-filtered-wnlm-Q-compare" width="700px" />

### Scaling of the logit link function

There is very little impact on the $2\times$ difference in scaling of the logit link function because most of the observed cases are near the upper bound of the link function (which is 1), where the scaling has the smallest effect on the predictions.

Without filtering genes with poor fit by logi.S

<img src="figure/logi-S-logi2-S-compare-1.png" title="plot of chunk logi-S-logi2-S-compare" alt="plot of chunk logi-S-logi2-S-compare" width="700px" />

With filtering

<img src="figure/logi-S-filtered-logi2-S-compare-1.png" title="plot of chunk logi-S-filtered-logi2-S-compare" alt="plot of chunk logi-S-filtered-logi2-S-compare" width="700px" />

### Data transformations for the normal linear model: $Q$ vs $R$

Overall good agreement.

<img src="figure/wnlm-R-wnlm-Q-compare-1.png" title="plot of chunk wnlm-R-wnlm-Q-compare" alt="plot of chunk wnlm-R-wnlm-Q-compare" width="700px" />

### Weighting for the normal linear model

Overall good agreement, suggesting relatively small impact of weighting.

<img src="figure/unlm-Q-wnlm-Q-compare-1.png" title="plot of chunk unlm-Q-wnlm-Q-compare" alt="plot of chunk unlm-Q-wnlm-Q-compare" width="700px" />
