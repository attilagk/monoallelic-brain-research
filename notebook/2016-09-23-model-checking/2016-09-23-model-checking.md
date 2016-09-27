

## Preliminaries

Load functions


```r
source("~/projects/monoallelic-brain/src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
source("2016-09-23-model-checking.R")
```


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
```

Get data

```r
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
```

### Checking known outliers

Import outliers of the CMC study provided by Mette and Jessica and check if any has been included in our previous analyses.  The answer is: they have all been excluded.


```r
read.csv("../../data/outliers.csv", colClasses = "character")$DissectionID %in% row.names(E)
```

```
##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [12] FALSE
```

## Evaluation of model fit

Below is a visual evaluation of model fit using all 6 model types (unlm.S, wnlm.S, unlm.R, wnlm.R, logi.S, logi2.S) on each of the 30 data sets that represent the 30 selected genes for regression analysis.  This visual evaluation is based on three kind of plots, which check the following aspects of model fit:

1. normality of residuals
1. homogieneity of error (homoscedasticity)
1. influence of individual cases

These points will be detailed in the respective sections.

Start by fitting all models:


```r
# exclude unweighed aggregates UA.8 and UA from fitting
to.fit.ids <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
#sel.models <- c("logi.S", "wnlm.R"); names(sel.models) <- sel.models
sel.models <- NULL
M <- do.all.fits(Y[to.fit.ids], G = E, preds = e.vars, sel.models = sel.models)
```

Calculate diagnostics:


```r
diagnostics <- lapply(names(M), function(m) get.diagnostics(M[[m]][gene.ids], mtype = m))
diagnostics <- do.call(rbind, diagnostics)
```

### Normality of residuals

Two kind of standardized residual is plotted against the quantiles of the standard normal distribution:

1. standardized **d**eviance residual, $r_{\mathrm{d}i}$
1. a combination $r^\ast_i$ of the standardized deviance and **P**earson residuals (the latter denoted as $r_{\mathrm{p}i}$)

The combination is defined as in Statistical Models (A.C. Davison), p477: $r^\ast_i = r_{\mathrm{d}i} + r_{\mathrm{d}i}^{-1} \log (r_{\mathrm{p}i} / r_{\mathrm{d}i})$.  Note that $r_{\mathrm{d}i}$, $r_{\mathrm{d}i}$ and $r^\ast_i$ differ for generalized linear models such as logi.S and logi2.S but for linear models they all equal the usual definition of residual i.e. the observed minus fitted value $y_i - \bar{y}_i$.


```r
myqqnorm <- function(mtype)
    qqmath(~ res.std.dev + res.combined | gene, data = diagnostics,
           ylim = c(-4, 4), abline = c(0, 1), pch = "+", subset = model.type %in% mtype,
           main = mtype, xlab = "normal quantiles", ylab = "standardized residual")
myqqnorm("unlm.S")
```

<img src="figure/qqnorm-unlm-S-1.png" title="plot of chunk qqnorm-unlm-S" alt="plot of chunk qqnorm-unlm-S" width="700px" />

<img src="figure/qqnorm-wnlm-S-1.png" title="plot of chunk qqnorm-wnlm-S" alt="plot of chunk qqnorm-wnlm-S" width="700px" />

<img src="figure/qqnorm-unlm-R-1.png" title="plot of chunk qqnorm-unlm-R" alt="plot of chunk qqnorm-unlm-R" width="700px" />

<img src="figure/qqnorm-wnlm-R-1.png" title="plot of chunk qqnorm-wnlm-R" alt="plot of chunk qqnorm-wnlm-R" width="700px" />

<img src="figure/qqnorm-logi-S-1.png" title="plot of chunk qqnorm-logi-S" alt="plot of chunk qqnorm-logi-S" width="700px" />

<img src="figure/qqnorm-logi2-S-1.png" title="plot of chunk qqnorm-logi2-S" alt="plot of chunk qqnorm-logi2-S" width="700px" />

### Homogeneity of error (homoscedasticity)

$\sqrt{r_{\mathrm{d}i}}$ is plotted against the fitted value $\bar{y}_i$.  The black lines correspond to the data passed through a LOESS filter.  All models assume homogeneous error that is no systematic variation with the fitted value. 


```r
myhomoscedas <- function(mtype, xlim = c(0.65, 1.05))
    xyplot(sqrt(abs(res.std.dev)) ~ fitted | gene, data = diagnostics, xlim = xlim,
           panel = function(...) {
               panel.xyplot(...)
               panel.smoother(..., method = "loess", col = "black", se = FALSE)
           },
           subset = model.type %in% mtype, pch = "+", main = mtype, xlab = "fitted value",
           ylab = expression(sqrt(std.deviance.resid)))

myhomoscedas("unlm.S")
```

<img src="figure/homoscedas-unlm-S-1.png" title="plot of chunk homoscedas-unlm-S" alt="plot of chunk homoscedas-unlm-S" width="700px" />

<img src="figure/homoscedas-wnlm-S-1.png" title="plot of chunk homoscedas-wnlm-S" alt="plot of chunk homoscedas-wnlm-S" width="700px" />

<img src="figure/homoscedas-unlm-R-1.png" title="plot of chunk homoscedas-unlm-R" alt="plot of chunk homoscedas-unlm-R" width="700px" />

<img src="figure/homoscedas-wnlm-R-1.png" title="plot of chunk homoscedas-wnlm-R" alt="plot of chunk homoscedas-wnlm-R" width="700px" />

<img src="figure/homoscedas-logi-S-1.png" title="plot of chunk homoscedas-logi-S" alt="plot of chunk homoscedas-logi-S" width="700px" />

<img src="figure/homoscedas-logi2-S-1.png" title="plot of chunk homoscedas-logi2-S" alt="plot of chunk homoscedas-logi2-S" width="700px" />

### Influence of individual cases

Here Cook's distance $C_i$ is plotted against $h_{ii} / (1 - h_{ii})$, where the leverage $h_{ii}$ of case (individual) $i$ is the $i$-th element of the diagonal of the (weighted) projection matrix $H$.  The definition is

$$
\begin{equation}
C_i = r_{\mathrm{p}i}^2 \frac{h_{ii}}{p (1 - h_{ii})}
\end{equation}
$$
where $p$ is the number of parameters (see Statistical Models, p477).

So Cook's distance combines the influence of case $i$ based on both the response (i.e. the residual $r_{\mathrm{p}i}$) and  the predictors (more precisely the leverage $h_{ii}$, which for lm only depends on the predictors and for glm depends also on the response to some degree).  Good fit is indicated by equally low $C_i$ (and $h_{ii}$) for all cases implying that all cases carry equal amount of information on the regression parameters $\beta$.


```r
myinfluence <- function(mtype, xlim = c(-0.2, 4), ylim = c(-0.04, 0.8))
    xyplot(cooks.dist ~ leverage / (1 - leverage) | gene, data = diagnostics, xlim = xlim, ylim = ylim,
           subset = model.type %in% mtype, pch = "+", main = mtype, ylab = "Cook's distance")

myinfluence("unlm.S")
```

<img src="figure/influence-unlm-S-1.png" title="plot of chunk influence-unlm-S" alt="plot of chunk influence-unlm-S" width="700px" />

<img src="figure/influence-wnlm-S-1.png" title="plot of chunk influence-wnlm-S" alt="plot of chunk influence-wnlm-S" width="700px" />

<img src="figure/influence-unlm-R-1.png" title="plot of chunk influence-unlm-R" alt="plot of chunk influence-unlm-R" width="700px" />

<img src="figure/influence-wnlm-R-1.png" title="plot of chunk influence-wnlm-R" alt="plot of chunk influence-wnlm-R" width="700px" />

<img src="figure/influence-logi-S-1.png" title="plot of chunk influence-logi-S" alt="plot of chunk influence-logi-S" width="700px" />

<img src="figure/influence-logi2-S-1.png" title="plot of chunk influence-logi2-S" alt="plot of chunk influence-logi2-S" width="700px" />

### Identifying outliers

Here the median Cook's distance $C_i$ taken across all 30 genes (with possibly missing values) for any given case/individual $i$ is used to quantify the overall impact of that case.  The focus is on logi.S, since that model seems in general the most useful in making inferences.  The top 20 cases are shown; especially the top 2 appear to be extreme outliers, whose exclusion might substantially improve model fit for some genes.


```r
cooks.d <- lapply(M, function(m)
                  as.matrix(data.frame(lapply(m[gene.ids], get.influence, rownames(E)))))
med.cooks.d <- lapply(cooks.d,
                      function(d) sort(apply(d, 1, median, na.rm = TRUE), decreasing = TRUE))
dotplot(rev(med.cooks.d$logi.S[1:20]), main = "logi.S", xlab = "median Cook's distance")
```

<img src="figure/cooks-dist-1.png" title="plot of chunk cooks-dist" alt="plot of chunk cooks-dist" width="700px" />