## Preparations


```r
source('2016-04-22-glm-for-s-statistic.R')
source('../2016-03-02-ifats-regression-analysis/2016-03-02-ifats-regression-analysis.R')
```


### Data

Data files:

```r
files <- list(S="pop_skew_3June15.txt",
              N="pop_cov_3June15.txt",
              X="DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv",
              X2="samples.csv")
```

File `S` contains the matrix $S=[s_{ig}]$ of the observed $s_{ig}$ statistic for all individuals $i$ and selected genes $g$. (The capitalized $S_{ig}$ denotes the corresponding random variable).  Similarly, `N` contains a matrix of total read counts $N_{ig}$, where "total" reflects summing over both alleles.  `X` contains the design matrix $X=[X_{ir}]$, where $r$ indexes *explanatory variables*. `X2` contains a subset of explanatory variables found in `X` but, unlike `X`, `X2` contains `RNAseq_ID`s that are used for mapping from $X$ to $S$ (and to $N$) and are absent from `X2`.

Selected gene sets of size $8, 13, 16$, sequentially nested in each other:

```r
             # 8 genes analyzed by Ifat
genes <- c("PEG3", "INPP5F", "SNRPN", "PWAR6", "ZDBF2", "MEG3", "ZNF331", "GRB10",
             # 5 more genes analyzed by AGK 3/2/16
             "PEG10", "SNHG14", "NAP1L5", "KCNQ1OT1", "MEST",
             # 3 more genes present in data files
             "IGF2", "NLRP2", "UBE3A")
# response variables
genes.or.gsets <- c(genes, "avg8", "avg13", "avg16", "pool8", "pool13", "pool16")
names(genes.or.gsets) <- genes.or.gsets
```

Manipulations of $S$ and $N$ give rise to $response$ variables $Y=[y_{ig}]$ of
various regression models (below). These manipulations fall in the following
categories:

1. *filtering*: whether the $s_{ig}$ statistic is excluded for individual $i$ and gene $g$ if the total read count $N_{ig}\le 50$
1. *transformation*: whether $s_{ig}$ is transformed into a rank $r_{ig}$ for a given $g$ across all individuals $i$
1. *aggregation*: whether $s_{ig}$ or $r_{ig}$ are aggregated across genes based on gene sets of size $8, 13, 16$; aggregation may be achieved via
    1. *averaging* responses across genes
    1. *pooling* responses from multiple genes

Ifat's earlier work used only one set of manipulations: filtering + rank transformation + aggregation by averaging across the set of $8$ genes.  The resulting response variable was termed $\mathrm{LOI\_R}$ for "loss of imprinting ratio".  This term is avoided here because imprinting (or allelic imbalance) may increase with age for certain genes or sets of genes, and also because $r_{ig}$ is more concise than $\mathrm{LOI\_R}_{ig}$.

The explanatory variables $x_r \in X$ are

```r
expl.var <- c("`Age.of.Death`",
               "Institution",
               "Gender",
               "`PMI..in.hours.`",
               "Dx",
               "`DLPFC_RNA_isolation..RIN`", "`DLPFC_RNA_isolation..RIN.2`",
               "`DLPFC_RNA_report..Clustered.Library.Batch`",
               "`Ancestry.EV.1`", "`Ancestry.EV.2`", "`Ancestry.EV.3`", "`Ancestry.EV.4`", "`Ancestry.EV.5`" )
```
Note that there are additional variables in `X` but those are not included in the models, following Ifat's earlier work.

### Regression models

Three regression models are considered here: **nlm** (normal linear), **logi** (logistic) and **logi2** (scaled logistic). All three fit into the *generalized linear model* (glm) framework, and
characterized by

1. the *linear predictor* $\eta = \sum_r x_r \beta_r$, where $\beta_r$ are regression *coefficients* mediating the *effects* of $X$ on the response
1. the *link function*: a one to one mapping of $\eta$ onto the mean response $\mu\equiv \mathrm{E} Y$
1. $P(Y_i|\eta_i)$, the *conditional distribution of the response* given the predictor for observation (individual) $i$

For the logi and logi2 models the response, for each observation (individual) is distributed binomially and the denominator is used as weight for that observation.  Using $S_{ig}$ as response suggests using the corresponding observed $n_{ig}$ as weights.  In contrast the nlm model has linear link function and normal (Gaussian) response distribution; under this model the weight is uniformly 1 across observations (individuals) for single and pooled genes.  For averaged genes the weights are, as previously, also taken as uniformly 1.  In this case, however, a more rigorous treatment would be defining the weight at observation $i$ as the number $|G|_i$ of averaged genes at that $i$ since that number varies due to non-uniformly missing data and filtering.

The link functions are given by the following equations
$$
\begin{equation}
\mu = \eta \qquad \text{linear; nlm}
\end{equation}
$$
$$
\begin{equation}
\mu = \frac{1}{1 + e^{-\eta}} \qquad \text{logistic; logi}
\end{equation}
$$
$$
\begin{equation}
\mu = \frac{1}{2} + \frac{1}{2 (1 + e^{-\eta})} \qquad \text{scaled logistic; logi2}
\end{equation}
$$

These functions are illustrated by the following plot; nlm: solid green, logi: solid red, and logi2: dashed blue.  They reflect simple (univariate) regression models with $S_{i\mathrm{PEG3}}$ as response and `Age.of.Death` as the only explanatory variable.  Thus there are two regression coefficients in all three cases, and these were estimated with iterative weighted least squares using `R`'s `glm` function.  logi and logi2 are almost indistinguishable in the range of all observed ages but become quite different around 300 years (much longer than the human lifespan).

![plot of chunk s-stat-simple-regr](figure/s-stat-simple-regr-1.png) 

Below printed are two measures of fit: residual deviance (top block) and AIC (bottom block).


```
##        nlm.R       nlm2.R        nlm.S       logi.S      logi2.S 
## 4.348327e+05 4.348327e+05 2.028023e-01 5.484768e+03 1.127678e+04
```

```
## $nlm.R
## [1] 4740.081
## 
## $nlm2.R
## NULL
## 
## $nlm.S
## [1] -2432.414
## 
## $logi.S
## [1] 7116.834
## 
## $logi2.S
## [1] 13207.78
```

Although the models might appear freely combinable with both data transformations, it is not clear how logistic model(s) might be fitted after rank transformation of $S_{ig}$.  Because how should total read counts $n_{ig}$ (the weights) be transformed?  But the converse case is straight-forward: to fit nlm to untransformed $S_{ig}$ despite its apparent heteroscedasticity and nonlinearity (see [Ifat's plots][ifat] and the data exploration below).

Therefore the three model families give rise to four types of models depending on whether the response is rank-transformed or not.  These are summarized in the table below:

|                   |    nlm.R          |    nlm.S          |    logi.S         |       logi2.S     |
|:-----------------:|:-----------------:|:-----------------:|:-----------------:|:-----------------:|
|    link function  |    linear         |    linear         |    logistic       |  scaled logistic  |
|response distrib.  |normal (Gaussian)  |normal (Gaussian)  |   binomial        |     binomial      |
|         response  |  $R_{ig}$         |          $S_{ig}$ |   $S_{ig}$        |   $S_{ig}$        |
|   weights         |      1 (uniform)  |      1 (uniform)  |   $n_{ig}$        |   $n_{ig}$        |

So, considering these 4 types combined with 16 genes separately, or with the two kinds of data aggregations (averaging and pooling) using 3 gene sets, each with or without filtering, there are $4 \times (16 + 2 \times 3) \times 2 = 176$ fitted models in the present analysis.

### Implementation

The multifaceted goals of the present analysis called for a completely *new implementation* of data manipulations and models because the earlier implementation (by Ifat) lacked the necessary modularity.  The script files for the new and old implementation:

```r
source('2016-04-22-glm-for-s-statistic.R')
source('../2016-03-02-ifats-regression-analysis/2016-03-02-ifats-regression-analysis.R')
```

The new code provides the function `nlm` for the nlm model, `logi` for logi, and `logi2` for logi2.  Additionally, the `nlm2` function also implements nlm, but instead of calling `R`'s `glm` function (which is also called by `logi` and `logi2`), it calls `lm`.   The definition of `logi2` is nearly identical to  `logi`; the only difference is using `C2` as response variable instead of `C`.  Wheres `C` contains the original observed "higher" and "lower" read counts, `C2` is a transformation of those that corresponds to the (inverse) scaling function of the logistic functiion in logi2.


```r
f <- list(
          nlm.R = function(y, d) { # normal linear model with rank R as response
              glm(formula = mk.form(paste0("R_", y)), family = gaussian, data = d)
          } ,
          nlm2.R = function(y, d) { # as above but with the R's lm function instead of glm
              lm(formula = mk.form(paste0("R_", y)), data = d)
          } ,
          nlm.S = function(y, d) { # normal linear model with S statistic as response
              glm(formula = mk.form(paste0("S_", y)), family = gaussian, data = d)
          } ,
          nlm2.S = function(y, d) { # as above but with the R's lm function instead of glm
              lm(formula = mk.form(paste0("S_", y)), data = d)
          } ,
          logi.S = function(y, d) { # logistic regression with the S statistic as response
              glm(formula = mk.form(paste0("C_", y)), family = binomial, data = d)
          } ,
          logi2.S = function(y, d) { # as above but with rescaled and offset logistic link function
              glm(formula = mk.form(paste0("C2_", y)), family = binomial, data = d)
          } )
```
As expected, all results with `nlm2` were found to be identical to `nlm` (not shown).

## Exploratory analysis

The following plots explore the relationship between a given response variable and 12 of the 13 explanatory variables (the 13th one, Ancestry.EV.5 showed no systematic relationship and is omitted here so that plots can be arranged in $3 \times 4$ arrays).  In this document only filtered data are presented because filtering had no visually appreciable effect on the plot.

### Responses averaged across genes

When the response is the average $\bar{S}_{i}=\left( \sum_g n_g \right)^{-1} \sum_{g=1}S_{ig} n_{ig}$ (below), a qualitatively similar, inverse, relationship emerges between $\bar{S}$ statistic and age as seen in [Ifat's plots][ifat], which depicted 13 genes separately.   Some of the remaining 11 explanatory variables, like Institution or PMI.in.hours seem to effect $\bar{S}$, while others like Gender don't.

![plot of chunk s-avg16-all-expl-var](figure/s-avg16-all-expl-var-1.png) 

When ranks are averaged as $\bar{R}_i=\sum_gR_{ig}$, the response appears more homoscedastic (its dispersion appears unaffected by explanatory variables) but some of the systematic effects that are seen without transformation (e.g. PMI.in.hours) are diminished.

![plot of chunk r-avg16-all-expl-var](figure/r-avg16-all-expl-var-1.png) 

### Responses pooled across genes

The next set of plots shows data points for all genes pooled together rather than averaged together so each plot has many more points than before.  Below are plots with $s_{ig}$ as response for $g\in$ all 16 genes.

![plot of chunk s-pool16-all-expl-var](figure/s-pool16-all-expl-var-1.png) 

Below are plots with $r_{ig}$ (observed ranks) as response.  Compared to the corresponding results obtained with averaging there is clearly less systematic variation of the response with explanatory variables, which hints at their differential effects on various genes.

![plot of chunk r-pool16-all-expl-var](figure/r-pool16-all-expl-var-1.png) 

### Response vs age for each gene separately

The differential effect of age on the response of various genes is illustrated now.
Plotted below is the observed $s_{ig}$ statistic versus age for each gene $g$ separately.  Dots are data points and the solid line is the smoothed data with the Lowess filter.  The lower percentile of $s_{ig}$, for each $g$, has been trimmed off to enhance clarity.  Judged from the smooth curves the effect of age appears quite small.  Still, gene-to-gene variation is apparent.

![plot of chunk s-stat-age-12genes](figure/s-stat-age-12genes-1.png) 

Below are similar plots to the ones above but now with $R$ as response (i.e. with rank-transformation).  The among genes variation is quite clear.

![plot of chunk r-stat-age-12genes](figure/r-stat-age-12genes-1.png) 

### Conditional dependence on age given other explanatory variables

The apparent dependence of the response on age is marginal in the sense that other explanatory variables are disregarded in the previous plots.  However, those other variables may induce spurious dependence between the response and age even if those are conditionally independent.  Therefore, the next two sets of plots are conditioned on one of three explanatory variables, `Institution`, the RNA quality measure `DLPFC_RNA_isolation..RIN`, and `Ancestry.EV.2` reporting on the population structure(?) that were found to exert highly significant effect in Ifat's earlier regression analysis.  `DLPFC_RNA_isolation..RIN.2` had also highly significant effect, but it correlates so tightly with `DLPFC_RNA_isolation..RIN` (Pearson corr. coef 1) that it carries no additional information.

Importantly, these plots below suggest that, at least for PEG3, the response's dependence on age seems genuine and not entirely due to confounding effects of other variables.

![plot of chunk coplot-r-peg3](figure/coplot-r-peg3-1.png) 

```
## 
##  Missing rows: 8, 13, 16, 17, 32, 54, 57, 59, 60, 66, 67, 68, 77, 79, 84, 100, 116, 122, 131, 142, 150, 161, 177, 189, 192, 194, 197, 205, 208, 212, 216, 219, 226, 234, 246, 251, 252, 264, 277, 279, 284, 285, 289, 300, 302, 311, 312, 321, 322, 339, 341, 347, 349, 350, 358, 359, 366, 369, 372, 376, 377, 384, 394, 409, 418, 423, 429, 432, 439, 445, 446, 451, 472, 475, 480, 490, 492, 508, 512, 513, 519, 520, 524, 532, 536, 554, 575
```

![plot of chunk coplot-r-peg3](figure/coplot-r-peg3-2.png) 

```
## 
##  Missing rows: 8, 13, 16, 17, 32, 54, 57, 59, 60, 66, 67, 68, 77, 79, 84, 100, 116, 122, 131, 142, 150, 161, 177, 189, 192, 194, 197, 205, 208, 212, 216, 219, 226, 234, 246, 251, 252, 264, 277, 279, 284, 285, 289, 300, 302, 311, 312, 321, 322, 339, 341, 347, 349, 350, 358, 359, 366, 369, 372, 376, 377, 384, 394, 409, 418, 423, 429, 432, 439, 445, 446, 451, 472, 475, 480, 490, 492, 508, 512, 513, 519, 520, 524, 532, 536, 554, 575
```

![plot of chunk coplot-r-peg3](figure/coplot-r-peg3-3.png) 

```
## 
##  Missing rows: 8, 13, 16, 17, 32, 54, 57, 59, 60, 66, 67, 68, 77, 79, 84, 100, 116, 122, 131, 142, 150, 161, 177, 189, 192, 194, 197, 205, 208, 212, 216, 219, 226, 234, 246, 251, 252, 264, 277, 279, 284, 285, 289, 300, 302, 311, 312, 321, 322, 339, 341, 347, 349, 350, 358, 359, 366, 369, 372, 376, 377, 384, 394, 409, 418, 423, 429, 432, 439, 445, 446, 451, 472, 475, 480, 490, 492, 508, 512, 513, 519, 520, 524, 532, 536, 554, 575
```

## Regression analysis

### Comparing implementations

Let's compare the normal linear model obtained with Ifat's implementation to that with the new implementation!  For consistency with Ifat's [previous results][ifat] the data manipulations are: filtering, rank-transformation, and averaging over $8$ selected genes.

```r
sapply( c("deviance", "aic", "coefficients"), function(s) all.equal(m.ifat.avg8[[s]], m$avg8$nlm.R[[s]]))
```

```
##                                deviance 
##  "Mean relative difference: 0.04045548" 
##                                     aic 
## "Mean relative difference: 0.004622527" 
##                            coefficients 
##  "Mean relative difference: 0.02769083"
```
These results show a close but not perfect match between the old and new implementation.  The slight discrepancy seems to be due to a rounding step in the old code, for which I found no justification and therefore omitted from the new implementation.

### Effect of filtering

Filtering out "noisy" observations (those with low total read counts) is expected to improve model fit to data.  When too few points are retained after filtering, overfitting may occur.  This is illustrated by the plots below, where filtering induces extremely large improvement for genes with scarce data, like IGF2 or NLRP2.  Improvement is expressed as *minus relative AIC* defined as
$- (\mathrm{AIC}_{f1} - \mathrm{AIC}_{f0}) / \mathrm{AIC}_{f1}$, where $f1$ indicates filtering and $f0$ no filtering.

![plot of chunk s-stat-filtering](figure/s-stat-filtering-1.png) 

All results presented in this document were obtained with filtering unless otherwise noted.

### To what extent is the apparent effect of age confounded?

The previous conditional plots of $R_{i\mathrm{PEG3}}$ vs age qualitatively examined the extent of confounding by three additional variables such as Institution.  The regression coefficient $\beta_{\mathrm{age}}$ under some model reports the effect of age on the response.  In *simple regression* $\beta_{\mathrm{age}}$ is estimated from a model where the linear predictor contains only age so that $\eta = \beta_0 + x_\mathrm{age} \beta_{\mathrm{age}}$.  Then the additional variables $x_r, \; r=2,...$ behave as latent variables confounding the effect $\beta_\mathrm{age}$ of age.  *Multiple regression* includes these additional variables in the linear predictor removes their confounding effect from $\beta_\mathrm{age}$.

Here I compare $\beta_\mathrm{age}$ estimated from simple and multiple regression under the nlm.S, logi.S and logi2.S model types.  First let's return to the earlier plot of fitted simple regression models to the $S_{i\mathrm{PEG3}}$ and $x_\mathrm{age}$ data.  Given a model family such as nlm, the mean $\mu$ of $S_{i\mathrm{PEG3}}$ is a function of $x_\mathrm{age}$ and as well as $\beta_0, \beta_\mathrm{age}$ estimated from simple regression, shown by the left graph.  This function is changed when $\beta_0, \beta_\mathrm{age}$ are estimated from multiple regression together with 21 other regression coefficients $\beta_r$ (right panel).  For instance, $\beta_\mathrm{RIN}$ and $\beta_\mathrm{RIN2}$ mediate the effects of RNA quality measures $x_\mathrm{RIN}$ and $x_\mathrm{RIN2}$.  Two fitted curves are shown for each model: the solid lines correspond to the case when $x_\mathrm{RIN}=\bar{x}_\mathrm{RIN}$ and $x_\mathrm{RIN2}=\bar{x}_\mathrm{RIN2}$ that is the average of the observed $x_{i\mathrm{RIN}}$ and $x_{i\mathrm{RIN2}}$ across all individuals $i$, respectively; on the other hand, the dashed lines correspond to $x_\mathrm{RIN}=\mathrm{min}_i x_{i \mathrm{RIN}}$ and $x_\mathrm{RIN2}=\mathrm{min}_i x_{i \mathrm{RIN2}}$.

The main point, though, is that for any given model the $\beta_\mathrm{age}$ estimated from multiple regression is less negative in that from simple regression.  This appears visually as shallower fitted curves.  The interpretation is that a large confounding effect by other variables was removed by multiple regression.

![plot of chunk s-stat-cmp-simple-multiple-regr-peg3](figure/s-stat-cmp-simple-multiple-regr-peg3-1.png) 

The changes in the fitted curves are qualitatively similar for all genes and aggregated gene sets.  ZNF331 shows the most profound effect of age on $S_{ig}$ after the removal of the confounding effects by multiple regression.
![plot of chunk s-stat-cmp-simple-multiple-regr-znf331](figure/s-stat-cmp-simple-multiple-regr-znf331-1.png) 

### Estimates of the effect of age and confounding

The following set of plots compares estimates of $\beta_\mathrm{age}$ from the simple (dark grey) and multiple (light grey) regression under various model families.  These results can be summarized as

1. age has a negative effect on the response (loss of imprinting)
1. based on multiple regression, the estimated effect size and significance of age is smaller, which indicates substantial confounding in simple regression 
1. for genes with scarce data the effect size can be large without being significantly different from zero

#### nlm.R

![plot of chunk r-stat-beta-age-nlm](figure/r-stat-beta-age-nlm-1.png) 

#### nlm.S

![plot of chunk s-stat-beta-age-nlm](figure/s-stat-beta-age-nlm-1.png) 

#### logi.S

![plot of chunk s-stat-beta-age-logi](figure/s-stat-beta-age-logi-1.png) 

#### logi2.S

![plot of chunk s-stat-beta-age-logi2](figure/s-stat-beta-age-logi2-1.png) 

### Comparing models using AIC

The figure compares all three model families using AIC; nlm is used both
![plot of chunk s-stat-regr-auc](figure/s-stat-regr-auc-1.png) 

## Conclusion


[ifat]: https://docs.google.com/presentation/d/1YvpA1AJ-zzir1Iw0F25tO9x8gkSAzqaO4fjB7K3zBhE/
