## Preparations

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
```

Manipulations of $S$ and $N$ give rise to $response$ variables $Y=[y_{ig}]$ of
various regression models (below). These manipulations fall in the following
categories:

1. *filtering*: whether the $s_{ig}$ statistic is excluded for individual $i$ and gene $g$ if the total read count $N_{ig}\le 50$
1. *transformation*: whether $s_{ig}$ is transformed into a rank $r_{ig}$ for a given $g$ across all individuals $i$
1. *aggregation*: whether $s_{ig}$ or $r_{ig}$ are aggregated across genes based on gene sets of size $8, 13, 16$

Ifat's earlier work used only one set of manipulations: filtering + rank transformation + aggregation by averaging across the set of $8$ genes.  The resulting response variable was termed $\mathrm{LOI\_R}$ for "loss of imprinting ratio".  This term is avoided here because imprinting (or allelic imbalance) may increase with age for certain genes or sets of genes, and also because $r_{ig}$ is more concise than $\mathrm{LOI\_R}_{ig}$.

The explanatory variables in $X$ are

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

Although, in principle, the models could be freely combined with various data filtering, transformation and aggregation, the focus here is to use the best fitting combinations of models and transformations.  Therefore the optimal response is $R_{ig}$ for nlm (i.e. transformation) and $S_{ig}$ for logi and logi2 (no transformation).

For the logi and logi2 models the response, for each observation (individual) is distributed binomially and the denominator is used as weight for that observation.  Using $S_{ig}$ as response suggests using the corresponding observed $n_{ig}$ as weight.

The table summarizes the models' properties.

|                   |    nlm            |    logi           |       logi2       |
|:-----------------:|:-----------------:|:-----------------:|:-----------------:|
|    link function  |    linear         |    logistic       |  scaled logistic  |
|response distrib.  |normal (Gaussian)  |   binomial        |     binomial      |
| optimal response  |  $R_{ig}$         |   $S_{ig}$        |   $S_{ig}$        |
|   weights         |      1 (uniform)  |   $n_{ig}$        |   $n_{ig}$        |

### Implementation

The multifaceted goals of the present analysis called for a completely *new implementation* of data manipulations and models because the earlier implementation (by Ifat) lacked the necessary modularity.  The script files for the new and old implementation:

```r
source('2016-04-22-glm-for-s-statistic.R')
source('../2016-03-02-ifats-regression-analysis/2016-03-02-ifats-regression-analysis.R')
```


The new script provides the function `nlm` for the nlm model, `logi` for logi, and `logi2` for logi2.  Additionally, the `nlm2` function also implements nlm, but instead of calling `R`'s `glm` function (which is also called by `logi` and `logi2`), it calls `lm`.   The definition of `logi2` is nearly identical to  `logi`; the only difference is using `C2` as response variable instead of `C`.  Wheres `C` contains the original observed "higher" and "lower" read counts, `C2` is a transformation of those that corresponds to the (inverse) scaling function of the logistic functiion in logi2.


```r
f <- list(
          nlm = function(y, d) { # normal linear model with rank R as response
              glm(formula = mk.form(paste0("R_", y)), family = gaussian, data = d)
          } ,
          nlm2 = function(y, d) { # as above but with the R's lm function instead of glm
              lm(formula = mk.form(paste0("R_", y)), data = d)
          } ,
          logi = function(y, d) { # logistic regression with the S statistic as response
              glm(formula = mk.form(paste0("C_", y)), family = binomial, data = d)
          } ,
          logi2 = function(y, d) { # as above but with rescaled and offset logistic link function
              glm(formula = mk.form(paste0("C2_", y)), family = binomial, data = d)
          } )
```
All results with `nlm2` were fournd to be identical to `nlm` (not shown) so the `nlm2`-related results will not be presented furthermore.

## Results

### Comparing implementations

Let's compare the normal linear model obtained with Ifat's implementation to that with the new implementation!  For consistency with Ifat's [previous results][ifat] the data manipulations are: filtering, rank-transformation, and averaging over $8$ selected genes.

```r
sapply( c("deviance", "aic", "coefficients"), function(s) all.equal(m.ifat.g8[[s]], m$g8$nlm[[s]]))
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

The table lists the relative change in AIC induced by omission of filtering, defined as $(\mathrm{AIC}_{f=1} - \mathrm{AIC}_{f=0}) / \mathrm{AIC}_{f=1}$, where $f=1$ indicates filtering.

```
##               nlm      logi     logi2
## PEG3     -0.01880 -3.52e-03 -3.53e-03
## INPP5F   -0.00204 -3.60e-04 -4.62e-04
## SNRPN    -0.22500 -4.35e-02 -4.85e-02
## PWAR6    -0.00558 -6.90e-05 -7.54e-05
## ZDBF2    -0.12900 -3.73e-02 -3.70e-02
## MEG3     -0.10800 -2.23e-02 -2.21e-02
## ZNF331   -0.35300 -2.86e-01 -3.66e-01
## GRB10    -1.11000 -7.95e-01 -9.46e-01
## PEG10    -0.08210 -3.75e-02 -3.93e-02
## SNHG14   -0.12300 -3.07e-02 -3.10e-02
## NAP1L5   -0.22300 -1.88e-01 -2.43e-01
## KCNQ1OT1 -1.01000 -5.37e-01 -6.48e-01
## MEST     -0.44400 -2.85e-01 -3.08e-01
## IGF2     -3.55000 -8.59e+00 -1.54e+01
## NLRP2    -5.43000 -3.74e+00 -5.71e+00
## UBE3A    -1.97000 -3.03e+00 -4.77e+00
## g8        0.01430 -9.99e-05 -2.17e-04
## g13       0.01710 -3.75e-03 -4.34e-03
## g16       0.02860  1.43e-02  1.59e-02
```
Filtering has small effect in most cases.  The exceptions are those genes for which filtering removed many points such as NLRP2 and IGF2 (169 and 149 points removed, respectively) in contrast with genes like PEG3 (only 8 points removed).

All results below were obtained with filtering.

### Comparing models

Below are relative AIC of logi compared to nlm, defined as $(\mathrm{AIC}_{\mathrm{logi}} - \mathrm{AIC}_{\mathrm{nlm}}) / \mathrm{AIC}_{\mathrm{logi}}$.  Note that the response is $S_{ig}$ for logi and $R_{ig}$ for nlm.

```
##     PEG3   INPP5F    SNRPN    PWAR6    ZDBF2     MEG3   ZNF331    GRB10 
##  0.18100  0.37500 -0.04590  0.00951 -0.25100 -0.05240 -0.32300 -1.11000 
##    PEG10   SNHG14   NAP1L5 KCNQ1OT1     MEST     IGF2    NLRP2    UBE3A 
## -0.36700 -0.60400 -1.20000 -1.20000 -0.65000  8.90000 -0.21600  9.09000 
##       g8      g13      g16 
##  0.72200  0.72500  0.72800
```
For 12 out of 16 single genes nlm performs better than logi.  But for all
aggregated (averaged) gene sets logi clearly performs better because it takes a
weighted average of genes using the total read counts $n_{ig}$, which nlm
ignores.

The relative AIC of logi2 compared to logi:

```
##     PEG3   INPP5F    SNRPN    PWAR6    ZDBF2     MEG3   ZNF331    GRB10 
##   0.4470   0.4470   0.4580   0.4650   0.4610   0.4610   0.3800   0.3210 
##    PEG10   SNHG14   NAP1L5 KCNQ1OT1     MEST     IGF2    NLRP2    UBE3A 
##   0.4090   0.4370   0.3800   0.3760   0.3200   0.0521   0.3080   0.0670 
##       g8      g13      g16 
##   0.4750   0.4700   0.4700
```
This shows that logi2 fits all data sets better than logi altough for genes with sparse data 

The figure compares all three model families using (the absolute) AIC
![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

logi2 clearly outperforms logi and nlm when fitted to the any of the three aggregated datasets.  This tendency holds for individual genes as well, at least for genes with ample data (large number $n_{ig}$ of total read counts).  For genes with sparse data (e.g. NLRP2 and IGF2 mentioned before) nlm and logi2 fit similarly well; but for these genes cases AIC is much smaller (under any given model) underscoring the sparsity of data and weakening conclusions.

Both logi and logi2 benefit from data aggregation to a great degree.  This is attributable to the regression weights gained from total read counts, which the unweighted nlm ignores.

## Conclusion


[nlm g13]: {% post_url 2016-03-02-ifats-regression-analysis %}
[ifat]: https://docs.google.com/presentation/d/1YvpA1AJ-zzir1Iw0F25tO9x8gkSAzqaO4fjB7K3zBhE/
