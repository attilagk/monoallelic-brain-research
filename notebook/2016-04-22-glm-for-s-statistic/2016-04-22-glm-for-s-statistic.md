## Introduction

## Preparations

### Data

Data files

```r
files <- list(S="pop_skew_3June15.txt",
              N="pop_cov_3June15.txt",
              X="DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv",
              X2="samples.csv")
```

Selected gene sets of size $8, 13, 16$, sequentially nested in each other

```r
             # 8 genes analyzed by Ifat
genes <- c("PEG3", "INPP5F", "SNRPN", "PWAR6", "ZDBF2", "MEG3", "ZNF331", "GRB10",
             # 5 more genes analyzed by AGK 3/2/16
             "PEG10", "SNHG14", "NAP1L5", "KCNQ1OT1", "MEST",
             # 3 more genes present in data files
             "IGF2", "NLRP2", "UBE3A")
```

Response variables, characterized by

1. *filtering*: whether the $S_{ig}$ statistic is excluded for individual $i$ and gene $g$ if the total read count $N_{ig}\le 50$
1. *transformation*: whether $S_{ig}$ is transformed into a rank $R_{ig}$ for a given $g$ across all individuals $i$
1. *aggregation*: whether $S_{ig}$ or $R_{ig}$ are aggregated across genes based on gene sets of size $8, 13, 16$

Explanatory variables $X$; only those that are included in the linear predictor (see below).

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

### Regression models

Three regression models are considered.  All three fit into the generalized linear model (glm) framework, and characterized by

1. the *linear predictor* $\eta = \sum_r x_r \beta_r$, where $r$ indexes explanatory variables $x_r$ and the corresponding regression *coefficients* $\beta_r$ mediating *effects*
1. the *link function*: a one to one mapping of $\eta$ onto the mean response $\mu\equiv \mathrm{E} Y$
1. $Y_i|\eta_i$, the *conditional distribution of the response* given the predictor for observation (individual) $i$

Logistic models are part of the generalized linear family, which use linear predictors but a non-linear link function and non-normal conditional distribution of the response given the predictors.  For logistic models the link function is the sigmoid-shaped logistic function and the conditional distribution is binomial, whose denominator is used as weight.  In the present case the denominator is the total read count $n_{ig}$, whereas the response is $S_{ig}$ without any transformation.  Thus, $S_{ig}$ is assumed to be binomially distributed.

### Implementation

Ifat's earlier implementation and my **reimplementation** of data filtering, transformation and fitting procedures, which was necessary for implementing newer regression models and to facilitate reproducibility.

```r
source('2016-04-22-glm-for-s-statistic.R')
source('../2016-03-02-ifats-regression-analysis/2016-03-02-ifats-regression-analysis.R')
```

```r
m.ifat.g8 <- fit.lm(transform.data(genes[1:13], genes[1:8]), do.thrs=FALSE)
```


The new script provides four implementations for three regression models.

```r
f <- list(
          nlm = function(y, d) { # normal linear model with rank R as response
              glm(formula = mk.form(paste0("R_", y)), family = gaussian, data = d)
          } ,
          nlm2 = function(y, d) { # as above but with the R's lm function instead of glm
              lm(formula = mk.form(paste0("R_", y)), data = d)
          } ,
          logi = function(y, d) { # logistic regression with the S statistic as response
              glm(formula = mk.form(paste0("S_", y)), family = binomial, data = d, weights = d[[paste0("N_", y)]])
          } ,
          logi2 = function(y, d) { # as above but with rescaled and offset logistic link function
              glm(formula = mk.form(paste0("S_", y, " * 2 - 1")), family = binomial, data = d, weights = d[[paste0("N_", y)]])
          } )
```

## Results

### Comparing implementations

Comparing the normal linear model obtained with Ifat's implementation to that with the new implementation (see `nlm` function above).

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
These results show a close but not perfect match between Ifat's and my implementation (the following numbers are given up to 4 significant digits).  The discrepancy is due to my omission of a rounding step and only slightly affects the results.

### Effect of filtering

The difference is larger (but still small) when we omit a filtering step that excludes $S_{ig}$ values for which the total read count $n_{ig}\le 50$ for a given individual $i$ and gene $g$ (third column).  For instance, the residual deviance (residual sum of squares in case of NLM), Akaike information criterion (AIC) and regression coefficient for `Age.of.Death` between are:

### Comparing models

### Effects on aggregated responses (gene sets)

### Variability of effects on individual genes

Fit a NLM for each gene separately (after transformation).  It is clear that there is much variability among genes, since for some of them (e.g. ZDBF2) the regression coefficient for age is more significantly different from 0 than when $\mathrm{LOI\_R}$ was based on 8 genes (see results above).  On the other hand for some other genes, (e.g. PEG3) the difference is far from significant.

```r
#nlm.loir <- lapply(genes, do.fit, family=gaussian, weights=as.data.frame(rep(1, nrow(tot.read.n))), filter.thrs=50)
#print(sapply(nlm.loir, function(x) summary(x)$coefficients["Age.of.Death", 4]), digits=3)
##print(sapply(nlm.loir, function(x) summary(x)$aic), digits=3)
```

## Conclusion
