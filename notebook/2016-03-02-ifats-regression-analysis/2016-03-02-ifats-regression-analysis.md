

## Goal

Reproduce and extend Ifat's previous regression analysis (documented [here][regression-slide], see also speaker notes for that slide) using the following sets of genes:

* "8 genes that looked promising on the 2X2 table fisher exact test"
* "all 13 genes from slide 3"


```r
genes13 <- c("PEG3", "INPP5F", "SNRPN", "PWAR6", "ZDBF2", "MEG3", "ZNF331",
           "GRB10", "PEG10", "SNHG14", "NAP1L5", "KCNQ1OT1", "MEST" )
```

The question is the **sensitivity** of the results of regression analysis to the gene set $G$.

## Ifat's script

To aid reproducibility, I turned Ifat's code into two functions (`transform.data` and `fit.lm`) stored in the source file below.

```r
source("2016-03-02-ifats-regression-analysis.R")
```

### Input files

Here I denote the test statitic for monoallelic expression as $S_{ig}$ for individual $i$ and gene $g$.  The corresponding matrix is $S = [ S_{ig}  ]$.  Based on Ifat's code, a data transformation (explained below) on $S$ gives rise to $\mathrm{LOI\_R}$.  $\mathrm{LOI\_R}$ plays the role of the response variable in the regression analysis.  The explanatory variables are arranged in the design matrix $X = [ x_{ij} ]$ for all individuals $i$ and variable type $j$.

|     content                                |       file                         |
|:-------------------------------------------|:-----------------------------------|
|   $S$                                      |   `pop_skew_3June15.txt`           |
| subject info, part of $X$                  |   `samples.csv`                    |
|  ???, part of $X$    |`DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv`|         


### Mathematical formulation

The `R` command for linear regression in Ifat's code is

```r
glm(LOI_R ~ (`Age.of.Death` > age.thrs) + Institution + ... + `Ancestry.EV.5`, 
    data=FULL)
```
which corresponds to the normal linear model
$$
\begin{equation} \mathrm{LOI\_R} = X \beta + \epsilon \end{equation}
$$
where $\epsilon$ is a vector of independent normal variables each with mean 0 and variance $\sigma^2$.

How is $S$ transformed into $\mathrm{LOI\_R}$?
Let's inspect Ifat's code, in which the `S2` matrix variable corresponds to the matrix $S$ of response variables, as above.

```r
# rank individuals for each gene
for (g in seq_along(genes)){
  gene=genes[g]
  gene_data<-as.numeric(round(S2[,gene],3))
  gene_rank<-rank(gene_data,ties.method = "min",na.last = "keep") # do the ranking
  R[,gene]<-gene_rank
  #...
  rank_max<-max(R[,gene], na.rm = TRUE)
}
#...
# average rank over genes for each individual
for (s in 1:579){
#...
  PR[s,gene]<-as.numeric(round(R[s,gene]/rank_max*100,0))
  #...
  LOI_R[s,1]<-as.numeric(sprintf("%.0f",mean(PR[s,],na.rm = TRUE)))
}
```

Let $m$ be the number of individuals and $G$ the set of selected genes.  Then
the above code chunk means the definition $\mathrm{LOI\_R}=100 \times m^{-1} \times (\bar{R}_1,...,\bar{R}_m)$, where
$$
\begin{eqnarray}
\bar{R}_i &=& |G|^{-1} \sum_{g\in{G}} R_{ig}
\end{eqnarray}
$$

### Consequence on confidence intervals

Let  $\beta_j$ denote a regression parameter of interest for some explanatory variable $x_j$, say for the age of death. Let $\hat{\beta}_j$ be the least squares (i.e. maximum likelihood) *estimate* of $\beta_j$.  We want to test the null hypothesis that $\beta_j=0$, which we interpret as no influence of age on allelic exclusion. To test that hypothesis at significance level $\alpha$, we need to calculate $\mathrm{CI}_{1-\alpha}$, the $1-\alpha$ confidence interval for $\beta_j$ and see if $0\in\mathrm{CI}_{1-\alpha}$, in which case the null hypothesis is supported. The [calculation][CI-wiki] of $\mathrm{CI}_{1-\alpha}$ follows from normal distribution theory and [weighted linear least squares][weighted], which together imply that the statistic
$$
\begin{equation}
T = \frac{\hat{\beta}_j - \beta_j}{\sqrt{w^{-1}S^2(X^\top X)^{-1}_{jj}}}
\end{equation}
$$
has Student's t-distribution with $m - p$ degrees of freedom, where $m$ is the number of individuals and $p$ is the number of parameters (one more than the number of explanatory variables).  The denominator is the *standard error* for $\hat{\beta}_j$ and includes the following terms:

* $w$, a uniform weight on each of the $m$ observations
* $S^2$, the residual sum of squares divided by $m-p$, which is an unbiased estimator for $\sigma^2$
* the inverse covariance matrix based on the design matrix $X$

Weight $w$ is central here.  Since for each individual $i$ the average rank $\bar{R}_i$ is the average of $|G|$ number of genes, all $m$ rank observations receive weight $w=|G|$.

With observed $S^2=s^2$ the confidence interval is $\mathrm{CI}_{1-\alpha} = \hat{\beta}_j \pm  t(\alpha/2) \sqrt{|G|^{-1}s^2(X^\top X)^{-1}_{jj}}$, where $t(q)$ is the $q$ quantile of the t-distribution with $m - p$ degrees of freedom.  Thus, averaging across all genes in the selected set $G$ shrinks the standard error by a factor of $|G|^{-1/2}$ and consequently shrinks the $1-\alpha$ confidence interval by the same factor, which is about 0.35 for the 8 gene set and 0.28 for the 13 gene set.  The shrinkage in turn leads to more significant p-values.

But the `R` code chunk for linear regression (above) shows that the default $w=1$ was used mistakenly instead of $w=|G|$ and thus the shrinkage effect of averaging $|G|$ genes was ignored.  Consequently, the mistake has lead to a $\sqrt{|G|}$-fold **overestimation** of the standard error and of the $1-\alpha$ confidence interval, and ultimately to a **less significant** p-value.

## Reanalysis

The design matrix $X$ has actually two versions: one where age is a continuous variable and another one where age is a binary variable indicating whether an individual deceased before or after a threshold (set to 70 years).  The following commands carry out the regression for both versions of $X$ and
both gene sets (the set of 8 and that of 13 genes):


```r
m8.thrs <- fit.lm(transform.data(genes13, genes13[1:8]), do.thrs=TRUE, age.thrs=70)
m13.thrs <- fit.lm(transform.data(genes13, genes13), do.thrs=TRUE, age.thrs=70)
m8 <- fit.lm(transform.data(genes13, genes13[1:8]), do.thrs=FALSE)
m13 <- fit.lm(transform.data(genes13, genes13), do.thrs=FALSE)
```

## Results

The detailed, quantitative, results are below.
Qualitatively, they show that using 13 instead of only 8 genes...

1. has little impact on parameter estimates,
2. slightly decreases most standard errors,
3. but slightly weakens significance (increases p-values).

The last two points together suggest that age and other explanatory variables tend to impact the originally selected 8 genes more strongly than the remaining 5 genes.

### With age threshold (70 years)

#### 8 out of 13 genes


```
## 
## Call:
## glm(formula = LOI_R ~ (Age.of.Death > age.thrs) + Institution + 
##     Gender + PMI..in.hours. + Dx + DLPFC_RNA_isolation..RIN + 
##     DLPFC_RNA_isolation..RIN.2 + DLPFC_RNA_report..Clustered.Library.Batch + 
##     Ancestry.EV.1 + Ancestry.EV.2 + Ancestry.EV.3 + Ancestry.EV.4 + 
##     Ancestry.EV.5, data = FULL)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -54.103  -10.955    1.312   10.935   49.660  
## 
## Coefficients:
##                                              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                                -120.98120   48.29967  -2.505 0.012537 *  
## Age.of.Death > age.thrsTRUE                  -8.48248    1.83862  -4.613 4.92e-06 ***
## InstitutionPenn                             -14.66759    2.21989  -6.607 9.19e-11 ***
## InstitutionPitt                              17.07766    2.32461   7.346 7.31e-13 ***
## GenderMale                                    1.21119    1.57096   0.771 0.441042    
## PMI..in.hours.                               -0.20652    0.07179  -2.877 0.004173 ** 
## DxControl                                     2.00790    2.77643   0.723 0.469864    
## DxSCZ                                         2.21546    2.84674   0.778 0.436757    
## DLPFC_RNA_isolation..RIN                     49.24825   13.11683   3.755 0.000192 ***
## DLPFC_RNA_isolation..RIN.2                   -3.19817    0.88394  -3.618 0.000324 ***
## DLPFC_RNA_report..Clustered.Library.BatchA   -6.37843    3.97267  -1.606 0.108936    
## DLPFC_RNA_report..Clustered.Library.BatchB   -9.21359    3.70588  -2.486 0.013204 *  
## DLPFC_RNA_report..Clustered.Library.BatchC   -5.35732    4.09303  -1.309 0.191115    
## DLPFC_RNA_report..Clustered.Library.BatchD   -8.87429    3.73932  -2.373 0.017973 *  
## DLPFC_RNA_report..Clustered.Library.BatchE   -5.43433    3.74031  -1.453 0.146815    
## DLPFC_RNA_report..Clustered.Library.BatchF   -4.62088    3.79535  -1.218 0.223927    
## DLPFC_RNA_report..Clustered.Library.BatchG   -0.79537    5.10595  -0.156 0.876269    
## DLPFC_RNA_report..Clustered.Library.BatchH   -9.60252    5.75671  -1.668 0.095869 .  
## Ancestry.EV.1                                -9.14623   19.19145  -0.477 0.633850    
## Ancestry.EV.2                                 6.64675   21.02653   0.316 0.752036    
## Ancestry.EV.3                               -13.42012   18.47021  -0.727 0.467788    
## Ancestry.EV.4                                32.20447   19.02029   1.693 0.090986 .  
## Ancestry.EV.5                                -2.52111   17.94634  -0.140 0.888331    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 292.3244)
## 
##     Null deviance: 273860  on 577  degrees of freedom
## Residual deviance: 162240  on 555  degrees of freedom
##   (1 observation deleted due to missingness)
## AIC: 4946.6
## 
## Number of Fisher Scoring iterations: 2
```

#### All 13 genes


```
## 
## Call:
## glm(formula = LOI_R ~ (Age.of.Death > age.thrs) + Institution + 
##     Gender + PMI..in.hours. + Dx + DLPFC_RNA_isolation..RIN + 
##     DLPFC_RNA_isolation..RIN.2 + DLPFC_RNA_report..Clustered.Library.Batch + 
##     Ancestry.EV.1 + Ancestry.EV.2 + Ancestry.EV.3 + Ancestry.EV.4 + 
##     Ancestry.EV.5, data = FULL)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -47.090   -8.634    1.034    9.652   43.776  
## 
## Coefficients:
##                                             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                                -94.60394   41.86086  -2.260 0.024210 *  
## Age.of.Death > age.thrsTRUE                 -5.70543    1.59352  -3.580 0.000373 ***
## InstitutionPenn                            -14.92868    1.92396  -7.759 4.12e-14 ***
## InstitutionPitt                             14.73272    2.01472   7.313 9.21e-13 ***
## GenderMale                                  -0.89002    1.36154  -0.654 0.513582    
## PMI..in.hours.                              -0.11965    0.06222  -1.923 0.054987 .  
## DxControl                                    2.02784    2.40630   0.843 0.399747    
## DxSCZ                                        2.95694    2.46724   1.198 0.231241    
## DLPFC_RNA_isolation..RIN                    42.44456   11.36823   3.734 0.000208 ***
## DLPFC_RNA_isolation..RIN.2                  -2.74199    0.76610  -3.579 0.000375 ***
## DLPFC_RNA_report..Clustered.Library.BatchA  -7.54408    3.44308  -2.191 0.028861 *  
## DLPFC_RNA_report..Clustered.Library.BatchB -10.63588    3.21185  -3.311 0.000988 ***
## DLPFC_RNA_report..Clustered.Library.BatchC  -6.82454    3.54739  -1.924 0.054889 .  
## DLPFC_RNA_report..Clustered.Library.BatchD -10.05955    3.24083  -3.104 0.002007 ** 
## DLPFC_RNA_report..Clustered.Library.BatchE  -6.88943    3.24169  -2.125 0.034006 *  
## DLPFC_RNA_report..Clustered.Library.BatchF  -6.49657    3.28940  -1.975 0.048763 *  
## DLPFC_RNA_report..Clustered.Library.BatchG  -2.27714    4.42528  -0.515 0.607055    
## DLPFC_RNA_report..Clustered.Library.BatchH -10.97661    4.98928  -2.200 0.028216 *  
## Ancestry.EV.1                              -11.23119   16.63305  -0.675 0.499809    
## Ancestry.EV.2                               -2.80328   18.22349  -0.154 0.877801    
## Ancestry.EV.3                               -3.01643   16.00795  -0.188 0.850606    
## Ancestry.EV.4                               34.65121   16.48471   2.102 0.036001 *  
## Ancestry.EV.5                              -22.72423   15.55392  -1.461 0.144582    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 219.5801)
## 
##     Null deviance: 204953  on 577  degrees of freedom
## Residual deviance: 121867  on 555  degrees of freedom
##   (1 observation deleted due to missingness)
## AIC: 4781.2
## 
## Number of Fisher Scoring iterations: 2
```

### Without age threshold

#### 8 out of 13 genes


```
## 
## Call:
## glm(formula = LOI_R ~ Age.of.Death + Institution + Gender + PMI..in.hours. + 
##     Dx + DLPFC_RNA_isolation..RIN + DLPFC_RNA_isolation..RIN.2 + 
##     DLPFC_RNA_report..Clustered.Library.Batch + Ancestry.EV.1 + 
##     Ancestry.EV.2 + Ancestry.EV.3 + Ancestry.EV.4 + Ancestry.EV.5, 
##     data = FULL)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -53.366  -11.375    1.408   11.331   49.625  
## 
## Coefficients:
##                                              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                                -113.40030   48.82905  -2.322 0.020573 *  
## Age.of.Death                                 -0.15857    0.05235  -3.029 0.002567 ** 
## InstitutionPenn                             -15.09527    2.24123  -6.735 4.10e-11 ***
## InstitutionPitt                              17.71864    2.41792   7.328 8.29e-13 ***
## GenderMale                                    1.74473    1.58324   1.102 0.270939    
## PMI..in.hours.                               -0.18658    0.07246  -2.575 0.010283 *  
## DxControl                                     2.43197    2.81880   0.863 0.388638    
## DxSCZ                                         2.09822    2.89410   0.725 0.468758    
## DLPFC_RNA_isolation..RIN                     48.59891   13.26109   3.665 0.000271 ***
## DLPFC_RNA_isolation..RIN.2                   -3.15831    0.89384  -3.533 0.000444 ***
## DLPFC_RNA_report..Clustered.Library.BatchA   -5.04292    4.02801  -1.252 0.211110    
## DLPFC_RNA_report..Clustered.Library.BatchB   -8.81720    3.75220  -2.350 0.019129 *  
## DLPFC_RNA_report..Clustered.Library.BatchC   -4.20825    4.14459  -1.015 0.310376    
## DLPFC_RNA_report..Clustered.Library.BatchD   -8.49083    3.78426  -2.244 0.025244 *  
## DLPFC_RNA_report..Clustered.Library.BatchE   -4.21118    3.78414  -1.113 0.266255    
## DLPFC_RNA_report..Clustered.Library.BatchF   -3.60987    3.83744  -0.941 0.347269    
## DLPFC_RNA_report..Clustered.Library.BatchG   -0.04775    5.16078  -0.009 0.992621    
## DLPFC_RNA_report..Clustered.Library.BatchH   -9.26188    5.82917  -1.589 0.112656    
## Ancestry.EV.1                               -12.66334   19.45112  -0.651 0.515294    
## Ancestry.EV.2                                12.77744   21.22616   0.602 0.547442    
## Ancestry.EV.3                               -14.80515   18.66940  -0.793 0.428107    
## Ancestry.EV.4                                34.11570   19.23850   1.773 0.076727 .  
## Ancestry.EV.5                                -2.53086   18.13910  -0.140 0.889086    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 298.5986)
## 
##     Null deviance: 273860  on 577  degrees of freedom
## Residual deviance: 165722  on 555  degrees of freedom
##   (1 observation deleted due to missingness)
## AIC: 4958.9
## 
## Number of Fisher Scoring iterations: 2
```

#### All 13 genes


```
## 
## Call:
## glm(formula = LOI_R ~ Age.of.Death + Institution + Gender + PMI..in.hours. + 
##     Dx + DLPFC_RNA_isolation..RIN + DLPFC_RNA_isolation..RIN.2 + 
##     DLPFC_RNA_report..Clustered.Library.Batch + Ancestry.EV.1 + 
##     Ancestry.EV.2 + Ancestry.EV.3 + Ancestry.EV.4 + Ancestry.EV.5, 
##     data = FULL)
## 
## Deviance Residuals: 
##     Min       1Q   Median       3Q      Max  
## -47.751   -9.055    0.919    9.552   43.508  
## 
## Coefficients:
##                                             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                                -90.13691   42.22684  -2.135 0.033232 *  
## Age.of.Death                                -0.08268    0.04527  -1.826 0.068356 .  
## InstitutionPenn                            -15.22339    1.93819  -7.854 2.09e-14 ***
## InstitutionPitt                             15.60847    2.09099   7.465 3.25e-13 ***
## GenderMale                                  -0.40588    1.36917  -0.296 0.767005    
## PMI..in.hours.                              -0.10115    0.06266  -1.614 0.107051    
## DxControl                                    2.18285    2.43767   0.895 0.370926    
## DxSCZ                                        2.67131    2.50279   1.067 0.286285    
## DLPFC_RNA_isolation..RIN                    41.69451   11.46805   3.636 0.000303 ***
## DLPFC_RNA_isolation..RIN.2                  -2.69139    0.77298  -3.482 0.000537 ***
## DLPFC_RNA_report..Clustered.Library.BatchA  -6.79726    3.48338  -1.951 0.051519 .  
## DLPFC_RNA_report..Clustered.Library.BatchB -10.47717    3.24486  -3.229 0.001316 ** 
## DLPFC_RNA_report..Clustered.Library.BatchC  -6.17352    3.58419  -1.722 0.085549 .  
## DLPFC_RNA_report..Clustered.Library.BatchD  -9.89373    3.27258  -3.023 0.002617 ** 
## DLPFC_RNA_report..Clustered.Library.BatchE  -6.16445    3.27249  -1.884 0.060125 .  
## DLPFC_RNA_report..Clustered.Library.BatchF  -5.88960    3.31858  -1.775 0.076489 .  
## DLPFC_RNA_report..Clustered.Library.BatchG  -1.82395    4.46298  -0.409 0.682929    
## DLPFC_RNA_report..Clustered.Library.BatchH -10.93461    5.04100  -2.169 0.030496 *  
## Ancestry.EV.1                              -15.10906   16.82112  -0.898 0.369458    
## Ancestry.EV.2                                2.95160   18.35615   0.161 0.872312    
## Ancestry.EV.3                               -3.81257   16.14510  -0.236 0.813408    
## Ancestry.EV.4                               36.74410   16.63725   2.209 0.027614 *  
## Ancestry.EV.5                              -22.85950   15.68650  -1.457 0.145607    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## (Dispersion parameter for gaussian family taken to be 223.31)
## 
##     Null deviance: 204953  on 577  degrees of freedom
## Residual deviance: 123937  on 555  degrees of freedom
##   (1 observation deleted due to missingness)
## AIC: 4791
## 
## Number of Fisher Scoring iterations: 2
```

[regression-slide]: https://docs.google.com/presentation/d/1YvpA1AJ-zzir1Iw0F25tO9x8gkSAzqaO4fjB7K3zBhE/edit#slide=id.p9
[CI-wiki]: https://en.wikipedia.org/wiki/Student%27s_t-distribution#Confidence_intervals
[weighted]: https://en.wikipedia.org/wiki/Linear_least_squares_(mathematics)#Weighted_linear_least_squares
