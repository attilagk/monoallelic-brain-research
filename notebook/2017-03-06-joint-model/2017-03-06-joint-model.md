---
layout: default
tags: [regression, anova ]
---

An ad-hoc forward model selection strategy is used to evaluate the effect of the terms of the linear predictor in mixed effects models.  A strongly reduced starting model is prompted by earlier analysis using linear models conditioned on genes.  The model space is explored using that prior knowledge combined with information criteria (BIC and AIC) and likelihood ratio test.  The biologically interesting finding is that Ancestry.1, Ancestry.3 and Age have significant gene-specific effects but no gene-independent effects.  Dx (schizophrenia) has neither gene-specific nor gene-independent effect.


```
## Loading required package: Matrix
```

Selected genes (inferred to be imprinted)


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
```

Prepare data frame including all 30 selected genes


```r
dat <- merge.data(gene.ids = gene.ids)
```

## Model selection

The strategy is to start from a model $$M1$$ that includes the strongest predictor terms based on previous analysis using separately modeling genes, then confirm the significant effect of these terms.  After this, single---mainly technical---terms are added separately (not yet accumulatively) to $$M1$$ to asses the effect of these.  The terms that improve the fit according to some criterion will be jointly added to $$M1$$ resulting in $$M2$$.  A similar procedure involving two sets of biological terms takes model selection further resulting in $$M3$$, $$M4$$, etc.

### Starting with model $$M1$$


```r
fo <- Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) + (scale(Age) | Gene)
(M1 <- lmer(fo, data = dat))
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) +  
##     (scale(Age) | Gene)
##    Data: dat
## REML criterion at convergence: 20464.87
## Random effects:
##  Groups                 Name        Std.Dev. Corr 
##  Institution:Individual (Intercept) 0.3175        
##  Gene                   (Intercept) 1.2355        
##                         scale(Age)  0.1388   -0.51
##  RNA_batch              (Intercept) 0.1581        
##  Institution            (Intercept) 0.2547        
##  Residual                           0.7938        
## Number of obs: 8213, groups:  
## Institution:Individual, 579; Gene, 30; RNA_batch, 9; Institution, 3
## Fixed Effects:
## (Intercept)   scale(RIN)  
##     3.11696      0.06006
```


```r
av <- list()
av$RIN <- anova(update(M1, . ~ . - scale(RIN)), M1)
av$RNA_batch <- anova(update(M1, . ~ . - (1 | RNA_batch)), M1)
av$Institution <- anova(update(M1, . ~ . - (1 | Institution)), M1)
av$Individual <- anova(update(M1, . ~ . - (1 | Institution:Individual)), M1)
av$Age.Gene <- anova(update(M1, . ~ . - (scale(Age) | Gene) + (1 | Gene)), M1)
av$Gene <- anova(update(M1, . ~ . - (scale(Age) | Gene)), update(M1, . ~ . - (scale(Age) | Gene) + (1 | Gene)))
```


```r
summarize.anova <- function(av)
    do.call(rbind, lapply(av, function(x)
                          data.frame(Delta.AIC = diff(x$AIC),
                                     Delta.BIC = diff(x$BIC),
                                     Chisq = x$Chisq[2],
                                     df = x[["Chi Df"]][2],
                                     p.Chi = x[["Pr(>Chisq)"]][2])))
print.av(summarize.anova(av))
```

```
##             Delta.AIC Delta.BIC    Chisq    df      p.Chi
## RIN              -8.0      -1.0     10.0     1    1.5e-03
## RNA_batch       -25.7     -18.7     27.7     1    1.4e-07
## Institution     -61.2     -54.2     63.2     1    1.9e-15
## Individual     -453.8    -446.8    455.8     1   4.0e-101
## Age.Gene       -192.2    -178.2    196.2     2    2.5e-43
## Gene          -8700.3   -8693.3   8702.3     1        0.0
```

#### Notable results

* every term improves model fit by all three criteria (AIC, BIC, $$\chi^2$$ test)
* in particular, the *(Age\|Gene)* leads to the 3rd largest improvement (behind *(1\|Gene)* and *(1\|Individual)*)
    * however, this result needs to be confirmed in the context of a more general model 

### Extending $$M1$$ to $$M2$$ with technical terms


```r
av.1 <- list()
# interactions between random-effect factors
av.1$Gene.Instit <- anova(M1, update(M1, . ~ . + (1 | Gene:Institution)))
av.1$Gene.RNA_batch <- anova(M1, update(M1, . ~ . + (1 | Gene:RNA_batch)))
av.1$RNA_batch.Instit <- anova(M1, update(M1, . ~ . + (1 | RNA_batch:Institution)))
# interactions between random-effect covariates and factors
# RIN
av.1$RIN.Instit <- anova(M1, update(M1, . ~ . - (1 | Institution) + (scale(RIN) | Institution)))
av.1$RIN.RNA_batch <- anova(M1, update(M1, . ~ . - (1 | RNA_batch) + (scale(RIN) | RNA_batch)))
av.1$RIN.Gene <- anova(M1, update(M1, . ~ . - (1 | Gene) + (scale(RIN) | Gene)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, : Model is nearly unidentifiable: large eigenvalue ratio
##  - Rescale variables?
```

```r
# PMI
av.1$PMI.Instit <- anova(M1, update(M1, . ~ . - (1 | Institution) + (scale(PMI) | Institution)))
av.1$PMI.RNA_batch <- anova(M1, update(M1, . ~ . - (1 | RNA_batch) + (scale(PMI) | RNA_batch)))
av.1$PMI.Gene <- anova(M1, update(M1, . ~ . - (1 | Gene) + (scale(PMI) | Gene)))
av.1$PMI <- anova(M1, update(M1, . ~ . + scale(PMI)))
# Age revisited
av.1$Age.Instit <- anova(M1, update(M1, . ~ . - (1 | Institution) + (scale(Age) | Institution)))
av.1$Age.RNA_batch <- anova(M1, update(M1, . ~ . - (1 | RNA_batch) + (scale(Age) | RNA_batch)))
av.1$Age <- anova(M1, update(M1, . ~ . + scale(Age)))
```


```r
print.av(summarize.anova(av.1))
```

```
##                  Delta.AIC Delta.BIC    Chisq    df      p.Chi
## Gene.Instit         -136.5    -129.5    138.5     1    5.8e-32
## Gene.RNA_batch         0.3       7.3      1.7     1    1.9e-01
## RNA_batch.Instit       2.0       9.0      0.0     1        1.0
## RIN.Instit             2.1      16.1      1.9     2    3.8e-01
## RIN.RNA_batch          3.4      17.4      0.6     2    7.3e-01
## RIN.Gene            -326.0    -305.0    332.0     3    1.2e-71
## PMI.Instit             3.7      17.7      0.3     2    8.5e-01
## PMI.RNA_batch          2.0      16.0      2.0     2    3.6e-01
## PMI.Gene             -11.4       9.7     17.4     3    6.0e-04
## PMI                    1.9       8.9      0.1     1    7.6e-01
## Age.Instit             1.2      15.2      2.8     2    2.4e-01
## Age.RNA_batch          2.3      16.4      1.7     2    4.4e-01
## Age                    1.8       8.8      0.2     1    6.5e-01
```

#### Notable results

* the *Gene*-specific random slope effect of certain predictors---e.g. *(RIN\|Gene)*---leads to large improvements in fit
    * this is not surprising since the improvement by the *(1\|Gene)* is the largest among all tested terms (see above)
    * but it is difficult to explain mechanistically the improvement by e.g. *(RIN\|Gene)* because how could *Gene* influence the effect of *RIN*?

### Extending $$M2$$ with biological terms


```r
M1b <- update(M1, . ~ . + (1 | Gene:Institution) - (scale(Age) | Gene))
(M2 <- update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(PMI) | Gene)))
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) +  
##     (1 | Gene:Institution) + (scale(Age) + scale(RIN) + scale(PMI) |  
##     Gene)
##    Data: dat
## REML criterion at convergence: 20019.65
## Random effects:
##  Groups                 Name        Std.Dev. Corr             
##  Institution:Individual (Intercept) 0.32736                   
##  Gene:Institution       (Intercept) 0.19668                   
##  Gene                   (Intercept) 1.22177                   
##                         scale(Age)  0.07065  -0.35            
##                         scale(RIN)  0.17682   0.21 -0.71      
##                         scale(PMI)  0.01870   0.24  0.82 -0.66
##  RNA_batch              (Intercept) 0.16600                   
##  Institution            (Intercept) 0.21699                   
##  Residual                           0.76397                   
## Number of obs: 8213, groups:  
## Institution:Individual, 579; Gene:Institution, 90; Gene, 30; RNA_batch, 9; Institution, 3
## Fixed Effects:
## (Intercept)   scale(RIN)  
##     3.10931      0.05196
```


```r
av.2 <- list()
# Ancestry
av.2$Ances1 <- anova(M2, update(M2, . ~ . + scale(Ancestry.1)))
av.2$Ances2 <- anova(M2, update(M2, . ~ . + scale(Ancestry.2)))
av.2$Ances3 <- anova(M2, update(M2, . ~ . + scale(Ancestry.3)))
av.2$Ances4 <- anova(M2, update(M2, . ~ . + scale(Ancestry.4)))
av.2$Ances5 <- anova(M2, update(M2, . ~ . + scale(Ancestry.5)))
av.2$Ances1.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(PMI) + scale(Ancestry.1) | Gene)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
## $checkConv, : unable to evaluate scaled gradient
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
## $checkConv, : Model failed to converge: degenerate Hessian with 1 negative
## eigenvalues
```

```r
av.2$Ances2.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(PMI) + scale(Ancestry.2) | Gene)))
av.2$Ances3.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(PMI) + scale(Ancestry.3) | Gene)))
av.2$Ances4.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(PMI) + scale(Ancestry.4) | Gene)))
av.2$Ances5.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(PMI) + scale(Ancestry.5) | Gene)))
```


```r
print.av(summarize.anova(av.2))
```

```
##             Delta.AIC Delta.BIC    Chisq    df      p.Chi
## Ances1           -0.2       6.8      2.2     1    1.3e-01
## Ances2            1.7       8.7      0.3     1    5.6e-01
## Ances3            1.9       8.9      0.1     1    7.7e-01
## Ances4            0.9       7.9      1.1     1    2.9e-01
## Ances5            2.0       9.0      0.0     1    8.4e-01
## Ances1.Gene     -56.8     -21.7     66.8     5    4.8e-13
## Ances2.Gene       5.1      40.2      4.9     5    4.3e-01
## Ances3.Gene      -3.7      31.3     13.7     5    1.7e-02
## Ances4.Gene       7.9      42.9      2.1     5    8.3e-01
## Ances5.Gene       9.1      44.1      0.9     5    9.7e-01
```

#### Notable results

* *Ancestry.n* components are not equally important: only *Ancestry.1* and *Ancestry.3* appears to matter
* their effects are *Gene*-specific

### Model $$M3$$ and $$M4$$

Based on the above findings the linear predictor for model $$M3$$ and $$M4$$ is defined as below, and $$M3$$ and $$M4$$ are fitted:


```r
(M3 <- update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene)))
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) +  
##     (1 | Gene:Institution) + (scale(Age) + scale(RIN) + scale(Ancestry.1) |  
##     Gene)
##    Data: dat
## REML criterion at convergence: 19956.89
## Random effects:
##  Groups                 Name              Std.Dev. Corr             
##  Institution:Individual (Intercept)       0.32862                   
##  Gene:Institution       (Intercept)       0.19890                   
##  Gene                   (Intercept)       1.22187                   
##                         scale(Age)        0.06926  -0.44            
##                         scale(RIN)        0.17915   0.20 -0.75      
##                         scale(Ancestry.1) 0.08387   0.24  0.04  0.10
##  RNA_batch              (Intercept)       0.16599                   
##  Institution            (Intercept)       0.21712                   
##  Residual                                 0.75881                   
## Number of obs: 8213, groups:  
## Institution:Individual, 579; Gene:Institution, 90; Gene, 30; RNA_batch, 9; Institution, 3
## Fixed Effects:
## (Intercept)   scale(RIN)  
##     3.14180      0.05632
```

```r
M4 <- update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(PMI) +
                           scale(Ancestry.1) + scale(Ancestry.3) | Gene) +
             scale(Ancestry.1))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
## $checkConv, : unable to evaluate scaled gradient
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
## $checkConv, : Model failed to converge: degenerate Hessian with 1 negative
## eigenvalues
```

While the fitting of the more simple model $$M3$$ converges that of the more complex $$M4$$ does not.  Therefore $$M3$$ will be used for subsequent analysis.

### Adding more biological terms

Now add terms including gender *Gender* and *Dx*:


```r
av.3 <- list()
# Gender
av.3$Gender.fix <- anova(M3, update(M3, . ~ . + Gender))
av.3$Gender.ran <- anova(M3, m <- update(M3, . ~ . + (1 | Gender)))
av.3$Gender.Gene <- anova(m, update(m, . ~ . + (1 | Gender:Gene)))
# Dx
av.3$Dx.fix <- anova(M3, update(M3, . ~ . + Dx))
av.3$Dx.ran <- anova(M3, m <- update(M3, . ~ . + (1 | Dx)))
av.3$Dx.Gene <- anova(m, update(m, . ~ . + (1 | Dx:Gene)))
```

Note that adding *Gender* and *Dx* in the following way (i.e.*(Gender+...\|Gene)*) is not the proper way because *Gender* and *Dx* are factors, whose effect cannot be modeled as a slope.  Despite this the `lme4` allows fitting models with such formulas.  But as the warnings below show the fit does not converge.


```r
# factor playing the role of slope doesn't make much sense
invisible(update(M1b, . ~ . + (Gender + scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene)))
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
## $checkConv, : unable to evaluate scaled gradient
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
## $checkConv, : Model failed to converge: degenerate Hessian with 1 negative
## eigenvalues
```

#### Results


```r
print.av(summarize.anova(av.3))
```

```
##             Delta.AIC Delta.BIC    Chisq    df      p.Chi
## Gender.fix        1.0       8.0      1.0     1    3.2e-01
## Gender.ran        2.0       9.0      0.0     1        1.0
## Gender.Gene      -5.8       1.2      7.8     1    5.2e-03
## Dx.fix            3.9      17.9      0.1     2    9.5e-01
## Dx.ran            2.0       9.0      0.0     1        1.0
## Dx.Gene           0.5       7.5      1.5     1    2.1e-01
```

Only *(1\|Gender:Gene)* improves fit appreciably, while *Dx* has no effect (neither random, nor fixed, nor in interaction with *Gene*).

### The impact of removing *(Age|Gene)* from M3

This is to investigate how much improvement *(Age\|Gene)* can lead to in the more complex *M3* (cf. improvement by *(Age\|Gene)* in *M1*).  The result below shows the improvement is less in *M3* but still highly significant.


```r
print.av(summarize.anova(list(anova(update(M1b, . ~ . + (scale(RIN) + scale(Ancestry.1) | Gene)),
              update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene))))))
```

```
## refitting model(s) with ML (instead of REML)
```

```
##   Delta.AIC Delta.BIC    Chisq    df      p.Chi
## 1     -23.4       4.7     31.4     4    2.5e-06
```

### *M5* and *M6*: more biological terms

Based on the findings above *M5* and *M5.Dx.Gene* are defined (the latter is renamed to *M6* later).  These are thus better fitting models than any of the previous ones.


```r
(M5 <- update(M1b, . ~ . + (1 | Gender:Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)))
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) +  
##     (1 | Gene:Institution) + (1 | Gender:Gene) + (scale(Age) +  
##     scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)
##    Data: dat
## REML criterion at convergence: 19921.19
## Random effects:
##  Groups                 Name              Std.Dev. Corr                   
##  Institution:Individual (Intercept)       0.32947                         
##  Gene:Institution       (Intercept)       0.19838                         
##  Gender:Gene            (Intercept)       0.05939                         
##  Gene                   (Intercept)       1.22077                         
##                         scale(Age)        0.06600  -0.45                  
##                         scale(RIN)        0.17777   0.19 -0.71            
##                         scale(Ancestry.1) 0.08317   0.24  0.04  0.04      
##                         scale(Ancestry.3) 0.04528  -0.31  0.37 -0.20 -0.88
##  RNA_batch              (Intercept)       0.16581                         
##  Institution            (Intercept)       0.21710                         
##  Residual                                 0.75584                         
## Number of obs: 8213, groups:  
## Institution:Individual, 579; Gene:Institution, 90; Gender:Gene, 60; Gene, 30; RNA_batch, 9; Institution, 3
## Fixed Effects:
## (Intercept)   scale(RIN)  
##     3.13115      0.05611
```

```r
(M5.Dx.Gene <- update(M5, . ~ . + (1 | Dx:Gene)))
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) +  
##     (1 | Gene:Institution) + (1 | Gender:Gene) + (scale(Age) +  
##     scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) +  
##     (1 | Dx:Gene)
##    Data: dat
## REML criterion at convergence: 19919.61
## Random effects:
##  Groups                 Name              Std.Dev. Corr                   
##  Institution:Individual (Intercept)       0.32953                         
##  Dx:Gene                (Intercept)       0.03872                         
##  Gene:Institution       (Intercept)       0.19903                         
##  Gender:Gene            (Intercept)       0.06114                         
##  Gene                   (Intercept)       1.22093                         
##                         scale(Age)        0.06523  -0.45                  
##                         scale(RIN)        0.17692   0.20 -0.70            
##                         scale(Ancestry.1) 0.08326   0.24  0.05  0.04      
##                         scale(Ancestry.3) 0.04468  -0.30  0.37 -0.20 -0.88
##  RNA_batch              (Intercept)       0.16608                         
##  Institution            (Intercept)       0.21711                         
##  Residual                                 0.75525                         
## Number of obs: 8213, groups:  
## Institution:Individual, 579; Dx:Gene, 90; Gene:Institution, 90; Gender:Gene, 60; Gene, 30; RNA_batch, 9; Institution, 3
## Fixed Effects:
## (Intercept)   scale(RIN)  
##     3.13076      0.05625
```

The following sequence of ANOVAs tests the effects of biological terms in these better fitting but also more complex models:


```r
av.5 <- list()
av.5$Age.Gene <-
    anova(update(M5, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)),
          M5)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.1.Gene <-
    anova(update(M5, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.3) | Gene)),
          M5)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.2.Gene <-
    anova(M5,
          update(M5, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.2) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.3.Gene <-
    anova(update(M5, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene)),
          M5)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.4.Gene <-
    anova(M5,
          update(M5, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.4) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.5.Gene <-
    anova(M5,
          update(M5, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.5) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Gender.Gene <-
    anova(update(M5, . ~ . - (1 | Gender:Gene)), M5)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Gender <- anova(M5, update(M5, . ~ . + (1 | Gender)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
# note that the next one is an addition to M5
av.5$Dx.Gene <-
    anova(M5, update(M5, . ~ . + (1 | Dx:Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Dx <- anova(M5, update(M5, . ~ . + (1 | Dx)))
```

```
## refitting model(s) with ML (instead of REML)
```

Now fit *M6* and perform ANOVAs similarly to *M5*


```r
M6 <- M5.Dx.Gene
av.6 <- list()
av.6$Age.Gene <-
    anova(update(M6, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)),
          M6)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.1.Gene <-
    anova(update(M6, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.3) | Gene)),
          M6)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.2.Gene <-
    anova(M6,
          update(M6, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.2) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.3.Gene <-
    anova(update(M6, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene)),
          M6)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.4.Gene <-
    anova(M6,
          update(M6, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.4) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.5.Gene <-
    anova(M6,
          update(M6, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.5) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Gender.Gene <-
    anova(update(M6, . ~ . - (1 | Gender:Gene)), M6)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Gender <- anova(M6, update(M6, . ~ . + (1 | Gender)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Dx.Gene <-
    anova(update(M6, . ~ . - (1 | Dx:Gene)), M6)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Dx <- anova(M6, update(M6, . ~ . + (1 | Dx)))
```

```
## refitting model(s) with ML (instead of REML)
```


```r
av.5$Gene <-
    anova(update(M5, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (-1 + scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)),
          M5)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Gene <-
    anova(update(M6, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (-1 + scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)),
          M6)
```

```
## refitting model(s) with ML (instead of REML)
```

Effect of *Age*, *Ancestry.n* independently from other variables.


```r
av.5$Age <- anova(M5, update(M5, . ~ . + scale(Age)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.1 <- anova(M5, update(M5, . ~ . + scale(Ancestry.1)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.2 <- anova(M5, update(M5, . ~ . + scale(Ancestry.2)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.3 <- anova(M5, update(M5, . ~ . + scale(Ancestry.3)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.4 <- anova(M5, update(M5, . ~ . + scale(Ancestry.4)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5$Ancestry.5 <- anova(M5, update(M5, . ~ . + scale(Ancestry.5)))
```

```
## refitting model(s) with ML (instead of REML)
```


```r
av.6$Age <- anova(M6, update(M6, . ~ . + scale(Age)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.1 <- anova(M6, update(M6, . ~ . + scale(Ancestry.1)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.2 <- anova(M6, update(M6, . ~ . + scale(Ancestry.2)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.3 <- anova(M6, update(M6, . ~ . + scale(Ancestry.3)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.4 <- anova(M6, update(M6, . ~ . + scale(Ancestry.4)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.6$Ancestry.5 <- anova(M6, update(M6, . ~ . + scale(Ancestry.5)))
```

```
## refitting model(s) with ML (instead of REML)
```

#### Results of ANOVA


```r
write.csv(summarize.anova(av.5), "../../results/anova-mixed-M5.csv")
print.av(summarize.anova(av.5))
```

```
##                 Delta.AIC Delta.BIC    Chisq    df      p.Chi
## Age.Gene            -18.9      16.2     28.9     5    2.5e-05
## Ancestry.1.Gene     -71.2     -36.2     81.2     5    4.6e-16
## Ancestry.2.Gene       6.3      48.3      5.7     6    4.5e-01
## Ancestry.3.Gene     -17.9      17.2     27.9     5    3.8e-05
## Ancestry.4.Gene       6.9      48.9      5.1     6    5.3e-01
## Ancestry.5.Gene      10.8      52.9      1.2     6    9.8e-01
## Gender.Gene          -5.7       1.3      7.7     1    5.5e-03
## Gender                2.0       9.0      0.0     1        1.0
## Dx.Gene               0.4       7.4      1.6     1    2.1e-01
## Dx                    2.0       9.0      0.0     1        1.0
## Gene               -126.8     -91.7    136.8     5    8.5e-28
## Age                   1.3       8.3      0.7     1    3.9e-01
## Ancestry.1            0.6       7.6      1.4     1    2.4e-01
## Ancestry.2            1.7       8.7      0.3     1    5.6e-01
## Ancestry.3            1.6       8.6      0.4     1    5.4e-01
## Ancestry.4            0.9       8.0      1.1     1    3.0e-01
## Ancestry.5            1.9       9.0      0.1     1    8.1e-01
```

```r
write.csv(summarize.anova(av.6), "../../results/anova-mixed-M6.csv")
print.av(summarize.anova(av.6))
```

```
##                 Delta.AIC Delta.BIC    Chisq    df      p.Chi
## Age.Gene            -17.5      17.5     27.5     5    4.5e-05
## Ancestry.1.Gene     -71.4     -36.4     81.4     5    4.2e-16
## Ancestry.2.Gene       6.4      48.4      5.6     6    4.7e-01
## Ancestry.3.Gene     -17.1      18.0     27.1     5    5.5e-05
## Ancestry.4.Gene       6.9      49.0      5.1     6    5.4e-01
## Ancestry.5.Gene      10.9      53.0      1.1     6    9.8e-01
## Gender.Gene          -6.4       0.6      8.4     1    3.8e-03
## Gender                2.0       9.0      0.0     1        1.0
## Dx.Gene               0.4       7.4      1.6     1    2.1e-01
## Dx                    2.0       9.0      0.0     1        1.0
## Gene               -126.1     -91.0    136.1     5    1.2e-27
## Age                   1.2       8.3      0.8     1    3.9e-01
## Ancestry.1            0.6       7.6      1.4     1    2.4e-01
## Ancestry.2            1.7       8.7      0.3     1    5.6e-01
## Ancestry.3            1.6       8.6      0.4     1    5.4e-01
## Ancestry.4            1.0       8.0      1.0     1    3.1e-01
## Ancestry.5            1.9       9.0      0.1     1    8.1e-01
```

#### Revision for Nature Communications


```r
M5.odd <- lmer(formula(M5), data = dat, subset = Gene %in% gene.ids[seq(from = 1, to = length(gene.ids), by = 2)])
M5.even <- lmer(formula(M5), data = dat, subset = Gene %in% gene.ids[seq(from = 2, to = length(gene.ids), by = 2)])
```


```r
av.5.odd <- list()
av.5.odd$Age.Gene <-
    anova(update(M5.odd, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)),
          M5.odd)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.1.Gene <-
    anova(update(M5.odd, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.3) | Gene)),
          M5.odd)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.2.Gene <-
    anova(M5.odd,
          update(M5.odd, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.2) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.3.Gene <-
    anova(update(M5.odd, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene)),
          M5.odd)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.4.Gene <-
    anova(M5.odd,
          update(M5.odd, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.4) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.5.Gene <-
    anova(M5.odd,
          update(M5.odd, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.5) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Gender.Gene <-
    anova(update(M5.odd, . ~ . - (1 | Gender:Gene)), M5.odd)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Gender <- anova(M5.odd, update(M5.odd, . ~ . + (1 | Gender)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
# note that the next one is an addition to M5.odd
av.5.odd$Dx.Gene <-
    anova(M5.odd, update(M5.odd, . ~ . + (1 | Dx:Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Dx <- anova(M5.odd, update(M5.odd, . ~ . + (1 | Dx)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Gene <-
    anova(update(M5.odd, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (-1 + scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)),
          M5.odd)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Age <- anova(M5.odd, update(M5.odd, . ~ . + scale(Age)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.1 <- anova(M5.odd, update(M5.odd, . ~ . + scale(Ancestry.1)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.2 <- anova(M5.odd, update(M5.odd, . ~ . + scale(Ancestry.2)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.3 <- anova(M5.odd, update(M5.odd, . ~ . + scale(Ancestry.3)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.4 <- anova(M5.odd, update(M5.odd, . ~ . + scale(Ancestry.4)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.odd$Ancestry.5 <- anova(M5.odd, update(M5.odd, . ~ . + scale(Ancestry.5)))
```

```
## refitting model(s) with ML (instead of REML)
```


```r
av.5.even <- list()
av.5.even$Age.Gene <-
    anova(update(M5.even, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)),
          M5.even)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.1.Gene <-
    anova(update(M5.even, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.3) | Gene)),
          M5.even)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.2.Gene <-
    anova(M5.even,
          update(M5.even, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.2) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.3.Gene <-
    anova(update(M5.even, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene)),
          M5.even)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.4.Gene <-
    anova(M5.even,
          update(M5.even, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.4) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.5.Gene <-
    anova(M5.even,
          update(M5.even, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.5) + scale(Ancestry.3) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Gender.Gene <-
    anova(update(M5.even, . ~ . - (1 | Gender:Gene)), M5.even)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Gender <- anova(M5.even, update(M5.even, . ~ . + (1 | Gender)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
# note that the next one is an addition to M5.even
av.5.even$Dx.Gene <-
    anova(M5.even, update(M5.even, . ~ . + (1 | Dx:Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Dx <- anova(M5.even, update(M5.even, . ~ . + (1 | Dx)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Gene <-
    anova(update(M5.even, . ~ . - (scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene) + (-1 + scale(Age) + scale(RIN) + scale(Ancestry.1) + scale(Ancestry.3) | Gene)),
          M5.even)
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Age <- anova(M5.even, update(M5.even, . ~ . + scale(Age)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.1 <- anova(M5.even, update(M5.even, . ~ . + scale(Ancestry.1)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.2 <- anova(M5.even, update(M5.even, . ~ . + scale(Ancestry.2)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.3 <- anova(M5.even, update(M5.even, . ~ . + scale(Ancestry.3)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.4 <- anova(M5.even, update(M5.even, . ~ . + scale(Ancestry.4)))
```

```
## refitting model(s) with ML (instead of REML)
```

```r
av.5.even$Ancestry.5 <- anova(M5.even, update(M5.even, . ~ . + scale(Ancestry.5)))
```

```
## refitting model(s) with ML (instead of REML)
```


```r
write.csv(summarize.anova(av.5.odd), "../../results/anova-mixed-M5.odd.csv")
print.av(summarize.anova(av.5.odd))
```

```
##                 Delta.AIC Delta.BIC    Chisq    df      p.Chi
## Age.Gene            -11.8      20.1     21.8     5    5.8e-04
## Ancestry.1.Gene     -40.1      -8.2     50.1     5    1.3e-09
## Ancestry.2.Gene       5.6      43.8      6.4     6    3.8e-01
## Ancestry.3.Gene     -13.3      18.5     23.3     5    2.9e-04
## Ancestry.4.Gene      11.6      49.9      0.4     6        1.0
## Ancestry.5.Gene       9.0      47.3      3.0     6    8.1e-01
## Gender.Gene          -2.2       4.2      4.2     1    4.0e-02
## Gender                2.0       8.4      0.0     1        1.0
## Dx.Gene               1.9       8.2      0.1     1    7.1e-01
## Dx                    2.0       8.4      0.0     1        1.0
## Gene                -61.2     -29.3     71.2     5    5.7e-14
## Age                  -0.0       6.4      2.0     1    1.6e-01
## Ancestry.1           -0.4       6.0      2.4     1    1.2e-01
## Ancestry.2            1.4       7.8      0.6     1    4.4e-01
## Ancestry.3            1.7       8.1      0.3     1    5.9e-01
## Ancestry.4           -1.0       5.4      3.0     1    8.2e-02
## Ancestry.5            1.3       7.6      0.7     1    3.9e-01
```

```r
write.csv(summarize.anova(av.5.even), "../../results/anova-mixed-M5.even.csv")
print.av(summarize.anova(av.5.even))
```

```
##                 Delta.AIC Delta.BIC    Chisq    df      p.Chi
## Age.Gene              5.1      36.4      4.9     5    4.3e-01
## Ancestry.1.Gene     -18.5      12.8     28.5     5    2.9e-05
## Ancestry.2.Gene       8.6      46.2      3.4     6    7.6e-01
## Ancestry.3.Gene       6.0      37.3      4.0     5    5.5e-01
## Ancestry.4.Gene       5.5      43.1      6.5     6    3.7e-01
## Ancestry.5.Gene       8.1      45.7      3.9     6    6.9e-01
## Gender.Gene           0.1       6.4      1.9     1    1.7e-01
## Gender                0.7       6.9      1.3     1    2.5e-01
## Dx.Gene              -0.0       6.2      2.0     1    1.6e-01
## Dx                    2.0       8.3      0.0     1        1.0
## Gene                -59.2     -27.9     69.2     5    1.5e-13
## Age                   2.0       8.2      0.0     1    8.6e-01
## Ancestry.1            1.8       8.1      0.2     1    6.6e-01
## Ancestry.2            1.9       8.1      0.1     1    7.1e-01
## Ancestry.3            1.6       7.9      0.4     1    5.4e-01
## Ancestry.4            2.0       8.3      0.0     1    9.5e-01
## Ancestry.5            1.7       8.0      0.3     1    6.0e-01
```

#### Summary

* significant *Gene*-specific random effect of *Age*, *Ancestry.1* *Ancestry.3*,  and *Gender*
* no effect of *Dx*


### Saving model formulas

Save all model formulas in the `results` directory


```r
write.csv(data.frame(lapply(list(M1 = M1, M1b = M1b, M2 = M2, M3 = M3, M4 = M4, M5 = M5, M6 = M6),
                            function(x) as.character(formula(x)))),
          file = "../../results/M-formulas.csv", row.names = FALSE)
```
