

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

The strategy is to start from a model $M1$ that includes the strongest predictor terms based on previous analysis using separately modeling genes, then confirm the significant effect of these terms.  After this, single---mainly technical---terms are added separately (not yet accumulatively) to $M1$ to asses the effect of these.  The terms that improve the fit according to some criterion will be jointly added to $M1$ resulting in $M2$

### Starting model $M1$


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
av$RIN <- anova(update(M1, . ~ . - RIN), M1)
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
summarize.anova(av)
```

```
##               Delta.AIC   Delta.BIC      Chisq df         p.Chi
## RIN             0.00000     0.00000    0.00000  0  1.000000e+00
## RNA_batch     -25.69463   -18.68116   27.69463  1  1.420564e-07
## Institution   -61.18655   -54.17308   63.18655  1  1.880277e-15
## Individual   -453.77916  -446.76569  455.77916  1 3.984989e-101
## Age.Gene     -192.22905  -178.20210  196.22905  2  2.451338e-43
## Gene        -8700.32456 -8693.31108 8702.32456  1  0.000000e+00
```

### Extending $M1$ to $M2$ with technical terms


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
summarize.anova(av.1)
```

```
##                     Delta.AIC   Delta.BIC        Chisq df        p.Chi
## Gene.Instit      -136.4699800 -129.456506 138.46997998  1 5.751836e-32
## Gene.RNA_batch      0.2955551    7.309029   1.70444492  1 1.917077e-01
## RNA_batch.Instit    2.0000000    9.013474   0.00000000  1 1.000000e+00
## RIN.Instit          2.0907634   16.117711   1.90923656  2 3.849591e-01
## RIN.RNA_batch       3.3818324   17.408780   0.61816755  2 7.341193e-01
## RIN.Gene         -326.0009803 -304.960560 332.00098032  3 1.176834e-71
## PMI.Instit          3.6657565   17.692704   0.33424348  2 8.460966e-01
## PMI.RNA_batch       1.9699532   15.996900   2.03004677  2 3.623940e-01
## PMI.Gene          -11.3562482    9.684172  17.35624821  3 5.969749e-04
## PMI                 1.9095940    8.923068   0.09040597  1 7.636617e-01
## Age.Instit          1.1806642   15.207611   2.81933578  2 2.442244e-01
## Age.RNA_batch       2.3483136   16.375261   1.65168640  2 4.378656e-01
## Age                 1.7975291    8.811003   0.20247091  1 6.527338e-01
```

### Extending $M2$ with biological terms


```r
M1b <- update(M1, . ~ . + (1 | Gene:Institution) - (scale(Age) | Gene))
(M2 <- update(M1b, . ~ . + (scale(Age) + scale(RIN) | Gene)))
```

```
## Linear mixed model fit by REML ['lmerMod']
## Formula: 
## Q ~ scale(RIN) + (1 | RNA_batch) + (1 | Institution) + (1 | Institution:Individual) +  
##     (1 | Gene:Institution) + (scale(Age) + scale(RIN) | Gene)
##    Data: dat
## REML criterion at convergence: 20023.88
## Random effects:
##  Groups                 Name        Std.Dev. Corr       
##  Institution:Individual (Intercept) 0.32727             
##  Gene:Institution       (Intercept) 0.19865             
##  Gene                   (Intercept) 1.22153             
##                         scale(Age)  0.06724  -0.39      
##                         scale(RIN)  0.17946   0.20 -0.71
##  RNA_batch              (Intercept) 0.16579             
##  Institution            (Intercept) 0.21666             
##  Residual                           0.76422             
## Number of obs: 8213, groups:  
## Institution:Individual, 579; Gene:Institution, 90; Gene, 30; RNA_batch, 9; Institution, 3
## Fixed Effects:
## (Intercept)   scale(RIN)  
##     3.10284      0.05333
```


```r
av.2 <- list()
# Ancestry
av.2$Ances1 <- anova(M2, update(M2, . ~ . + scale(Ancestry.1)))
av.2$Ances2 <- anova(M2, update(M2, . ~ . + scale(Ancestry.2)))
av.2$Ances3 <- anova(M2, update(M2, . ~ . + scale(Ancestry.3)))
av.2$Ances4 <- anova(M2, update(M2, . ~ . + scale(Ancestry.4)))
av.2$Ances5 <- anova(M2, update(M2, . ~ . + scale(Ancestry.5)))
av.2$Ances1.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene)))
av.2$Ances2.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(Ancestry.2) | Gene)))
av.2$Ances3.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(Ancestry.3) | Gene)))
av.2$Ances4.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(Ancestry.4) | Gene)))
av.2$Ances5.Gene <- anova(M2, update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(Ancestry.5) | Gene)))
```


```r
summarize.anova(av.2)
```

```
##               Delta.AIC  Delta.BIC       Chisq df        p.Chi
## Ances1       -0.2178647   6.795609  2.21786474  1 1.364216e-01
## Ances2        1.6489647   8.662438  0.35103533  1 5.535276e-01
## Ances3        1.9161077   8.929581  0.08389229  1 7.720904e-01
## Ances4        0.8485452   7.862019  1.15145483  1 2.832448e-01
## Ances5        1.9533949   8.966868  0.04660514  1 8.290795e-01
## Ances1.Gene -59.0249725 -30.971078 67.02497254  4 9.631377e-14
## Ances2.Gene   3.9732244  32.027119  4.02677564  4 4.023943e-01
## Ances3.Gene  -5.6506303  22.403264 13.65063034  4 8.497773e-03
## Ances4.Gene   5.6031803  33.657074  2.39681974  4 6.632021e-01
## Ances5.Gene   7.1409406  35.194835  0.85905936  4 9.303595e-01
```


```r
av.3 <- list()
# Gender
av.3$Gender <- anova(M2, update(M2, . ~ . + Gender))
av.3$Gender.Gene <- anova(M2, update(M2, . ~ . + (Gender + scale(Age) + scale(RIN) + scale(Ancestry.5) | Gene)))
# Dx
av.3$Dx <- anova(M2, update(M2, . ~ . + Dx))
av.3$Dx.Gene <- anova(M2, update(M1b, . ~ . + (Dx + scale(Age) + scale(RIN) + scale(Ancestry.5) | Gene)))
```


```r
summarize.anova(av.3)
```

```
## Error in lapply(av, function(x) data.frame(Delta.AIC = diff(x$AIC), Delta.BIC = diff(x$BIC), : object 'av.3' not found
```
