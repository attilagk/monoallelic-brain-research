

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

The strategy is to start from a model $M1$ that includes the strongest predictor terms based on previous analysis using separately modeling genes, then confirm the significant effect of these terms.  After this, single---mainly technical---terms are added separately (not yet accumulatively) to $M1$ to asses the effect of these.  The terms that improve the fit according to some criterion will be jointly added to $M1$ resulting in $M2$.  A similar procedure involving two sets of biological terms takes model selection further resulting in $M3$ and $M4$.

### Starting with model $M1$


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

* every term improves model fit by all three criteria (AIC, BIC, $\chi^2$ test)
* in particular, *Age* conditioned on *Gene* leads to the 3rd largest improvement (behind *Gene* and *Individual*)
    * however, this result needs to be confirmed in the context of a more general model 

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

* interaction of *Gene* with certain other predictors leads to large improvement in fit
    * this is not surprising since *Gene* has by far the largest impact among other main effects (see above)
    * but the fact that the improvement is due to accounting for interaction between a biological predictor (*Age*) and some technical predictor (*RIN*, *Insitution*, or *PMI*) is difficult to explain mechanistically
* *Age* shows essentially no interaction with *Institution* or *RNA_batch*
    * this is expected: the effect of *Age* should not depend on e.g. the *Institution*; so the result shows the advantage of joint modeling because previous results with separate, gene-based, modeling suggested strong interaction between *Age* and *Institution*

### Extending $M2$ with biological terms


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

* ancestry components are not equally important: only *Ancestry.1* and *Ancestry.3* appears to matter
* their main effect is weaker than their interaction with *Gene*
    * this is expected: the effect of *Age* should not depend on e.g. the *Institution*; so the result shows the advantage of joint modeling because previous results with separate, gene-based, modeling suggested strong interaction between *Age* and *Institution*

### Model $M3$ and $M4$

Based on the above findings the linear predictor for model $M3$ and $M4$ is defined as below, and $M3$ and $M4$ are fitted:


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

While the fitting of the more simple model $M3$ converges that of the more complex $M4$ does not.  Therefore $M3$ will be used for subsequent analysis.

### Adding more biological terms

Now add *Gender* and *Dx* either as fixed effect, a random effect (both as main effect and as an interaction term with *Gene*)


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

Note that adding *Gender* and *Dx* in the following way is not the proper way because they are factors.  Moreover, the fit does not converge.


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

Only *Gender*---but not *Dx*---improves the fit a little and only in interaction with *Gene*

### An important result


```r
av.4 <- anova(update(M1b, . ~ . + (scale(RIN) + scale(Ancestry.1) | Gene)),
              update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(Ancestry.1) | Gene)))
```

```
## refitting model(s) with ML (instead of REML)
```


```r
print.av(summarize.anova(list(av.4)))
```

```
##   Delta.AIC Delta.BIC    Chisq    df      p.Chi
## 1     -23.4       4.7     31.4     4    2.5e-06
```

### Saving model formulas

Save all model formulas in the `results` directory


```r
write.csv(data.frame(lapply(list(M1 = M1, M1b = M1b, M2 = M2, M3 = M3, M4 = M4),
                            function(x) as.character(formula(x)))),
          file = "../../results/M-formulas.csv", row.names = FALSE)
```

