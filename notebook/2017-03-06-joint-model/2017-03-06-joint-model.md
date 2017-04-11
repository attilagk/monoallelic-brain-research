
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
signif(summarize.anova(av), digits = 3)
```

```
##             Delta.AIC Delta.BIC  Chisq df     p.Chi
## RIN             -8.04     -1.03   10.0  1  1.53e-03
## RNA_batch      -25.70    -18.70   27.7  1  1.42e-07
## Institution    -61.20    -54.20   63.2  1  1.88e-15
## Individual    -454.00   -447.00  456.0  1 3.98e-101
## Age.Gene      -192.00   -178.00  196.0  2  2.45e-43
## Gene         -8700.00  -8690.00 8700.0  1  0.00e+00
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
signif(summarize.anova(av.1), digits = 3)
```

```
##                  Delta.AIC Delta.BIC    Chisq df    p.Chi
## Gene.Instit       -136.000   -129.00 138.0000  1 5.75e-32
## Gene.RNA_batch       0.296      7.31   1.7000  1 1.92e-01
## RNA_batch.Instit     2.000      9.01   0.0000  1 1.00e+00
## RIN.Instit           2.090     16.10   1.9100  2 3.85e-01
## RIN.RNA_batch        3.380     17.40   0.6180  2 7.34e-01
## RIN.Gene          -326.000   -305.00 332.0000  3 1.18e-71
## PMI.Instit           3.670     17.70   0.3340  2 8.46e-01
## PMI.RNA_batch        1.970     16.00   2.0300  2 3.62e-01
## PMI.Gene           -11.400      9.68  17.4000  3 5.97e-04
## PMI                  1.910      8.92   0.0904  1 7.64e-01
## Age.Instit           1.180     15.20   2.8200  2 2.44e-01
## Age.RNA_batch        2.350     16.40   1.6500  2 4.38e-01
## Age                  1.800      8.81   0.2020  1 6.53e-01
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
signif(summarize.anova(av.2), digits = 3)
```

```
##             Delta.AIC Delta.BIC   Chisq df    p.Chi
## Ances1         -0.246      6.77  2.2500  1 1.34e-01
## Ances2          1.670      8.68  0.3320  1 5.65e-01
## Ances3          1.920      8.93  0.0832  1 7.73e-01
## Ances4          0.864      7.88  1.1400  1 2.87e-01
## Ances5          1.960      8.97  0.0410  1 8.40e-01
## Ances1.Gene   -56.800    -21.70 66.8000  5 4.83e-13
## Ances2.Gene     5.090     40.20  4.9100  5 4.27e-01
## Ances3.Gene    -3.750     31.30 13.7000  5 1.73e-02
## Ances4.Gene     7.850     42.90  2.1500  5 8.28e-01
## Ances5.Gene     9.070     44.10  0.9260  5 9.68e-01
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
signif(summarize.anova(av.3), digits = 3)
```

```
##             Delta.AIC Delta.BIC    Chisq df   p.Chi
## Gender.fix      1.020      8.03 9.85e-01  1 0.32100
## Gender.ran      2.000      9.01 7.49e-10  1 1.00000
## Gender.Gene    -5.800      1.22 7.80e+00  1 0.00524
## Dx.fix          3.890     17.90 1.07e-01  2 0.94800
## Dx.ran          2.000      9.01 1.19e-09  1 1.00000
## Dx.Gene         0.453      7.47 1.55e+00  1 0.21400
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
signif(summarize.anova(list(av.4)), digits = 3)
```

```
##   Delta.AIC Delta.BIC Chisq df    p.Chi
## 1     -23.4      4.65  31.4  4 2.53e-06
```

### Saving model formulas

Save all model formulas in the `results` directory


```r
write.csv(data.frame(lapply(list(M1 = M1, M1b = M1b, M2 = M2, M3 = M3, M4 = M4),
                            function(x) as.character(formula(x)))),
          file = "../../results/M-formulas.csv", row.names = FALSE)
```

