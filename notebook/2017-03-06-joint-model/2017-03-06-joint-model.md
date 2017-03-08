

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
summarize.anova(av)
```

```
##               Delta.AIC    Delta.BIC      Chisq df         p.Chi
## RIN            -8.03994    -1.026466   10.03994  1  1.531822e-03
## RNA_batch     -25.69463   -18.681158   27.69463  1  1.420564e-07
## Institution   -61.18655   -54.173078   63.18655  1  1.880277e-15
## Individual   -453.77916  -446.765691  455.77916  1 3.984989e-101
## Age.Gene     -192.22905  -178.202103  196.22905  2  2.451338e-43
## Gene        -8700.32456 -8693.311084 8702.32456  1  0.000000e+00
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
summarize.anova(av.2)
```

```
##               Delta.AIC  Delta.BIC       Chisq df        p.Chi
## Ances1       -0.2462851   6.767188  2.24628508  1 1.339356e-01
## Ances2        1.6684863   8.681960  0.33151372  1 5.647691e-01
## Ances3        1.9168282   8.930302  0.08317178  1 7.730443e-01
## Ances4        0.8643438   7.877817  1.13565620  1 2.865721e-01
## Ances5        1.9590356   8.972509  0.04096436  1 8.396067e-01
## Ances1.Gene -56.7626802 -21.695312 66.76268018  5 4.826726e-13
## Ances2.Gene   5.0936690  40.161037  4.90633097  5 4.274184e-01
## Ances3.Gene  -3.7489132  31.318455 13.74891321  5 1.728553e-02
## Ances4.Gene   7.8501326  42.917500  2.14986743  5 8.280435e-01
## Ances5.Gene   9.0735245  44.140892  0.92647551  5 9.682754e-01
```

#### Notable results

* ancestry components are not equally important: only *Ancestry.1* and *Ancestry.3* appears to matter
* their main effect is weaker than their interaction with *Gene*
    * this is expected: the effect of *Age* should not depend on e.g. the *Institution*; so the result shows the advantage of joint modeling because previous results with separate, gene-based, modeling suggested strong interaction between *Age* and *Institution*


```r
M3 <- update(M1b, . ~ . + (scale(Age) + scale(RIN) + scale(PMI) +
                           scale(Ancestry.1) + scale(Ancestry.3) | Gene) +
             Ancestry.1)
```


```r
av.3 <- list()
# Gender
av.3$Gender <- anova(M3, update(M3, . ~ . + Gender))
av.3$Gender.Gene <-
    anova(M3, update(M1b, . ~ . + (Gender + scale(Age) + scale(RIN) + scale(PMI) +
                                   scale(Ancestry.1) + scale(Ancestry.3) | Gene) + Ancestry.1))
```

```
## Warning in commonArgs(par, fn, control, environment()): maxfun < 10 *
## length(par)^2 is not recommended.
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

```
## Warning in commonArgs(par, fn, control, environment()): maxfun < 10 *
## length(par)^2 is not recommended.
```

```r
# Dx
av.3$Dx <- anova(M3, update(M3, . ~ . + Dx))
av.3$Dx.Gene <-
    anova(M3, update(M1b, . ~ . + (Dx + scale(Age) + scale(RIN) + scale(PMI) +
                                   scale(Ancestry.1) + scale(Ancestry.3) | Gene) + Ancestry.1))
```

```
## Warning in commonArgs(par, fn, control, environment()): maxfun < 10 *
## length(par)^2 is not recommended.
```

```
## Warning in optwrap(optimizer, devfun, getStart(start, rho$lower, rho$pp), :
## convergence code 1 from bobyqa: bobyqa -- maximum number of function
## evaluations exceeded
```

```
## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control
## $checkConv, : Model failed to converge with max|grad| = 0.012066 (tol =
## 0.002, component 1)
```

```
## Warning in commonArgs(par, fn, control, environment()): maxfun < 10 *
## length(par)^2 is not recommended.
```

```
## Warning in optwrap(optimizer, devfun, x@theta, lower = x@lower, calc.derivs
## = TRUE, : convergence code 1 from bobyqa: bobyqa -- maximum number of
## function evaluations exceeded
```


```r
summarize.anova(av.3)
```

```
##              Delta.AIC  Delta.BIC      Chisq df      p.Chi
## Gender       0.9596084   7.973082  1.0403916  1 0.30773042
## Gender.Gene -1.9494366  47.144878 15.9494366  7 0.02558235
## Dx           3.8394762  17.866423  0.1605238  2 0.92287461
## Dx.Gene     17.6465272 122.848630 12.3534728 15 0.65210083
```
