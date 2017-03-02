Hello, World!


```r
library(lme4)
library(lattice)
lattice.options(default.args = list(as.table = TRUE))
lattice.options(default.theme = "standard.theme")
opts_chunk$set(dpi = 144)
opts_chunk$set(out.width = "700px")
opts_chunk$set(dev = c("png", "pdf"))
```


```r
xyplot(Reaction ~ Days | Subject, data = sleepstudy, layout = c(6, 3))
```

<img src="figure/unnamed-chunk-2-1.png" title="plot of chunk unnamed-chunk-2" alt="plot of chunk unnamed-chunk-2" width="700px" />


```r
fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
fm3 <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
```

The magenta lines represent the fitted line under `fm1`, the green lines under `fm3` (the observed data remain cyan).


```r
df <- cbind(sleepstudy, data.frame(yhat.fm1 = predict(fm1), yhat.fm3 = predict(fm3)))
xyplot(Reaction + yhat.fm1 + yhat.fm3 ~ Days | Subject, data = df, type = "l", ylab = "Reaction", layout = c(6, 3),
       auto.key = list(text = c("observed", "predicted by fm1", "predicted by fm3"), points = FALSE, lines = TRUE))
```

<img src="figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" width="700px" />

Likelihood ratio test comparing the `fm1` model to its constrained version `fm3`:


```r
anova(fm1, fm3)
```

```
## refitting model(s) with ML (instead of REML)
```

```
## Data: sleepstudy
## Models:
## fm3: Reaction ~ Days + (1 | Subject)
## fm1: Reaction ~ Days + (Days | Subject)
##     Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
## fm3  4 1802.1 1814.8 -897.04   1794.1                             
## fm1  6 1763.9 1783.1 -875.97   1751.9 42.139      2  7.072e-10 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Thus the data supports `fm1` highly significantly better than `fm3`.  Based on this we can conclude that the dependence of `Reaction` on `Days` is heterogeneous across `Subject`s.  In other words: the effect of `Days` depends on the `Subject`.
