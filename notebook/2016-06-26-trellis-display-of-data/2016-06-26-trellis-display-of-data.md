## Prepare data


```r
library(lattice)
library(latticeExtra)
library(ggplot2)
library(GGally)
source("../../src/import-data.R")
source("../../src/fit-glms.R")
```



## Dependence of $S$ on certain variables

### Dependence on gene, age, and institution

Implementation of the same plot both with the `lattice` and the `ggplot2` package.


```r
P <- list()
# lattice implementation
P$s.age.inst$lattice <-
    xyplot(S ~ Age | Gene, data = S.long, groups = Institution,
           panel = function(x, y, ...) {
               panel.xyplot(x, y, pch = 21, ...)
               panel.smoother(x, y, col = "black", lwd = 2, ...)
           },
           auto.key = list(title = "Insitution", space = "right"),
           aspect = "fill", layout = c(6, 5))
# ggplot2 implementation
g <- ggplot(data = S.long, aes(x = Age, y = S))
g <- g + geom_point(pch = "o", aes(color = Institution))
g <- g + geom_smooth(method = "loess", color = "black")
g <- g + facet_wrap(~ Gene)
P$s.age.inst$ggplot2 <- g
plot(P$s.age.inst$lattice)
```

![plot of chunk S-age-smooth](figure/S-age-smooth-1.png)

```r
plot(P$s.age.inst$ggplot2)
```

![plot of chunk S-age-smooth](figure/S-age-smooth-2.png)

### Dependence on gene, age, and total read count $N$

![plot of chunk S-age-tot-read-count](figure/S-age-tot-read-count-1.png)![plot of chunk S-age-tot-read-count](figure/S-age-tot-read-count-2.png)


```r
M <- do.all.fits(Y[ 0:1 - length(Y) ], # omit the last two components: "UA.8" and "UA"
                 E, preds = e.vars, sel.models = c("wnlm.R", "logi.S"),
                 x = TRUE, y = TRUE) # store model matrix 'x' and response 'y'
```

```
## Warning: glm.fit: algorithm did not converge
```

```
## Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
```

## Associations between explanatory variables

![plot of chunk rin-rin2](figure/rin-rin2-1.png)![plot of chunk rin-rin2](figure/rin-rin2-2.png)

![plot of chunk evar-scatterplot-matrix](figure/evar-scatterplot-matrix-1.png)![plot of chunk evar-scatterplot-matrix](figure/evar-scatterplot-matrix-2.png)

### Regression coefficients


```r
Betas <- lapply(M, get.estimate.CI)
A.long <- lapply(M, function(m) reshape.1(l.anova(m), type = "anova"))
Ef.long <- lapply(M, function(m) reshape.1(l.effects(m), type = "effects"))
```

![plot of chunk reg-coef-logi.S](figure/reg-coef-logi.S-1.png)

![plot of chunk reg-coef-wnlm.R](figure/reg-coef-wnlm.R-1.png)

![plot of chunk anova-logi.S](figure/anova-logi.S-1.png)

![plot of chunk anova-wnlm.R](figure/anova-wnlm.R-1.png)

![plot of chunk effects-logi.S](figure/effects-logi.S-1.png)

![plot of chunk effects-wnlm.R](figure/effects-wnlm.R-1.png)
