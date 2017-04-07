## Prepare data





## Dependence of $S$ on certain variables

### Dependence on gene, age, and institution

Implementation of the same plot both with the `lattice` and the `ggplot2` package.


```r
P <- list()
# lattice implementation
P$s.age.inst$lattice <-
    xyplot(S ~ Age | Gene, data = Y.long,
           subset = Gene %in% gene.ids,
           groups = Institution,
           panel = function(x, y, ...) {
               panel.xyplot(x, y, pch = 21, cex = 0.3, ...)
               panel.smoother(x, y, col = "black", lwd = 2, ...)
           },
           auto.key = list(title = "institution", columns = 3),
           par.settings = list(add.text = list(cex = 0.8)),
           ylab = "read count ratio, S",
           xlab = "age",
           aspect = "fill", layout = c(6, 5))
# ggplot2 implementation
g <- ggplot(data = Y.long, aes(x = Age, y = S))
g <- g + geom_point(pch = "o", aes(color = Institution))
g <- g + geom_smooth(method = "loess", color = "black")
g <- g + facet_wrap(~ Gene)
P$s.age.inst$ggplot2 <- g
plot(P$s.age.inst$lattice)
```

<img src="figure/S-age-smooth-1.png" title="plot of chunk S-age-smooth" alt="plot of chunk S-age-smooth" width="700px" />

```r
plot(P$s.age.inst$ggplot2)
```

<img src="figure/S-age-smooth-2.png" title="plot of chunk S-age-smooth" alt="plot of chunk S-age-smooth" width="700px" />

### Gene, age, and gender


```r
P <- list()
# lattice implementation
P$s.age.gender$lattice <-
    xyplot(S ~ Age | Gene, data = Y.long, groups = Gender,
           subset = Gene %in% gene.ids,
           panel = function(x, y, ...) {
               panel.xyplot(x, y, pch = 21, ...)
               #panel.smoother(x, y, col = "black", lwd = 2, ...)
           },
           par.settings = list(add.text = list(cex = 0.8),
                               superpose.symbol = list(cex = 0.5,
                                                       fill = trellis.par.get("superpose.symbol")$fill[c(2, 1)],
                                                       col = trellis.par.get("superpose.symbol")$col[c(2, 1)])),
           auto.key = list(title = "gender", columns = 2),
           ylab = "read count ratio, S",
           xlab = "age",
           aspect = "fill", layout = c(6, 5))
# ggplot2 implementation
g <- ggplot(data = Y.long, aes(x = Age, y = S))
g <- g + geom_point(pch = "o", aes(color = Gender))
g <- g + geom_smooth(method = "loess", color = "black")
g <- g + facet_wrap(~ Gene)
P$s.age.gender$ggplot2 <- g
plot(P$s.age.gender$lattice)
```

<img src="figure/S-age-gender-1.png" title="plot of chunk S-age-gender" alt="plot of chunk S-age-gender" width="700px" />

```r
#plot(P$s.age.gender$ggplot2)
```


```r
P$s.age$lattice <-
    xyplot(S ~ Age | Gene, data = Y.long,
           subset = Gene %in% gene.ids,
           par.settings = list(add.text = list(cex = 0.8),
                               strip.background = list(col = "gray90"),
                               plot.symbol = list(pch = 21, cex = 0.5, col = "black", fill = "gray", alpha = 0.5)),
           auto.key = list(title = "gender", columns = 2),
           ylab = "read count ratio, S",
           xlab = "age",
           aspect = "fill", layout = c(6, 5))
plot(P$s.age$lattice)
```

<img src="figure/S-age-1.png" title="plot of chunk S-age" alt="plot of chunk S-age" width="700px" />


```r
update(P$s.age$lattice, layout = c(5, 6))
```

<img src="figure/S-age-b-1.png" title="plot of chunk S-age-b" alt="plot of chunk S-age-b" width="700px" />

<img src="figure/R-age-gender-1.png" title="plot of chunk R-age-gender" alt="plot of chunk R-age-gender" width="700px" />

<img src="figure/Q-age-gender-1.png" title="plot of chunk Q-age-gender" alt="plot of chunk Q-age-gender" width="700px" />

### Dependence on gene, age, and total read count $N$

<img src="figure/S-age-tot-read-count-1.png" title="plot of chunk S-age-tot-read-count" alt="plot of chunk S-age-tot-read-count" width="700px" /><img src="figure/S-age-tot-read-count-2.png" title="plot of chunk S-age-tot-read-count" alt="plot of chunk S-age-tot-read-count" width="700px" />

## Associations between explanatory variables

### Deterministic association: RIN and RIN2

<img src="figure/rin-rin2-1.png" title="plot of chunk rin-rin2" alt="plot of chunk rin-rin2" width="300" /><img src="figure/rin-rin2-2.png" title="plot of chunk rin-rin2" alt="plot of chunk rin-rin2" width="300" />

### Stochastic (statistical) associations

Both "scatter plot matrices" show the same set of pairwise associations (top: `lattice`, bottom: `ggplot2` and `GGally` packages).

<img src="figure/evar-scatterplot-matrix-1.png" title="plot of chunk evar-scatterplot-matrix" alt="plot of chunk evar-scatterplot-matrix" width="700px" /><img src="figure/evar-scatterplot-matrix-2.png" title="plot of chunk evar-scatterplot-matrix" alt="plot of chunk evar-scatterplot-matrix" width="700px" />
