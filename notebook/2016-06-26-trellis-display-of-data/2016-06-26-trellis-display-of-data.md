## Prepare data


```r
library(lattice)
library(latticeExtra)
library(ggplot2)
source("../../src/import-data.R")
source("../../src/fit-glms.R")
```


```r
E <- get.predictors()
Y <- get.readcounts(gene.ids) # gene.ids is predifined
S.long <- reshape(cbind(lapply(Y, getElement, "S"), E[ e.vars ]), # the data to be reshaped; e.vars is predifined
                  v.names = "S", varying = list(names(Y)), # the component that varies with genes ('time')
                  timevar = "Gene", times = names(Y), # the 'time' variable, i.e. gene symbol
                  ids = rownames(E), idvar = "sample.id", # observation names: RNA sample ids
                  direction = "long")
S.long$Gene <- factor(S.long$Gene, levels = names(Y), ordered = TRUE)
# in what follows exclude the unweighted averages 'UA.8' and 'UA' because total read count 'N' is not available for those
Y.long <- cbind(stack(lapply(Y[ c(gene.ids, "WA.8", "WA") ], getElement, "S")),
                stack(lapply(Y[ c(gene.ids, "WA.8", "WA") ], getElement, "N")),
                E[ rep(seq_len(nrow(E)), length(gene.ids) + 2), e.vars ])
names(Y.long)[1:3] <- c("S", "Gene", "N")
Y.long$Gene <- factor(Y.long$Gene, levels = names(Y)[ 0:1 - length(Y) ], ordered = TRUE)
Y.long[[4]] <- NULL
```




```r
plot(plot.s.age.inst)
```

![plot of chunk S-age-smooth](figure/S-age-smooth-1.png)



```r
plot(plot.s.age.n)
```

```
## Warning: Removed 6319 rows containing missing values (geom_point).
```

![plot of chunk S-age-tot-read-count](figure/S-age-tot-read-count-1.png)


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
