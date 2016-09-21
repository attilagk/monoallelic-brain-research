
```r
library(ggplot2)
library(lattice)
library(latticeExtra)
opts_chunk$set(dpi = 144)
opts_chunk$set(out.width = "700px")
opts_chunk$set(dev = c("png", "pdf"))
lattice.options(default.args = list(as.table = TRUE))
lattice.options(default.theme = "standard.theme")
```

Load functions


```r
source("~/projects/monoallelic-brain/src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
```


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
```

Get data: observations on predictors (explanatory variables) and on the higher and lower read count from selected genes (more details in a previous post)

```r
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
```


```r
set.seed(1976)
perm.obs <- sample.int(nrow(E))
names(e.vars) <- e.vars
EP <- lapply(e.vars, function(v) { E1 <- E; E1[[v]] <- E[[v]][perm.obs]; return(E1) })
EP$Unpermuted <- E
```


```r
# exclude unweighed aggregates UA.8 and UA from fitting
to.fit.ids <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
M <- lapply(EP, function(e) do.all.fits(Y[to.fit.ids], G = e, preds = e.vars, sel.models = c("logi.S", "wnlm.R")))
```


```r
Betas.Unpermuted <- lapply(M$Unpermuted,
                           function(m) { x <- get.estimate.CI(m); x <- x[ ! x$Coefficient %in% "(Intercept)", ] })
```


```r
Betas.Permuted.logi.S <- lapply(e.vars,
                                function(e) {
                                    x <- get.estimate.CI(M[[e]][["logi.S"]])
                                    x <- x[ grep(e, row.names(x), value = TRUE), ]
                                    #x <- x[ ! x$Coefficient %in% "(Intercept)", ]
                                })
```


```r
Betas.Permuted.wnlm.R <- lapply(e.vars,
                                function(e) {
                                    x <- get.estimate.CI(M[[e]][["wnlm.R"]])
                                    x <- x[ grep(e, row.names(x), value = TRUE), ]
                                    #x <- x[ ! x$Coefficient %in% "(Intercept)", ]
                                })
```


```r
my.segplot(data = Betas.Unpermuted$logi.S, xlim = my.xlim, main = "Unpermuted under logi.S")
```

<img src="figure/unpermuted-logi-S-1.png" title="plot of chunk unpermuted-logi-S" alt="plot of chunk unpermuted-logi-S" width="700px" />


```r
my.segplot(data = do.call(rbind, Betas.Permuted.logi.S), xlim = my.xlim, main = "Permuted under logi.S")
```

<img src="figure/permuted-logi-S-1.png" title="plot of chunk permuted-logi-S" alt="plot of chunk permuted-logi-S" width="700px" />


```r
my.segplot(data = Betas.Unpermuted$wnlm.R, main = "Unpermuted under wnlm.R")
```

<img src="figure/unpermuted-wnlm-R-1.png" title="plot of chunk unpermuted-wnlm-R" alt="plot of chunk unpermuted-wnlm-R" width="700px" />


```r
my.segplot(data = do.call(rbind, Betas.Permuted.wnlm.R), main = "Permuted under wnlm.R")
```

<img src="figure/permuted-wnlm-R-1.png" title="plot of chunk permuted-wnlm-R" alt="plot of chunk permuted-wnlm-R" width="700px" />

