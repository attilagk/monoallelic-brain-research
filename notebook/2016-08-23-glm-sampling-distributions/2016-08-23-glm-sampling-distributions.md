## Data preparation


```
## Loading required package: RColorBrewer
```

Load data importer functions:

```r
source("../../src/import-data.R")
source("../../src/fit-glms.R")
source("2016-08-23-glm-sampling-distributions.R")
```

Do the import:


```r
genes2fit <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(genes2fit) <- genes2fit
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = genes2fit, count.thrs = 0)
```

Fit models:


```r
# exclude unweighed aggregates UA.8 and UA from fitting
ids2fit <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
M <- list()
f.ids <- list()
M$multi <- do.all.fits(Y[ids2fit], preds = e.vars)
f.ids$multi <- as.data.frame(lapply(M$multi, function(m) ! sapply(m, is.null)))
f.ids$multi["TMEM261P1", c("logi.S", "logi2.S")] <- FALSE
f.ids$multi[c("WA.8", "WA"), ] <- FALSE # gene aggregates are not needed here
M$simple <- do.all.fits(Y[ids2fit], preds = "Age")
f.ids$simple <- as.data.frame(lapply(M$simple, function(m) ! sapply(m, is.null)))
f.ids$simple["TMEM261P1", c("logi.S", "logi2.S")] <- FALSE
f.ids$simple[c("WA.8", "WA"), ] <- FALSE # gene aggregates are not needed here
```

## Results


```r
plots <- list()
plots$GRB10 <- plot.mdls.1gene(gene = "GRB10", l.l.M = M$simple) 
plots$PEG3 <- plot.mdls.1gene(gene = "PEG3", l.l.M = M$simple) 
```

<img src="figure/GRB10-1.png" title="plot of chunk GRB10" alt="plot of chunk GRB10" width="700px" />

<img src="figure/PEG3-1.png" title="plot of chunk PEG3" alt="plot of chunk PEG3" width="700px" />
