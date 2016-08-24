## Data preparation



Load data importer functions:

```r
source("../../src/import-data.R")
source("../../src/fit-glms.R")
source("2016-08-21-likelihood-surface.R")
```

```
## Warning in file(filename, "r", encoding = encoding): cannot open file
## '2016-08-21-likelihood-surface.R': No such file or directory
```

```
## Error in file(filename, "r", encoding = encoding): cannot open the connection
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
