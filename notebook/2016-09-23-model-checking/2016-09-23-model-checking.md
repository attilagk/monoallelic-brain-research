

Load functions


```r
source("~/projects/monoallelic-brain/src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
```


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
```

Get data

```r
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = gene.ids, count.thrs = 0)
```

Import outliers of the CMC study provided by Mette and Jessica and check if any has been included in our previous analyses.  The answer is: they have all been excluded.


```r
read.csv("../../data/outliers.csv", colClasses = "character")$DissectionID %in% row.names(E)
```

```
##  [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
## [12] FALSE
```

Fit logi.S and wnlm.R models


```r
# exclude unweighed aggregates UA.8 and UA from fitting
to.fit.ids <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
sel.models <- c("logi.S", "wnlm.R"); names(sel.models) <- sel.models
M <- do.all.fits(Y[to.fit.ids], G = E, preds = e.vars, sel.models = sel.models)
```

