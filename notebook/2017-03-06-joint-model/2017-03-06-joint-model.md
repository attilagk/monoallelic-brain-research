Hello, World!




```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
E <- get.predictors()[e.vars]
Y <- get.readcounts(gene.ids)[gene.ids]
dat <- do.call(rbind,
               lapply(gene.ids,
                      function(g) cbind(y <- Y[[g]], data.frame(id = rownames(y), Gene = g), E)))
```
