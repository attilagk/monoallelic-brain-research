### Import read count data

Load data importer functions:

```r
source("../../src/import-data.R")
source("2016-07-28-test-expr-balance.R")
```

Import $S_{ig}$ for all 1.5584 &times; 10<sup>4</sup> genes for which the csv file is nonempty (see previous post).


```r
gene.summary <-
    read.csv("../../data/readcount/summary-all-genes.csv")
rownames(gene.summary) <- gene.summary$Symbol
gene.summary$file.size <-
    read.csv("../../data/readcount/fsize-genes.csv", row.names = 2)[rownames(gene.summary), , drop = TRUE]
gene.ids <- with(gene.summary, as.character(Symbol)[ file.size > 0 ])
Y <- get.readcounts(gene.ids, g.subsets = list(), sel.var = c("S", "N"))
S <- data.frame(lapply(Y, getElement, "S"))
N <- data.frame(lapply(Y, getElement, "N"))
```


```r
min.n.obs <- 10
ok.genes <- names(S)[sapply(S, function(y) sum(! is.na(y)) >= min.n.obs)]
ED <- lapply(S[ok.genes], ecdf)
gene.order <- order(sapply(ED, function(f) f(0.9)))
sorted.genes <- as.list(ok.genes[gene.order])
names(sorted.genes) <- sorted.genes
cum.frac.obs <- data.frame(lapply(ED[names(sorted.genes)], function(f) f(10:6 / 10)))
row.names(cum.frac.obs) <- as.character(10:6 / 10)
andys.test <-
    data.frame(lapply(sorted.genes,
                      function(g)
                          sum(S[[g]] <= 0.6 & CI.p(S[[g]], N[[g]])$upper < 0.7, na.rm = TRUE) / sum(! is.na(S[[g]]))))
row.names(andys.test) <- "andys.test"
cum.frac.obs <- rbind(cum.frac.obs, andys.test)
rm(andys.test)
```


```r
#
frac.obs <- data.frame(lapply(cum.frac.obs, function(y) - diff(c(y, 0))))
```
