## Motivation

Imprinted gene **clusters**

* delineation of clusters
* clusters and age effect (gain, loss or no change)

## Data import and preparation



Load data importer functions:

```r
source("../../src/import-data.R")
source("../../src/utils.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
```

### Genome-wide data


```r
gene.summary <-
    read.csv("../../data/readcount/summary-all-genes.csv")
rownames(gene.summary) <- gene.summary$Symbol
gene.summary$file.size <-
    read.csv("../../data/readcount/fsize-genes.csv", row.names = 2)[rownames(gene.summary), , drop = TRUE]
gene.ids <- with(gene.summary, as.character(Symbol)[ file.size > 0 ])
genes2fit <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(genes2fit) <- genes2fit
```


```r
Y <- get.readcounts(gene.ids, g.subsets = list(), sel.var = c("H", "N"))
S <- data.frame(lapply(gene.ids, function(g) Y[[g]]$H / Y[[g]]$N), check.names = FALSE)
names(S) <- gene.ids
N <- data.frame(lapply(Y, getElement, "N"), check.names = FALSE)
rm(Y)
```

#### Filtering

1. read count-based filter with threshold $t_\mathrm{rc}=15$
1. individual-based filter with threshold $t_\mathrm{ind}=25$

The code was copied from an earlier post and is hidden here.



```r
gene.summary$imprinting.status <- factor(gene.summary$imprinted, ordered = TRUE)
levels(gene.summary$imprinting.status) <- c("candidate", "candidate, <1MB", "known")
gene.summary$Symbol <- factor(gene.summary$Symbol, levels = gene.summary$Symbol, ordered = TRUE)
gene.summary$chr <- factor(paste("chr", gene.summary$chr), levels = paste("chr", seq_along(levels(factor(gene.summary$chr)))), ordered = TRUE)

# sort according to location
gs.loc <- gene.summary[with(gene.summary, order(chr, start)), ]

# remove filtered genes
gs <- gene.summary[names(frac), ]
gs$score <- unlist(frac["1", ])
```

### Fitting models to selected genes

Import read count data but do **not** filter, to be consistent with the most recent regression analysis:

```r
genes2fit <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(genes2fit) <- genes2fit
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = genes2fit, count.thrs = 0)
```

Fitting all models to all retained gene-wise and aggregated read count data sets

```r
# exclude unweighed aggregates UA.8 and UA from fitting
ids2fit <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
M <- do.all.fits(Y[ids2fit], preds = e.vars)
f.ids <- as.data.frame(lapply(M, function(m) ! sapply(m, is.null)))
f.ids["TMEM261P1", c("logi.S", "logi2.S")] <- FALSE
f.ids[c("WA.8", "WA"), ] <- FALSE # gene aggregates are not needed here
```


```r
beta.age <- get.CI(M$logi.S[f.ids$logi.S], coef.name = "Age", conf.lev = 0.99)
segplot(Gene ~ Lower.CL + Upper.CL, data = beta.age, draw.bands = FALSE, centers = beta.hat, xlab = expression(beta[age]), main = expression(paste(hat(beta)[age], " and CI")))
```

![plot of chunk segplot](figure/segplot-1.png)

## Analysis


```r
t.par <- list(superpose.symbol = list(pch = c(21, 21, 21), alpha = c(0.3, 1, 1), fill = c("pink", "green", "lightblue"), col = c("red", "darkgreen", "blue")))
trellis.par.set(t.par)
xyplot(score ~ start | chr, data = gs, groups = imprinting.status, auto.key = list(columns = 3), layout = c(4, 6))
```

![plot of chunk score-genomic-location](figure/score-genomic-location-1.png)


```r
gs$beta.hat <- NA
gs[rownames(beta.age), "beta.hat"] <- beta.age$beta.hat
xyplot(beta.hat ~ start | chr, data = gs, groups = imprinting.status, auto.key = list(columns = 3), layout = c(4, 6), panel = function(...) { panel.abline(h = 0, col = trellis.par.get("reference.line")$col); panel.xyplot(...) })
```

![plot of chunk beta-genomic-location](figure/beta-genomic-location-1.png)
