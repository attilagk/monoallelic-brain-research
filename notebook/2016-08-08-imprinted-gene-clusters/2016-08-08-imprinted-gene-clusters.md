## Motivation

Imprinted genes form **clusters** of one or more genes.  Previous work by Ifat established the **imprinting status** of each gene as either known to be imprinted, *not* known to be imprinted but near an "known" gene, or else neither.  I refer to these three categories as "known", "candidate, <1MB" and "candidate", respectively.

The main question queries the mechanism of the age effect on imprinting (loss, gain or lack of effect).  In particular: does age regulate genes within some cluster in a concerted or an independent manner?

The current work studies this by performing two main steps

1. delineation of clusters
1. visualize clusters and age effect in terms of regression coefficients for age

## Data import and preparation



Load data importer functions:

```r
source("../../src/import-data.R")
source("../../src/utils.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
source("2016-08-08-imprinted-gene-clusters.R")
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
# imprinting segments in component 'seg': clusters and segments inbetween
gs.seg <- make.impr.segs(gene.summary)
# remove filtered genes
gs <- gs.seg[names(frac), ]
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

## Analysis

### Genomic location

The plot below shows the genomic location of all 16026 genes in the filtered data set and indicates their imprinting status with different colors.  Also shown is the gene score according to which genes have been ranked and called as monoallelically expressing or not.


```r
t.par <- list(superpose.symbol = list(pch = c(21, 21, 21), alpha = c(0.3, 1, 1), fill = c("pink", "green", "lightblue"), col = c("red", "darkgreen", "blue")))
trellis.par.set(t.par)
xyplot(score ~ start | chr, data = gs, groups = imprinting.status, auto.key = list(columns = 3), layout = c(4, 6),
       xlab = "genomic location", ylab = "score: 1 - ECDF at s = 0.9")
```

![plot of chunk score-genomic-location](figure/score-genomic-location-1.png)

The next plot presents the maximum likelihood estimate $\hat{\beta}_\mathrm{age}$ of the regression coefficient mediating age's effect in the logistic model `logi.S`.  Only those genes are shown that were fitted, of course.  Confidence intervals are not shown for clarity.


```r
beta.age <- get.CI(M$logi.S[f.ids$logi.S], coef.name = "Age", conf.lev = 0.99)
gs$beta.hat <- NA
gs[rownames(beta.age), "beta.hat"] <- beta.age$beta.hat
xyplot(beta.hat ~ start | chr, data = gs, groups = imprinting.status, auto.key = list(columns = 3), layout = c(4, 6), panel = function(...) { panel.abline(h = 0, col = trellis.par.get("reference.line")$col); panel.xyplot(...) })
```

![plot of chunk beta-genomic-location](figure/beta-genomic-location-1.png)

### Delineation of imprinted gene clusters

Let $K$ denote the set of genes $g$ such that $g\in K$ means the gene is either "known" (to be imprinted) or near some "known", i.e. "candidate, <1MB". Let $C$ be the set of "candidate" genes (further away from "known"s). Let $g_{i-1}$ and $g_i$ be neighboring genes on the same chromosome at the $i-1$-th and $i$-th site.  Delineation of clusters was done using the following simple rule: if $g_{i-1}\in C$ and $g_{i}\in K$ then a cluster starts at the $i$-th gene.  The cluster ends at the first $j\gt i$ such that $g_{j}\in K$ but $g_{i}\in C$.  The rule is implemented in the `make.impr.segs` function in `2016-08-08-imprinted-gene-clusters.R`.

With that rule 36 clusters have been defined:

![plot of chunk cluster-sizes](figure/cluster-sizes-1.png)



### Clusters and the age effect

After some uninteresting data manipulation (code hidden) **the main result** can be presented:
![plot of chunk segplot](figure/segplot-1.png)

## Conclusion

The last, main, result suggests that

1. age operates on genes within the same cluster independently rather than in concert
1. however, many results are "weak":
   * relatively low power indicated by several extremely wide confidence intervals
   * the fit of the logistic model is not very satisfactory (discussed earlier, see conditional analysis)
