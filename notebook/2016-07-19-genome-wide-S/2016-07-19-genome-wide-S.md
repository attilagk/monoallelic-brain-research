### Missing data due to bug

For each gene $g$ of 16026 genes the higher and lower read counts across $I_g\le 579$ individuals $\{H_{ig}, L_{ig} \,|\, i=1,...,I_g\}$ are stored in a `.csv` file.  However, nearly 3 % of those files are empty due to a bug that resulted in multiple headers in [Ifat's html tables][ifat] (see for instance gene AK3) and to the inability of my converter script to deal with such tables.  To see this:


```bash
DATADIR="../../data/readcount/"
SEDCMD='s/^.*attila attila\s\+\([[:digit:]]\+\).*\/\([^/]\+\)\.csv$/\1,\2/'
OUTFILE="${DATADIR}/fsize-genes.csv"
[ -f $OUTFILE ] || {
    # create header
    echo "file.size,gene.symbol" > $OUTFILE
    # get file size (in bytes) for each gene under 'genes' directory
    find "${DATADIR}/genes" -name '*.csv' | xargs ls -l | sed "$SEDCMD" >> $OUTFILE
}
cat <<EOF
Number of .csv files

all: $(grep --count '^[[:digit:]]\+,' $OUTFILE )
empty: $(grep --count '^0,' $OUTFILE )
EOF
```

```
## Number of .csv files
## 
## all: 16026
## empty: 442
```



### Import read count data

Load data importer functions:

```r
source("../../src/import-data.R")
source("2016-07-19-genome-wide-S.R")
```

The following expressions import $S_{ig}$ for all 1.5584 &times; 10<sup>4</sup> genes for which the csv file is nonempty.


```r
gene.summary <-
    read.csv("../../data/readcount/summary-all-genes.csv")
rownames(gene.summary) <- gene.summary$Symbol
gene.summary$file.size <-
    read.csv("../../data/readcount/fsize-genes.csv", row.names = 2)[rownames(gene.summary), , drop = TRUE]
gene.ids <- with(gene.summary, as.character(Symbol)[ file.size > 0 ])
```


```r
Y <- get.readcounts(gene.ids, g.subsets = list(), sel.var = "S")
Y <- data.frame(lapply(Y, getElement, "S"))
```

### Filtering genes based on number of observations

For more than half of even the genes $g$ with nonempty files the number $I_g$ of observations (the number of individuals/RNA samples with read count data on $g$) is zero.  In what follows, not only these genes are filtered out but also those with less than 10 observations, indicated by the vertical dashed line on the empirical ECDF plot below.


```r
min.n.obs <- 10
```
![plot of chunk ecdf-n-obs](figure/ecdf-n-obs-1.png)

### Data preparation

Computationally demanding calculations to prepare data for presentation:


```r
# filter genes given the minimum number of allowed observations 'min.n.obs'
ok.genes <- names(Y)[sapply(Y, function(y) sum(! is.na(y)) >= min.n.obs)]
# obtain the ECDF of S_ig for all given genes g
ED <- emp.distr.S(Y[ , ok.genes],
                  ss = seq(from = 0.5, to = 1, length.out = 101),
                  with.density = TRUE)
# order genes according to the value of the ECDF at 0.9
gene.order <- order(sapply(ED[[1]], function(f) f(0.9)))
ecdf.val.w <- ED$ecdf.val[ , gene.order]
density.w <- ED$density[ , gene.order]
ED.long <- reshape(ecdf.val.w, v.names = "ECDF", varying = names(ecdf.val.w),
                   timevar = "gene", times = factor(names(ecdf.val.w)),
                   idvar = "s", ids = ED$ss, direction = "long")
density.long <- reshape(density.w, v.names = "density", varying = names(density.w),
                   timevar = "gene", times = factor(names(density.w)),
                   idvar = "s", ids = ED$ss, direction = "long")
ED.long$density <- density.long$density
ED.long$gene <- factor(ED.long$gene, levels = names(ecdf.val.w), ordered = TRUE)
rm(list = c("ecdf.val.w", "density.w", "density.long"))
```





### Figure for manuscript

This figure is intended to:

* present a few representative genes characterized by the parental "imbalance score" $S$ statistic (two upper graphs)
  * e.g. PEG10 and ZNF331 represent monoallelic expression and AFAP1 biallelic expression
  1. density est.: kernel density estimates, whose interpretation is identical to that of histograms
  1. ECDF: empirical cumulative distribution function, sometimes called cumulative fraction
* present a compact genome-wide overview of parental imbalance (lower left graph)
  * imbalance is expressed in terms of the ECDF of $S$ using a color scheme, introduced in the upper 2nd and 3rd graphs
  * the thousands of genes are ranked from the most imbalanced (monoallelic, top) to the most balanced (biallelic, bottom)
  * the ranking is based on the ECDF evaluated at $s = 0.9$, as shown in the upper 2nd graph and the bottom right plots
* support the conclusion that $\approx 1 \%$ of all genes are appreciably imbalanced (monoallelically expressed)

![plot of chunk complex-plot](figure/complex-plot-1.png)

### Fixed ranking bug

An earlier version of the above figure had incorrect ranking.  The figure compares that to the present, correct, ranking:

![plot of chunk fixed-ranking](figure/fixed-ranking-1.png)

### The top ranking genes

The plot below shows the top 50 genes using the same ranking as the one above

![plot of chunk rank-by-ecdf-top](figure/rank-by-ecdf-top-1.png)

The next plot is meant to serve a consistency check with Ifat's corresponding plot, that is [Fig 1 of the previous manuscript][Fig 1]

![plot of chunk compare-to-ifats-fig](figure/compare-to-ifats-fig-1.png)



[ifat]: http://katahdin.mssm.edu/ifat/web/cm/home
[Fig 1]: https://docs.google.com/presentation/d/1YvpA1AJ-zzir1Iw0F25tO9x8gkSAzqaO4fjB7K3zBhE/edit#slide=id.p4
