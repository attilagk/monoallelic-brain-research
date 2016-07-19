Load data importer functions:

```r
source("../../src/import-data.R")
```

Finally!


```
## Loading required package: RColorBrewer
```


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


```r
gene.summary <-
    read.csv("../../data/readcount/summary-all-genes.csv")
rownames(gene.summary) <- gene.summary$Symbol
gene.summary$file.size <-
    read.csv("../../data/readcount/fsize-genes.csv", row.names = 2)[rownames(gene.summary), , drop = TRUE]
gene.ids <- with(gene.summary, as.character(Symbol)[ file.size > 0 ])
```

### Import read count data


```r
Y <- get.readcounts(gene.ids, g.subsets = list(), sel.var = "S")
Y <- data.frame(lapply(Y, getElement, "S"))
```


```r
min.n.obs <- 10
```
![plot of chunk ecdf-n-obs](figure/ecdf-n-obs-1.png)



![plot of chunk ecdf-levelplot](figure/ecdf-levelplot-1.png)

[ifat]: http://katahdin.mssm.edu/ifat/web/cm/home
