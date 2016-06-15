## Preparing data for analysis in `R`

### Read counts from Ifat's website

The source of read count data is the largest HTML table generated dynamically using the `All genes` option in [Ifat's website][summary].  The output is a `csv` file for each gene within a gene-specific subdirectory.  Due to the specifics of Ifat's html tables and the general difficulty of the `html` to `csv` conversion only the "higher" and "lower" read counts could be formatted in a usable way; the underlying **SNP-specific read counts could not be retrieved** in a clean way.  The same holds for the corresponding genotypes.


```bash
# directories
PDIR=$HOME/projects/monoallelic-brain # project directory
RCDIR=$PDIR/data/readcount/ # directory for readcounts
OUTDIR=$RCDIR/genes # subdir for gene-wise readcount

# data files
SUMMARYHTML=$RCDIR/summary-all-genes.html # should already exist
SUMMARYCSV=$RCDIR/summary-all-genes.csv # will be created

# script files
SCRIPT1=$PDIR/src/summary-html2csv
SCRIPT2=$PDIR/src/all-genes-readcounts-for-R 

# run scripts
test -f $SUMMARYHTML || exit 1
test -f $SUMMARYCSV || $SCRIPT1 $SUMMARYHTML > $SUMMARYCSV
test -d $OUTDIR || $SCRIPT2 $SUMMARYCSV 2> /dev/null
echo "number of gene-wise read count datasets:"
find $OUTDIR -type d | wc -l
echo "total size:"
du -sh $OUTDIR

```

```
## number of gene-wise read count datasets:
## 16027
## total size:
## 1.4G	/home/attila/projects/monoallelic-brain/data/readcount//genes
```

### Predictors

The data on predictors (explanatory variables, which can be either covariates or factors) was sent to me by Andy as email attached files.  The following code produces two files: `predictors.csv` is the main data table whereas `RNAseq_ID.DLPFC_RNA_ID.csv` is a mapping between RNAseqID (present in read count tables) and DLPFC_RNA_ID (present in the table of predictors).  Putting together those two is much easier in `R` than in shell scripts.


```bash
# directories
PDIR=$HOME/projects/monoallelic-brain # project directory
INDIR=$PDIR/data/ifat/age-dependence # input directory
OUTDIR=$PDIR/data/predictors/ # output directory

# input files from Ifat via Andy
INPRED="$INDIR/DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv"
INID="$INDIR/samples.csv"

# script files
SCRIPT=$PDIR/src/predictors-for-R

# run scripts
test -d $OUTDIR || $SCRIPT $INPRED $INID $OUTDIR
echo "output files:"
file $OUTDIR/*
```

```
## output files:
## /home/attila/projects/monoallelic-brain/data/predictors//predictors.csv:             ASCII text, with very long lines
## /home/attila/projects/monoallelic-brain/data/predictors//RNAseq_ID.DLPFC_RNA_ID.csv: ASCII text
```

## Importing into `R`

### Observed predictors and read counts

We will generate a data frame `E` representing matrix $E$ such that each column holds the observed values of some predictor (explanatory variable).  During the fitting of some regression model `R` will automatically construct the design matrix $X$ from $E$ (or a submatrix of $E$ if we omit some predictors) by making contrasts for categorical predictors (i.e. factors).

To get as set of responses for regression, we will also generate a list `Y` of data frames representing the set of matrices $Y = \{ Y_g : g\in\mathcal{G} \}$ so that each matrix $Y_g$ corresponds to gene $g$ in a selected set $\mathcal{G}$.  The columns of $Y_g$ are either the observed values $H_g$ (or $L_g$) of the "higher" (or "lower") read counts for $g$, respectively, or the observed values of derived statistics.  These derived statistics are the total read count $N_g = L_g + H_g$, the ratio $S_g = H_g / N_g$, whose rank-transformation yields $R_g$.  For all matrices $E, Y_{g_1}, Y_{g_2},...$ the rows correspond to $n$ observations, that is $n$ RNA samples, each from a different individual.

$\mathcal{G}$ is defined by the `gene.ids` below:

```r
             # 8 genes analyzed by Ifat
gene.ids <- c("PEG3", "INPP5F", "SNRPN", "PWAR6", "ZDBF2", "MEG3", "ZNF331", "GRB10",
             # 5 more genes analyzed by AGK 3/2/16
             "PEG10", "SNHG14", "NAP1L5", "KCNQ1OT1", "MEST",
             # 3 more genes present in data files
             "IGF2", "NLRP2", "UBE3A",
             # 'green' novel 1 MB imprinted genes; note that PWAR6 is already included above
             "TMEM261P1", "AL132709.5", "RP11-909M7.3", "SNORD116-20", "RP13-487P22.1", "hsa-mir-335", "PWRN1")
```

### Aggregations

We will redefine $Y$ after appending (horizontally concatenating) to it the set of matrices $\{Y_{a_k\mathcal{G}_l} : k=1,2 \;;\; l=1,2,... \}$, indexed by combinations of an aggregation method $a_l$ and a gene subset $\mathcal{G}_l\subset\mathcal{G}$.  $a_1$ is aggregation of $\{S_g : g\in\mathcal{G}_l\}$ (or the corresponding $\{R_g\}$) by weighted average (`WA`) using total read count weights $\{N_g\}$, whereas $a_2$ is unweighted average (`UA`).  Ifat's original analysis combined `UA` with the subset given by `gene.ids[1:8]`.  Thus, for observation $i$ the we have $\bar{S}_{i;k\mathcal{G}_j} = \left( \sum_{g\in\mathcal{G}_l} W_{ig} \right)^{-1} \sum_{g\in\mathcal{G}_l} W_{ig} S_{ig}$, where $W_{ig}=N_{ig}$ if $k=1$ (`WA`) and $W_{ig}=1$ if $k=2$ (`UA`).  Rank transformation of $\bar{S}_{i;k\mathcal{G}_j}$ yields $\bar{R}_{i;k\mathcal{G}_j}$.

### Filtering

Given threshold $t$ (`count.thrs`) and gene $g$, an observation $i$ is filtered out whenever $N_{ig}\le t$.  For weighted average over a subset $\mathcal{G}_l$ the filtering rule is $\sum_{g\in\mathcal{G}_l} N_{ig}\le t$; consequently for observation $i$ the statistic $S_{ig}$ may be filtered out from some individual genes $g\in\mathcal{G}_l$ but may still contribute to the weighted average $\bar{S}_{i;k\mathcal{G}_j}$.  This is not the case for unweighted averages, which are derived from gene-wise statistics only *after* those have been filtered.

### Results

#### New implementation

Using the new implementation...

```r
source("~/projects/monoallelic-brain/src/import-data.R")
```


```r
# default arguments given explictely to both function calls
E <- get.predictors(f.predictors = "~/projects/monoallelic-brain/data/predictors/predictors.csv",
                      f.rna.ids = "~/projects/monoallelic-brain/data/predictors/RNAseq_ID.DLPFC_RNA_ID.csv")
Y <- get.readcounts(gene.ids = gene.ids,
                           data.dir = "~/projects/monoallelic-brain/data/readcount/genes",
                           count.thrs = 50,
                           sel.obs = row.names(get.predictors()),
                           g.subsets = list(A.8 = gene.ids[1:8], A = gene.ids))
```

```
## Warning in max(y, na.rm = TRUE): no non-missing arguments to max; returning
## -Inf
```

The warning (if any) about the arguments to `max` may arise when all observations are filtered out for some gene like for TMEM261P1:

```r
sapply(Y, function(y) sum(! is.na(y)))
```

```
##          PEG3        INPP5F         SNRPN         PWAR6         ZDBF2 
##          2460          1980          1580          1930          1930 
##          MEG3        ZNF331         GRB10         PEG10        SNHG14 
##          2320          1285           970          1845          2375 
##        NAP1L5      KCNQ1OT1          MEST          IGF2         NLRP2 
##           915           955          1185            70           140 
##         UBE3A     TMEM261P1    AL132709.5  RP11-909M7.3   SNORD116-20 
##            95             0           665            55           925 
## RP13-487P22.1   hsa-mir-335         PWRN1          WA.8            WA 
##            35            20            50          2890          2890 
##          UA.8            UA 
##          1156          1156
```

#### Agreement with previous implementation

Comparing with my previous implementation (which was shown to give results consistent with Ifat's)...

```r
source("../2016-04-22-glm-for-s-statistic/2016-04-22-glm-for-s-statistic.R")
source("../2016-04-22-glm-for-s-statistic/2016-04-22-glm-for-s-statistic-run.R")
```

The old (left arguments) and new (right arguments) implementation agree perfectly on $N_g$ and $\bar{R}_{i;k\mathcal{G}_j}$.  For instance:

```r
c(identical(d$N_MEST, Y$MEST$N),
identical(d$R_MEST, Y$MEST$R),
identical(d$R_avg8, Y$UA.8$R))
```

```
## [1] TRUE TRUE TRUE
```

But there are slight differences regarding $S_g$ because the new implementation calculates it afresh from $L_g$ and $H_g$ whereas the old implementation imported rounded numbers from Ifat's `pop_skew_3June15.txt` file loosing some precision.

```r
c(identical(d$S_MEST, Y$MEST$S), all.equal(d$S_MEST, Y$MEST$S))
```

```
## [1] FALSE  TRUE
```
Moreover:

```r
c(identical(d$S_avg8, Y$WA.8$S), all.equal(d$S_avg8, Y$WA.8$S))
```

```
## [1] "FALSE"                                
## [2] "Mean relative difference: 0.001112457"
```

```r
c(identical(d$S_avg8, Y$UA.8$S), all.equal(d$S_avg8, Y$UA.8$S))
```

```
## [1] "FALSE"                                
## [2] "Mean relative difference: 0.009172348"
```
which shows that the `S_avg8` statistic corresponds to the `WA.8` weighted average $\bar{S}_{i}$ more closely than to the unweighted `UA.8`.  Analyzing the code of the old implementation in `../2016-04-22-glm-for-s-statistic/2016-04-22-glm-for-s-statistic.R` confirms this.

Finally, note that the names of predictors have been simplified in the new implementation:

```r
str(E[ , 1:15])
```

```
## 'data.frame':	579 obs. of  15 variables:
##  $ Institution  : Factor w/ 3 levels "MSSM","Penn",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ Gender       : Factor w/ 2 levels "Female","Male": 2 2 1 1 2 1 1 2 1 2 ...
##  $ Age          : int  42 58 28 36 52 78 49 62 60 51 ...
##  $ PMI          : num  22.3 19.5 22.8 17.3 22.2 16 20 20.8 24 21.3 ...
##  $ Dx           : Factor w/ 3 levels "AFF","Control",..: 1 1 1 1 1 1 1 1 1 1 ...
##  $ RIN          : num  6.9 7 6.9 6.9 7.6 7.4 8 8.4 7.7 7.9 ...
##  $ RIN2         : num  47.6 49 47.6 47.6 57.8 ...
##  $ RNA_lib_batch: Factor w/ 9 levels "0","A","B","C",..: 6 6 5 6 2 7 8 7 6 3 ...
##  $ Ancestry_EV.1: num  0.0214 0.0213 0.0202 0.0213 0.0213 ...
##  $ Ancestry_EV.2: num  0.00459 0.03477 -0.00671 0.02264 -0.00656 ...
##  $ Ancestry_EV.3: num  -0.003252 0.002797 0.000894 0.003056 0.006738 ...
##  $ Ancestry_EV.4: num  0.0381 -0.0202 0.0461 -0.0182 0.041 ...
##  $ Ancestry_EV.5: num  0.000824 -0.00345 -0.005654 -0.007156 -0.013422 ...
##  $ SV1          : num  0.02559 -0.00105 -0.00124 0.02031 -0.02787 ...
##  $ SV2          : num  -0.00651 -0.04842 -0.00523 -0.01151 -0.00416 ...
```

[summary]: http://katahdin.mssm.edu/ifat/web/cm/get_pop_freq.pl
