## Preparing data for analysis in `R`

### Read counts from Ifat's website

The source of read count data is the largest HTML table generated dynamically using the `All genes` option in [Ifat's website][summary].  The output is a `csv` file for each gene within a gene-specific subdirectory.  Due to the specifics of Ifat's html tables and the general difficulty of the `html` to `csv` conversion only the "higher" and "lower" read counts could be formatted in a usable way; the underlying SNP-specific read counts could not be neatly retrieved just as the corresponding genotypes.


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

Given `predictors.csv` and `RNAseq_ID.DLPFC_RNA_ID.csv` the data frame `predictors` has the following structure:


```r
source("~/projects/monoallelic-brain/src/import-data.R")
str(predictors[ , 1:15])
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

Each observation corresponds to an RNA sample (and the corresponding individual)


```r
head(row.names(predictors))
```

```
## [1] "CCM_MSSM_bp_DLPFC_10" "CCM_MSSM_bp_DLPFC_11" "CCM_MSSM_bp_DLPFC_12"
## [4] "CCM_MSSM_bp_DLPFC_13" "CCM_MSSM_bp_DLPFC_14" "CCM_MSSM_bp_DLPFC_15"
```

[summary]: http://katahdin.mssm.edu/ifat/web/cm/get_pop_freq.pl
