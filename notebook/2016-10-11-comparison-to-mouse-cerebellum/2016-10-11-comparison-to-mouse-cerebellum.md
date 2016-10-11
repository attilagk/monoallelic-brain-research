## Preparation



Load functions


```r
source("~/projects/monoallelic-brain/src/import-data.R")
source("~/projects/monoallelic-brain/src/fit-glms.R")
source("~/projects/monoallelic-brain/src/graphics.R")
```


```r
human.g.symbols <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(human.g.symbols) <- human.g.symbols
E <- get.predictors() # default arguments
Y <- get.readcounts(human.g.symbols = human.g.symbols, count.thrs = 0)
```

```
## Error in get.readcounts(human.g.symbols = human.g.symbols, count.thrs = 0): unused argument (human.g.symbols = human.g.symbols)
```

```r
# exclude unweighed aggregates UA.8 and UA from fitting
to.fit.ids <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
to.fit.ids <- grep("^WA.8$", to.fit.ids, value = TRUE, invert = TRUE)
M <- do.all.fits(Y[to.fit.ids], preds = e.vars, sel.models = c(wnlm.Q = "wnlm.Q", logi.S = "logi.S"))
```


```r
library(AnnotationDbi)
library("org.Hs.eg.db")
library("org.Mm.eg.db")
library("hom.Hs.inp.db")
```


```r
human.eg.ids <- unlist(as.list(org.Hs.egSYMBOL2EG[mappedkeys(org.Hs.egSYMBOL2EG)])[human.g.symbols])
human.g.symbols[! human.g.symbols %in% names(human.eg.ids)]
```

```
##      AL132709.5    RP11-909M7.3   RP13-487P22.1     hsa-mir-335 
##    "AL132709.5"  "RP11-909M7.3" "RP13-487P22.1"   "hsa-mir-335"
```

```r
mouse.all <- as.list(org.Mm.egSYMBOL2EG[mappedkeys(org.Mm.egSYMBOL2EG)])
names(mouse.all) <- toupper(names(mouse.all))
human.g.symbols[! human.g.symbols %in% names(human.eg.ids)]
```

```
##      AL132709.5    RP11-909M7.3   RP13-487P22.1     hsa-mir-335 
##    "AL132709.5"  "RP11-909M7.3" "RP13-487P22.1"   "hsa-mir-335"
```


```r
perez <- read.csv("../../data/elife-07860/elife-07860-supp1-v2.csv")
```


```r
# convert mouse gene names (i.e. symbols) to upper case...
# ...and check for non-matching human gene symbols
names(human.g.symbols[! human.g.symbols %in% toupper(as.character(perez$gene_name))])
```

```
##  [1] "TMEM261P1"     "SNHG14"        "AL132709.5"    "RP11-909M7.3" 
##  [5] "ZIM2"          "PWAR6"         "FAM50B"        "SNURF"        
##  [9] "KCNQ1OT1"      "SNORD116-20"   "RP13-487P22.1" "ZNF331"       
## [13] "hsa-mir-335"   "DIRAS3"        "PWRN1"         "NLRP2"
```
