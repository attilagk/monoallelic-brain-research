## Data preparation



Load data importer functions:

```r
source("../../src/import-data.R")
source("../../src/fit-glms.R")
source("2016-08-23-glm-sampling-distributions.R")
```

Do the import and fit all simple and multiple regression models as before (code hidden).





## Results

The plots below show the first 4 of the 6 combinations (logi.S, logi2.S, wnlm.S, wnlm.R, unlm.S, unlm.R) of some model family and some response variable ($S$ or $R$) that were used to fit data.  unlm.S and unlm.R are not shown because they are so similar to wnlm.S and wnlm.R, respectively.

For all predicted points of logistic models the total read count (used as the denominator of the binomial error function) was set to be the average of the observed total read counts (given the gene in question).

Theoretically possible values (non-negative age and S between 0.5 and 1 or S between 0 and 100) are demarcated by a dotted rectangle.  The logi2.S model respects both the upper and lower side of the rectangle but on the lower side it has asymptotically zero variance, which contradicts with the observed large variance of well-balanced (biallelically expressed) genes.

Unlike the normal linear models, both logistic models explain the systematic change of variance with average read count.  However, the explained change seems somewhat less than what is observed indicating overdispersion.  The weakness of logistic models is that the steepest, most informative, part of their sigmoidal predicted curve lies out of the observed range of data (massive extrapolation).

The rank transformation achieves nearly constant variance allowing a better fit of the normal liear model.  But much information is lost in the transformation, which substantially worsens signal-to-noise (explained vs unexplained variance).


```r
gp <- 
    lapply(c(GRB10 = "GRB10", ZNF331 = "ZNF331", KCNK9 = "KCNK9", PEG3 = "PEG3", PEG10 = "PEG10"),
             grid.predictions.1gene)
```

### Examples for down-regulation

<img src="figure/GRB10-1.png" title="plot of chunk GRB10" alt="plot of chunk GRB10" width="700px" />

<img src="figure/PEG3-1.png" title="plot of chunk PEG3" alt="plot of chunk PEG3" width="700px" />

<img src="figure/ZNF331-1.png" title="plot of chunk ZNF331" alt="plot of chunk ZNF331" width="700px" />

### Examples for up-regulation

<img src="figure/KCNK9-1.png" title="plot of chunk KCNK9" alt="plot of chunk KCNK9" width="700px" />

<img src="figure/PEG10-1.png" title="plot of chunk PEG10" alt="plot of chunk PEG10" width="700px" />

### Figure for manuscript

The figures illustrate the widening of prediction intervals under the logistic models but not the normal linear models.  The discrete nature of the binomial error distribution under logistic models can also be seen: the denominator of the binomial distribution is 61 (GRB10) and 16  (KCNK9); these were defined as the average of the total read count for a given gene.


```r
print(plots$GRB10, split = c(1, 1, 1, 2), more = TRUE)
print(plots$KCNK9, split = c(1, 2, 1, 2), more = FALSE)
```

<img src="figure/density-and-predicted-curve-1.png" title="plot of chunk density-and-predicted-curve" alt="plot of chunk density-and-predicted-curve" width="700px" height="1400" />
