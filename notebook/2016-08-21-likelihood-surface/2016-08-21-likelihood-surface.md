## Data preparation



Load data importer functions:

```r
source("../../src/import-data.R")
source("../../src/fit-glms.R")
source("2016-08-21-likelihood-surface.R")
```

Do the import:


```r
genes2fit <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(genes2fit) <- genes2fit
E <- get.predictors() # default arguments
Y <- get.readcounts(gene.ids = genes2fit, count.thrs = 0)
```

Fit models:


```r
# exclude unweighed aggregates UA.8 and UA from fitting
ids2fit <- grep("^UA(.8)?$", names(Y), value = TRUE, invert = TRUE)
M <- do.all.fits(Y[ids2fit], preds = e.vars)
f.ids <- as.data.frame(lapply(M, function(m) ! sapply(m, is.null)))
f.ids["TMEM261P1", c("logi.S", "logi2.S")] <- FALSE
f.ids[c("WA.8", "WA"), ] <- FALSE # gene aggregates are not needed here
```

## Results

### Association of $\beta_\mathrm{Age}$ to other coefficients

Relative log likelihood is shown for straight forward comparison between different models.  $\beta_\mathrm{Age}$ is associated with coefficients of both technical and biological predictors.  Association to the coefficient of RNA-quality RIN (RNA Integrity Number) is the strongest.  From these pairwise associations the network of causal dependencies does not follow directly.


```r
dat <- list()
dat$coefs <- # using gene = "PEG3" as default
    rbind(ll.grid(v.name.A = "InstitutionPenn", CI.lev.A = c.a <- 0.999999, CI.lev.B = c.b <- 0.9),
          ll.grid(v.name.A = "PMI", CI.lev.A = 0.9999, CI.lev.B = c.b),
          ll.grid(v.name.A = "RIN", CI.lev.A = 0.05, CI.lev.B = c.b),
          #ll.grid(v.name.A = "RNA_batchB", CI.lev.A = 0.95, CI.lev.B = c.b),
          ll.grid(v.name.A = "DxControl", CI.lev.A = 0.99, CI.lev.B = c.b),
          ll.grid(v.name.A = "GenderMale", CI.lev.A = c.a, CI.lev.B = c.b),
          ll.grid(v.name.A = "Ancestry.2", CI.lev.A = 1 - 1e-10, CI.lev.B = c.b))
fac <- factor(dat$coefs$v.name.A)
dat$coefs$v.name.A <-
    factor(dat$coefs$v.name.A, ordered = TRUE,
           levels = c("InstitutionPenn", "PMI", "RIN", "DxControl", "GenderMale", "Ancestry.2"))
```


```r
ll.surfaceplot(fm = formula(rel.log.L ~ beta.A * beta.B | v.name.A), df = dat$coefs,
               layout = c(3, 2),
               par.settings = list(regions = list(col = rev(trellis.par.get("regions")$col))),
               main = "rel. log likelihood, PEG3",
               xlab = expression(paste(beta, "[ ... ]")),
               ylab = expression(paste(beta, "[ Age ]")))
```

<img src="figure/ll-surf-coefs-1.png" title="plot of chunk ll-surf-coefs" alt="plot of chunk ll-surf-coefs" width="700px" />

### Consistency among genes

The same pattern of associations can be observed for other genes as well.


```r
dat$genes <- # using v.name.A = "InstitutionPenn" as default
    rbind(ll.grid(gene = "PEG3", CI.lev.A = c.a <- 0.9999, CI.lev.B = c.b <- 0.9),
          ll.grid(gene = "ZNF331", CI.lev.A = c.a, CI.lev.B = c.b),
          ll.grid(gene = "UBE3A", CI.lev.A = c.a, CI.lev.B = c.b))
```

<img src="figure/ll-surf-genes-1.png" title="plot of chunk ll-surf-genes" alt="plot of chunk ll-surf-genes" width="700px" />

