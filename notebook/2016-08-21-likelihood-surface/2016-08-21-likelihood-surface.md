## Data preparation


```
## Loading required package: RColorBrewer
```

Load data importer functions:

```r
source("../../src/import-data.R")
source("../../src/fit-glms.R")
source("../../src/likelihood-surface.R")
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

### Relative LL surface

Select a gene


```r
gene <- "PEG3"
```

The expressions below restrict the likelihood function, defined on a 22 dimensional parameter space, to various 2D subspaces (planar slices) defined by pairs of parameters; the remaining 20 parameters are held fixed at the maximum likelihood location $\cap{\beta}$.


```r
dat <- list()
l.M <- M$wnlm.Q; ll.fun <- ll.wnlm; args.fun <- args.wnlm
dat$coefs.wnlm.Q <- 
    rbind(ll.grid(l.M = l.M, v.name.A = "InstitutionPenn", CI.lev.A = c.a <- 1 - 1e-9, CI.lev.B = c.b <- 0.9, ll.fun = ll.fun, args.fun = args.fun),
          ll.grid(l.M = l.M, gene = gene, v.name.A = "PMI", CI.lev.A = 0.9999, CI.lev.B = c.b, ll.fun = ll.fun, args.fun = args.fun),
          ll.grid(l.M = l.M, gene = gene, v.name.A = "RIN", CI.lev.A = 0.5, CI.lev.B = c.b, ll.fun = ll.fun, args.fun = args.fun),
          ll.grid(l.M = l.M, gene = gene, v.name.A = "DxSCZ", CI.lev.A = 1 - 1e-6, CI.lev.B = c.b, ll.fun = ll.fun, args.fun = args.fun),
          ll.grid(l.M = l.M, gene = gene, v.name.A = "GenderMale", CI.lev.A = 1 - 1e-5, CI.lev.B = c.b, ll.fun = ll.fun, args.fun = args.fun),
          ll.grid(l.M = l.M, gene = gene, v.name.A = "Ancestry.2", CI.lev.A = 1 - 1e-11, CI.lev.B = c.b, ll.fun = ll.fun, args.fun = args.fun))
fac <- factor(dat$coefs.wnlm.Q$v.name.A)
dat$coefs.wnlm.Q$v.name.A <-
    factor(dat$coefs.wnlm.Q$v.name.A, ordered = TRUE,
           levels = c("InstitutionPenn", "PMI", "RIN", "DxSCZ", "GenderMale", "Ancestry.2"))
```



This is one example for a restricted likelihood function, defined on the subspace of $\{\beta_\mathrm{Age}\} \times \{\beta_\mathrm{Ancestry.2}\}$:


```r
ll.wireframe(par.settings = list(axis.line = list(col = "transparent")))
```

<img src="figure/explain-rll-wireframe-1.png" title="plot of chunk explain-rll-wireframe" alt="plot of chunk explain-rll-wireframe" width="600px" />

The same restricted function in a different representation:


```r
lp2 <- ll.surfaceplot(fm = formula(rel.log.L ~ beta.A * beta.B | v.name.A),
                      df = dat$coefs.wnlm.Q[dat$coefs.wnlm.Q$v.name.A == "Ancestry.2", ])
update(lp2, strip = FALSE, xlab = expression(paste(beta, "[ Ancestry.2 ]")), ylab = expression(paste(beta, "[ Age ]")), main = "rel. log likelihood")
```

<img src="figure/explain-rll-levelplot-1.png" title="plot of chunk explain-rll-levelplot" alt="plot of chunk explain-rll-levelplot" width="500px" />

### Association of $\beta_\mathrm{Age}$ to other coefficients

Relative log likelihood is shown for straight forward comparison between different models.  $\beta_\mathrm{Age}$ is associated with coefficients of both technical and biological predictors.  Association to the coefficient of RNA-quality RIN (RNA Integrity Number) is the strongest.  From these pairwise associations the network of causal dependencies does not follow directly.

#### wnlm.Q, PEG3


```r
ll.surfaceplot(fm = formula(rel.log.L ~ beta.A * beta.B | v.name.A), df = dat$coefs.wnlm.Q,
               layout = c(3, 2),
               par.settings = list(regions = list(col = rev(trellis.par.get("regions")$col))),
               main = "rel. log likelihood",
               xlab = expression(paste(beta, "[ ... ]")),
               ylab = expression(paste(beta, "[ Age ]")))
```

<img src="figure/ll-surf-coefs-wnlm-Q-1.png" title="plot of chunk ll-surf-coefs-wnlm-Q" alt="plot of chunk ll-surf-coefs-wnlm-Q" width="700px" />

#### logi.S, PEG3

<img src="figure/ll-surf-coefs-logi-S-1.png" title="plot of chunk ll-surf-coefs-logi-S" alt="plot of chunk ll-surf-coefs-logi-S" width="700px" />

### Consistency among genes

The same pattern of associations can be observed for other genes as well.  The model for PEG3 is supported by more observations and total read count than for UBE3A or ZNF331, and therefore the likelihood surface is more peaked (higher information content, larger curvature around the estimate).


```r
l.M <- M$wnlm.Q; ll.fun <- ll.wnlm; args.fun <- args.wnlm
dat$genes.wnlm.Q <- # using v.name.A = "InstitutionPenn" as default
    rbind(ll.grid(l.M = l.M, gene = "PEG3", CI.lev.A = 1-1e-3, CI.lev.B = 0.7, ll.fun = ll.fun, args.fun = args.fun),
          ll.grid(l.M = l.M, gene = "UBE3A", CI.lev.A = 1-1e-15, CI.lev.B = 1-1e-2, ll.fun = ll.fun, args.fun = args.fun),
          ll.grid(l.M = l.M, gene = "ZNF331", CI.lev.A = 1-1e-15, CI.lev.B = 1-2e-2, ll.fun = ll.fun, args.fun = args.fun))
```

<img src="figure/ll-surf-genes-wnlm-Q-1.png" title="plot of chunk ll-surf-genes-wnlm-Q" alt="plot of chunk ll-surf-genes-wnlm-Q" width="700px" />

