## Data preparation



Load data importer functions:

```r
source("../../src/import-data.R")
source("../../src/fit-glms.R")
source("2016-08-23-glm-sampling-distributions.R")
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
M <- list()
f.ids <- list()
M$multi <- do.all.fits(Y[ids2fit], preds = e.vars)
f.ids$multi <- as.data.frame(lapply(M$multi, function(m) ! sapply(m, is.null)))
f.ids$multi["TMEM261P1", c("logi.S", "logi2.S")] <- FALSE
f.ids$multi[c("WA.8", "WA"), ] <- FALSE # gene aggregates are not needed here
M$simple <- do.all.fits(Y[ids2fit], preds = "Age")
f.ids$simple <- as.data.frame(lapply(M$simple, function(m) ! sapply(m, is.null)))
f.ids$simple["TMEM261P1", c("logi.S", "logi2.S")] <- FALSE
f.ids$simple[c("WA.8", "WA"), ] <- FALSE # gene aggregates are not needed here
```


```r
plot.mdls.1gene <- function(gene = "GRB10", l.l.M = M$simple) {
    lp <- list()
    lp$logi.S <-
        update(plot.density(l.l.M$logi.S[[gene]], type = "logi", ylim = c(0.4, 1.1), is.R = FALSE),
               main = "logi.S")
    lp$logi.S.large <-
        update(plot.density(l.l.M$logi.S[[gene]], type = "logi", xlim = c(-500, 1500), ylim = c(-0.1, 1.1), is.R = FALSE),
               main = "logi.S")
    lp$wnlm.S <-
        update(plot.density(l.l.M$wnlm.S[[gene]], type = "nlm", ylim = c(0.4, 1.1), is.R = FALSE),
               main = "wnlm.S")
    lp$wnlm.R <-
        update(plot.density(l.l.M$wnlm.R[[gene]], type = "nlm", ylim = c(-9, 110), is.R = TRUE),
               main = "wnlm.R")
    lp$unlm.S <-
        update(plot.density(l.l.M$unlm.S[[gene]], type = "nlm", ylim = c(0.4, 1.1), is.R = FALSE),
               main = "unlm.S")
    lp$unlm.R <-
        update(plot.density(l.l.M$unlm.R[[gene]], type = "nlm", ylim = c(-9, 110), is.R = TRUE),
               main = "unlm.R")
    return(lp)
}
```


```r
plots <- list()
plots$GRB10 <- plot.mdls.1gene(gene = "GRB10", l.l.M = M$simple) 
print(plots$GRB10$logi.S, split = c(1, 1, 2, 2), more = TRUE)
print(plots$GRB10$logi.S.large, split = c(2, 1, 2, 2), more = TRUE)
print(plots$GRB10$wnlm.S, split = c(1, 2, 2, 2), more = TRUE)
print(plots$GRB10$wnlm.R, split = c(2, 2, 2, 2), more = FALSE)
```

<img src="figure/GRB10-1.png" title="plot of chunk GRB10" alt="plot of chunk GRB10" width="700px" />


```r
plots$PEG3 <- plot.mdls.1gene(gene = "PEG3", l.l.M = M$simple) 
print(plots$PEG3$logi.S, split = c(1, 1, 2, 2), more = TRUE)
print(plots$PEG3$logi.S.large, split = c(2, 1, 2, 2), more = TRUE)
print(plots$PEG3$wnlm.S, split = c(1, 2, 2, 2), more = TRUE)
print(plots$PEG3$wnlm.R, split = c(2, 2, 2, 2), more = FALSE)
```

<img src="figure/PEG3-1.png" title="plot of chunk PEG3" alt="plot of chunk PEG3" width="700px" />
