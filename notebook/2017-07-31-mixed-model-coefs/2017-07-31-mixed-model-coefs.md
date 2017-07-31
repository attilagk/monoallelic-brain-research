---
layout: default
tags: [regression, anova ]
---

Estimate regression coefficients of a mixed effect model.  This requires extracting the conditional modes of random effects.



Selected genes (inferred to be imprinted)


```r
gene.ids <- unlist(read.csv("../../data/genes.regression.new", as.is = TRUE))
names(gene.ids) <- gene.ids
dat <- merge.data(gene.ids = gene.ids)
```

Get model $$M5$$ formula and fit to data


```r
get.formula <- function(model.name = "M5") {
    x <- read.csv(file = "../../results/M-formulas.csv", stringsAsFactors = FALSE)[[model.name]]
    formula(do.call(paste, as.list(x[c(2, 1, 3)])))
}
M5 <- lmer(get.formula("M5"), data = dat)
M6 <- lmer(get.formula("M6"), data = dat)
```


