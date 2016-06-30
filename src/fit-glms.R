# explanatory variables (a.k.a. predictors)
e.vars <- c("Age",
               "Institution",
               "Gender",
               "PMI",
               "Dx",
               "RIN", "RIN2",
               "RNA_batch",
               "Ancestry.1", "Ancestry.2", "Ancestry.3", "Ancestry.4", "Ancestry.5" )

# Construct a formula from data and names of predictors and fit some glm
#
# Parameters
# y: observed responses: a vector (nlm) or a 2-column matrix (logistic models)
# X: observed predictors(see 'get.predictors' in 'import-data.R')
# e.v: names of predictors; their ordering typically affects some results (e.g. ANOVA)
# thrs: the least number tolerable obervations in response 'y' to perform the fit
# ...: arguments to be passed to glm such as model 'family' (default nlm)
#
# Value
# an objects of class glm (or NULL if n.obs < thrs)
do.fit <- function(response = Y[[1]]$S,
                   X = E,
                   e.v = names(E)[1:13],
                   thrs = 0, # setting to Inf tolerates all points
                   ...) {
    # check if number of observations is tolerable
    if(sum(! is.na(response)) < thrs) return(NULL)
    # append response variable the 'data' X
    X$Y <- response
    # construct formula and fit some model on 'data' X
    fm <- as.formula(paste("Y", "~", paste0(e.v, collapse = " + ")))
    glm(fm, data = X, ...)
}

# Perform all fits under selected models to read counts for each gene and each aggregate
#
# Parameters
# Z: read counts used as responses: a list of data frames from 'get.readcounts' (see 'import-data.R')
# G: observed explanatory variables in a data frame from 'get.predictors' (see 'import-data.R')
# preds: names of explanatory variables in a character vector; their ordering typically affects some results (e.g. ANOVA)
# min.obs: the least number tolerable obervations in response 'y' to perform the fit
# sel.models: fit selected models only: a numeric or character vector for the list of models l.fitters; NULL (default) selects all models
#
# Value
# a list of list of objects of class glm (or NULL if n.obs < min.obs); the
# outer list is named after models whereas the nested inner lists after genes
# and aggregates
#
# Details: this is a wrapper around 'do.fit'
do.all.fits <- function(Z = Y,
                        G = E,
                        preds = names(E)[1:13],
                        min.obs = 0,
                        sel.models = NULL,
                        ...) {
    # list of fitter functions, one for each combination of a model and a transformation
    l.fitters <- list(
                      # unweighted normal linear model with R
    unlm.R = function(z) {
        do.fit(response = z$R, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = NULL, ...)
    },
                      # weighted normal linear model with R
    wnlm.R = function(z) {
        do.fit(response = z$R, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = z$N, ...)
    },
                      # unweighted normal linear model with S
    unlm.S = function(z) {
        do.fit(response = z$S, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = NULL, ...)
    },
                      # weighted normal linear model with S
    wnlm.S = function(z) {
        do.fit(response = z$S, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = z$N, ...)
    },
                      # logistic model with S
    logi.S = function(z) {
        do.fit(response = cbind(z$H, z$L), X = G, e.v = preds, thrs = min.obs * 2, family = binomial, ...)
    },
                      # logistic model with affine transformed S
    logi2.S = function(z) {
        affine.transform.S <- function(y) {
            H2 <- as.integer((y$S * 2 - 1) * y$N)
            C <- cbind(H2[], y$N - H2[])
            C[ C < 0 & ! is.na(C) ] <- 0
            return(C)
        }
        do.fit(response = affine.transform.S(z), X = G, e.v = preds, thrs = min.obs * 2, family = binomial)
    }
    )
    # select fitters that correspond to selected models in 'sel.models'
    if(! is.null(sel.models))
        l.fitters <- l.fitters[ sel.models ]
    # Perform all fits:
    lapply(l.fitters, # Under each selected model...
           function(fitter) lapply(Z, # ...using the relevant predictor for each gene or aggregate...
                                   fitter)) # ...perform the fit!
}

# Calculate the mean difference of 'current' Y relative to 'target' X
#
# Parameters
# target: the reference vector
# current: the other vector
#
# Value
# a vector of mean relative differes
#
# Details
# mean difference of Y relative to X is avg(|X - Y|) / avg(|X|)
mean.rel.diff <- function(target, current, ...){
    x <- mean(abs(target - current), ...) # mean absolute difference
    return(x / mean(abs(target)))
}

# Get estimated coefficients and confidence intervals for a list of models
#
# Parameters
# l.models: list of models of class lm or glm
# coef.name: in case of multiple regression selects a coefficient
# conf.lev: confidence level
#
# Value
# A data frame with components (columns) lower.CI, upper.CI, beta.hat, and
# Gene.  Gene is a factor ordered by beta.hat, which is useful for plotting.
get.CI <- function(l.models, coef.name = "Age", conf.lev = 0.99, ...) {
    beta.hat <- sapply(l.models, function(m) coef(m)[coef.name])
    CI <- t(sapply(l.models, confint, parm = coef.name, level = conf.lev))
    df <- as.data.frame(cbind(CI, beta.hat))
    names(df)[1:2] <- c("lower.CI", "upper.CI")
    df$Gene <- reorder(factor(row.names(df)), df$beta.hat)
    return(df)
}

plot.CI <- function(df, package = "lattice", ...) {
    plotter <-
        list(lattice =
             function() {
                 segplot(Gene ~ lower.CI + upper.CI, data = df,
                         panel = function(x, y, ...) {
                             panel.grid(h = -1, v = -1)
                             panel.abline(v = 0, lty = "dashed", ...)
                             panel.segplot(x, y, ...)
                         },
                         draw.bands = FALSE, centers = beta.hat, xlab = expression(beta[age]),
                         main = expression(paste(hat(beta)[age], " and CI")))
             },
           ggplot2 =
               function() {
                   g <- ggplot(data = df, aes(x = Gene, y = beta.hat))
                   g <- g + coord_flip()
                   g <- g + geom_hline(yintercept = 0, linetype = 2)
                   g <- g + geom_pointrange(aes(ymin = lower.CI, ymax = upper.CI))
                   g <- g + labs(x = "", y = expression(beta[age]), title = expression(paste(hat(beta)[age], " and CI")))
                   g
               })
    plot(plotter[[ package ]]())
}

# creates a data frame of ANOVA deviances
#
# Parameters
# l.m: a list of models
#
# Value
# a data frame of ANOVA deviances where each column is a term and rows are genes
l.anova <- function(l.m) {
    y <- sapply(l.m, function(m) anova(m)[ , "Deviance" ])[ -1, ]
    row.names(y) <- attributes(terms(l.m[[1]]))$term.labels
    data.frame(t(y))
}

# creates a data frame of effects
#
# Parameters
# l.m: a list of models
# ref.m: a reference model with a complete set of coefficients
#
# Value
# a data frame of ANOVA deviances where each column is a term and rows are genes
#
# Details
# ref.m is needed to check if every other model has the same number of parameters
l.effects <- function(l.m, ref.m = 1) {
    ix <- sapply(l.m, function(m) length(coef(m)) == length(coef(l.m[[ ref.m ]])))
    y <- sapply(l.m[ ix ], function(m) effects(m)[ names(coef(l.m[[ref.m]]))[-1] ])
    data.frame(t(y))
}
