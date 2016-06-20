# explanatory variables (a.k.a. predictors)
e.vars <- c("Age",
               "Institution",
               "Gender",
               "PMI",
               "Dx",
               "RIN", "RIN2",
               "RNA_lib_batch",
               "Ancestry_EV.1", "Ancestry_EV.2", "Ancestry_EV.3", "Ancestry_EV.4", "Ancestry_EV.5" )

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
do.fit <- function(y = Y[[1]]$S,
                   X = E,
                   e.v = e.vars,
                   thrs = 0, # setting to Inf tolerates all points
                   ...) {
    # check if number of observations is tolerable
    if(sum(! is.na(y)) < thrs) return(NULL)
    # append explanatory variable y to the 'data' X
    X$Y <- y
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
                        preds = e.vars,
                        min.obs = 0,
                        sel.models = NULL) {
    # list of fitter functions, one for each combination of a model and a transformation
    l.fitters <- list(
                      # unweighted normal linear model with R
    unlm.R = function(z) {
        do.fit(y = z$R, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = NULL)
    },
                      # weighted normal linear model with R
    wnlm.R = function(z) {
        do.fit(y = z$R, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = z$N)
    },
                      # unweighted normal linear model with S
    unlm.S = function(z) {
        do.fit(y = z$S, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = NULL)
    },
                      # weighted normal linear model with S
    wnlm.S = function(z) {
        do.fit(y = z$S, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = z$N)
    },
                      # logistic model with S
    logi.S = function(z) {
        do.fit(y = cbind(z$H, z$L), X = G, e.v = preds, thrs = min.obs * 2, family = binomial)
    },
                      # logistic model with affine transformed S
    logi2.S = function(z) {
        affine.transform.S <- function(y) {
            H2 <- as.integer((y$S * 2 - 1) * y$N)
            C <- cbind(H2[], y$N - H2[])
            C[ C < 0 & ! is.na(C) ] <- 0
            return(C)
        }
        do.fit(y = affine.transform.S(z), X = G, e.v = preds, thrs = min.obs * 2, family = binomial)
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

mean.rel.diff <- function(target, current, ...){
    x <- mean(abs(target - current), ...) # mean absolute difference
    return(x / mean(abs(target)))
}

coefs4plot <- function(l.models, coef = "Age", sort.on = "estimate") {
    #sapply(l.models, function(m) summary(m)$coefficients[coef, 1])
    nulls <- sapply(l.models, is.null)
    x <- sapply(c(estimate=1, SE=2, t.val=3, p.val=4),
           function(k) sapply(l.models[ ! nulls ],
                              function(m) summary(m)$coefficients[coef, k]))
    # confidence intervals (composed of lower and upper confidence limits)
    CI <- cbind(x[ , "estimate"] - x[ , "SE"], x[ , "estimate"] + x[ , "SE"])
    colnames(CI) <- c("CL.lo", "CL.up")
    cbind(x, CI)[ sort(x[ , sort.on], index.return = TRUE)$ix, ]
}

plot.betas <- function(coefs) {
    nr <- nrow(coefs) # number of rows
    par(mar = c(5, 7, 4, 2))
    plot.new()
    plot.window(ylim = c(1, nr),
                xlim = c(min(coefs[ , "estimate"], na.rm = TRUE), max(coefs[ , "estimate"], na.rm = TRUE)))
    axis(1)
    axis(2, labels = rownames(coefs), at = seq_len(nr), las = 1)
    points(coefs[ , "estimate" ], seq_len(nr),
         pch = 15)
    invisible(lapply(seq_len(nr),
           function(g) lines(coefs[g, c("CL.lo", "CL.up")], c(g, g))))
    abline(v = 0, lty = "dashed")
    grid(nx = NA, ny = NULL)
    title(xlab = expression(hat(beta)[age] %+-% s.e.))
}
