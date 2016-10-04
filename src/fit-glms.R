# explanatory variables (a.k.a. predictors)
e.vars <- c("Age",
               "Institution",
               "Gender",
               "PMI",
               "Dx",
               "RIN",
               #"RIN2", # this was removed after much analysis
               "RNA_batch",
               "Ancestry.1", "Ancestry.2", "Ancestry.3", "Ancestry.4", "Ancestry.5" )

# Construct a formula from data and names of predictors and fit some glm
#
# Arguments
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
                   e.v = names(E)[1:12],
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
# Arguments
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
                        preds = names(E)[1:12],
                        min.obs = 0,
                        sel.models = NULL,
                        ...) {
    # list of fitter functions, one for each combination of a model and a transformation
    l.fitters <- list(
                      # unweighted normal linear model with Q
    unlm.Q = function(z) {
        do.fit(response = z$Q, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = NULL, ...)
    },
                      # weighted normal linear model with Q
    wnlm.Q = function(z) {
        do.fit(response = z$Q, X = G, e.v = preds, thrs = min.obs, family = gaussian, weights = z$N, ...)
    },
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
# Arguments
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
# Arguments
# l.models: list of models of class lm or glm
# coef.name: in case of multiple regression selects a coefficient
# conf.lev: confidence level
#
# Value
# A data frame with components (columns) Lower.CL, Upper.CL, beta.hat, and
# Gene.  Gene is a factor ordered by beta.hat, which is useful for plotting.
get.CI <- function(l.models,
                   all.coefs = names(coef(l.models[[ order(sapply(lapply(l.models, coef), length), decreasing=TRUE)[1] ]])),
                   coef.name = "Age",
                   conf.lev = 0.99,
                   ...) {
    beta.hat <- sapply(l.models, function(m) coef(m)[coef.name])
    CI <- t(sapply(l.models, confint, parm = coef.name, level = conf.lev))
    df <- as.data.frame(cbind(CI, beta.hat))
    names(df)[1:2] <- c("Lower.CL", "Upper.CL")
    cbind(df,
          list(Gene = factor(row.names(df), levels = row.names(df), ordered = TRUE)),
          list(Coefficient = factor(coef.name, levels = all.coefs, ordered = TRUE)))
}

# low level function extracting estimate and CI for a single model object
#
# Arguments
# m: the model object
# conf.lev: confidence level
# skip.CI: do NOT calculate CI (saves much time)
#
# Value
# A data frame containing the estimate, the lower and upper CL and some other
# info.
get.single.estimate.CI <-
    function(m, conf.lev, skip.CI = FALSE) {
        d <- cbind(data.frame(Estimate = cf <- coef(m)),
                   data.frame(Coefficient =
                              factor(n.cf <- names(cf), levels = n.cf, ordered = TRUE)))
        if(! skip.CI) {
            d <- cbind(d, data.frame(tryCatch(confint(m, level = conf.lev),
                                              error = function(e) matrix(NA, nrow = 1, ncol = 2))))
            names(d)[3:4] <- c("Lower.CL", "Upper.CL")
        }
        return(d)
    }

# extract estimates and confidence intervals for regression coefficients
#
# Arguments
# l.M: a list of models
# conf.lev: confidence level
#
# Value
# a data frame in a long form for easy plotting
#
# Details
# Data are set to NA for models that haven't converged.
get.estimate.CI <- function(l.M, conf.lev = 0.99) {
    df <- do.call(make.groups,
                  lapply(l.M, get.single.estimate.CI, conf.lev))
    names(df)[5] <- c("Gene")
    not.converged <- df$Gene %in% names(l.M)[ ! sapply(l.M, `[[`, "converged") ]
    df[ not.converged, c("Estimate", "Lower.CL", "Upper.CL") ] <- NA
    return(df)
}


# get all coefficients for a given predictor
#
# Arguments
# m: a model object
# e.v: the name of the predictor
#
# Value:
# One or more coefficient names that are associated with the predictor.  For
# covariates (i.e. continuous variables) it is simply the name of the
# predictor without any modification.  For factors it is the name of one or
# more levels after the reference level (in the sense of contrast) was
# removed.
predictor2coefs <- function(m, e.v) {
    if(! e.v %in% names(m$xlevels)) return(e.v) # covariates (1 parameter for predictor e.v)
    paste0(e.v, m$xlevels[[e.v]][-1])
}


# aggregate estimates and CI across multiple permutations, genes, permuted variables, and models
#
# Arguments
# perms: a data frame of integer vectors that are random permutations of the original case/individual numbers
# gene.ids: the id (symbol) for the selected genes
# e.vars: names of all predictors (explanatory variables)
# sel.vars: names of selected predictors based on which permutation is done
# sel.models: selected model types, e.g. wnlm.R, logi.S,...
# E: data frame of predictor data
# Y: read count data
# conf.lev: confidence level for CI
# skip.CI: do NOT calculate CI (saves much time)
#
# Value
# A thin long data frame whose columns are not only the estimate and the lower
# and upper CL but also annotations such as the genes, model type, etc.
aggregate.CI.permut2 <- function(perms, gene.ids, e.vars, sel.vars = e.vars,
                                 sel.models = list(wnlm.R = "wnlm.R", logi.S = "logi.S"),
                                 E, Y, conf.lev = 0.99, skip.CI = FALSE) {
    helper <- function(perm.name, gene, e.v, type) {
        # permute cases for predictor e.v
        perm <- perms[[perm.name]]
        X <- E
        X[[e.v]] <- E[[e.v]][perm]
        # fit model on X and the given gene's read count data
        m <- do.all.fits(Z = Y[gene], G = X, preds = e.vars, sel.models = type)[[type]][[gene]]
        # a data frame for the estimate and CI
        d <- get.single.estimate.CI(m, conf.lev = conf.lev, skip.CI = skip.CI)[predictor2coefs(m, e.v), ]
        # annotate conditions
        d$Permutation <- factor(perm.name, levels = names(perms), ordered = TRUE)
        d$Gene <- factor(gene, levels = gene.ids, ordered = TRUE)
        d$Perm.Var <- factor(e.v, levels = e.vars, ordered = TRUE)
        d$Model <- factor(type, levels = sel.models, ordered = TRUE)
        # replace estimate and CI with NA if fit has not converged
        d[ ! m$converged, "Estimate" ] <- NA
        if(! skip.CI)
            d[ ! m$converged, c("Lower.CL", "Upper.CL") ] <- NA
        return(d)
    }
    # aggregation...
    do.call(rbind,
            lapply(sel.models, # ...across models
                   function(type)
                       do.call(rbind,
                               lapply(sel.vars, # ...across permuted variables
                                      function(e.v)
                                          do.call(rbind,
                                                  lapply(gene.ids, # ...across genes
                                                         function(gene)
                                                             do.call(rbind,
                                                                     lapply(names(perms), # ...across permutations
                                                                            helper, gene, e.v, type))))))))
}

# calculate p-value for H_0: param beta = 0 based on permutations
#
# Arguments
# beta.U: the unpermuted (i.e. observed) estimate for beta
# betas: vector of permuted estimates
# two.tailed: TRUE for H_0: beta = 0; FALSE for H_0: beta >(or <) 0 given that beta.U >(or <) 0
# do.norm: whether normalize the number of permuted betas as extreme as beta.U (TRUE for p-value)
#
# Value
# The estimated p-value for H_0 or, if do.norm = FALSE, the number of 
# estimated betas from permutations as extreme the estimate based on the
# observed data.
#
# Details
# The calculation uses H_0: beta = median(betas) instead of H_0 beta = 0 to
# avoid p-values > 1.  In the limit of the number of permutations it is
# expected that median(betas) -> 0.
get.p.val <- function(beta.U, betas, two.tailed = TRUE, do.norm = TRUE) {
    L.tail <- function() sum(beta.U >= betas, na.rm = TRUE)
    R.tail <- function() sum(beta.U <= betas, na.rm = TRUE)
    res <- if(beta.U > median(betas, na.rm = TRUE)) R.tail() else L.tail()
    res <- if(two.tailed) 2 * res
    if(do.norm) res / length(betas) else res
}

# an earlier, less complete aggregator function similar to aggregate.CI.permut2
#
# Arguments
# type: model type
# M: a list of list of list of models structured hierarchically by permuted predictors, model types and genes (from outer to inner lists)
# e.v: the name for all predictors
# conf.lev: the confidence level of the CI
#
# Value
# a long thin data frame
get.estimate.CI.permut <- function(type = "logi.S", M, e.v, conf.lev = 0.99) {
    betas <- lapply(e.v,
                    function(e) {
                        x <- get.estimate.CI(M[[e]][[type]], conf.lev = conf.lev)
                        x <- x[ grep(e, row.names(x), value = TRUE), ]
                    })
    betas$RIN2 <- NULL #  regex "RIN" already matched both RIN and RIN2
    return(do.call(rbind, betas))
}


# Earlier implementation of 'get.estimate.CI'.  Essentially the same behavior
# but somewhat less elegant program.
get.estimate.CI.old <- function(l.M, conf.lev = 0.99) {
    # get the names of all coefficients from one of the genes with a full set of coefs
    all.coefs <- names(coef(l.M[[ order(sapply(lapply(l.M, coef), length), decreasing=TRUE)[1] ]]))
    l.M.names <- names(l.M)
    helper <- function(n) {
        df <- data.frame(cbind(coef(l.M[[ n ]]), confint(l.M[[ n ]], level = conf.lev)))
        names(df) <- c("Estimate", "Lower.CL", "Upper.CL")
        cbind(df,
              list(Coefficient = factor(row.names(df), levels = all.coefs, ordered = TRUE)),
              list(Gene = factor(n, levels = l.M.names, ordered = TRUE)))
    }
    df <- Reduce(rbind, lapply(l.M.names, helper))
    not.converged <- df$Gene %in% names(l.M)[ ! sapply(l.M, `[[`, "converged") ]
    df[ not.converged, ] <- NA
    return(df)
}

# creates a long data frame of effects from a list of models
#
# Arguments
# l.m: a list of models
# coef.names: the name of all coefficients across all genes
#
# Value
# a long data frame with 3 columns: Effect, Coefficient, Gene. The latter two are ordered factors
l.effects <- function(l.m, coef.names = names(coef(l.m[[1]]))) {
    df <- do.call(make.groups,
                  lapply(l.m,
                         function(m) {
                             ef <- effects(m)
                             data.frame(Effect = ef <- ef[ names(ef) != "" ],
                                              Coefficient = factor(names(ef), levels = coef.names, ordered = TRUE))
                         }))
    names(df)[3] <- "Gene"
    return(df)
}

# Earlier implementation of 'l.effects' with slightly different behavior.  The
# difference concerns those genes for which some of the coefficients lack the
# corresponding component in the vector produced by calling 'effects'.  In
# this old implementation these genes are removed but in the new 'l.effects'
# they are retained.
l.effects.old <- function(l.m, ref.m = 1) {
    ix <- sapply(l.m, function(m) length(coef(m)) == length(coef(l.m[[ ref.m ]])))
    y <- sapply(l.m[ ix ], function(m) effects(m)[ names(coef(l.m[[ref.m]]))[-1] ])
    data.frame(t(y))
}

# creates a long data frame of effects from a list of lists of models
#
# Arguments
# l.l.m: a list of lists of models (see Details)
# coef.names: the name of all coefficients across all genes
#
# Value
# A long data frame with 4 columns: Effect, Coefficient, Gene, Order. The
# latter three are ordered factors.
#
# Details
# The structure of 'l.l.m' shows the following pattern: the main (outer) list
# corresponds to the *order* of predictors whereas the inner lists correspond
# to the set of genes/gene aggregates.  The latter is obtained with
# 'do.all.fits'.
l.l.effects <- function(l.l.m, coef.names = names(coef(l.l.m[[1]][[1]]))) {
    df <- do.call(make.groups, lapply(l.l.m, l.effects, coef.names = coef.names))
    names(df)[4] <- "Order"
    return(df)
}

# creates a data frame of ANOVA deviances
#
# Arguments
# l.m: a list of models
#
# Value
# a data frame of ANOVA deviances where each column is a term and rows are genes
l.anova <- function(l.m) {
    y <- sapply(l.m, function(m) anova(m)[ , "Deviance" ])[ -1, ]
    row.names(y) <- attributes(terms(l.m[[1]]))$term.labels
    data.frame(t(y))
}

# reshape a data frame of ANOVA/Effects info to long format for easy plotting
#
# Arguments
# A: a "wide" data frame of deviances/effects, whose rows are genes and columns predictors/coefficients
# type: type of informatino either "anova" for ANOVA deviances or "effects" for effects
#
# Value
# 'A' reshaped into a long format

# Details
# In the reshaped 'A' all deviances/effects are accumulated into a single
# component (i.e. column), and two new components are added: (1) Gene and (2),
# depending on 'type', Predictor/Coefficient; both of these components are
# ordered factors.
reshape.1 <- function(A, type = "anova") {
    types <- list(anova = c("Deviance", "Predictor"), effects = c("Effect", "Coefficient"))
    reshape(A,
            v.names = types[[type]][1], varying = list(names(A)),
            timevar = types[[type]][2], times = factor(names(A), levels = names(A), ordered = TRUE),
            idvar = "Gene", ids = factor(row.names(A), levels = row.names(A), ordered = TRUE),
            direction = "long")
}

# reshape a *list* of data frames of ANOVA/Effects info to long format for easy plotting
#
# Arguments
# l.A: a list of data frames in "wide" format; see parameter 'A' in reshape.1
# type: either "anova" or "effects"; see parameter 'A' in reshape.1
#
# Value
# 'l.A' reshaped into a *single* data frame in a long format
#
# Details
# Each data frame in 'l.A' is assumed to correspond to some *order* in which
# predictors are added successively in ANOVA.  Like in the case of 'reshape.1'
# all deviances/effects are accumulated into a single
# component (i.e. column), and a Gene as well as a Predictor/Coefficient
# component is added.  But unlike in 'reshape.1' and 'Order' component (also
# an ordered factor) is also created.
reshape.2 <- function(l.A, type = "anova") {
    Reduce(rbind, lapply(names(l.A), function(n)
                         cbind(reshape.1(l.A[[n]], type = type),
                                   list(Order = factor(n, levels = names(l.A), ordered = TRUE)))))
}

# Make index for a subset of observations
#
# Arguments
# X: a data frame of observed values of predictors
# predictor: a component (column) in 'X' interpreted as an explanatory variable
# lvls: the character vector of one or more levels of 'predictor' to be selected
#
# Value
# a character vector to be used as row name indeces for the subset to be extracted
sset.obs <- function(predictor, lvls, X = E) {
    row.names(X)[X[[predictor]] %in% lvls]
}

# Fit a single model to the same subset for each dataset
#
# Arguments
# sel.pred: a selected predictor (character)
# levs: one or more levels of 'sel.pred' to be used for subsetting
# preds, Z, G: see the corresponding parameters of 'do.all.fits'
# s.mod: a *single* selected model (unlike in 'do.all.fits', where 'sel.models' may name multiple models)
#
# Value
# a list of model (l1m) objects (see 'do.all.fits'), where each component
# corresponds to a specific dataset (i.e. a gene)
do.all.fits.sset <- function(sel.pred, levs,
                             Z = Y, G = E,
                             preds = e.vars[! e.vars %in% sel.pred],
                             s.mod = "logi2.S") {
    ss <- sset.obs(sel.pred, levs, G)
    Z <- lapply(Z, `[`, ss, TRUE)
    G <- G[ss, preds]
    do.all.fits(Z, G, preds = preds, sel.models = s.mod)[[s.mod]]
}

# Fit a single model to a *l*ist of subsets for each dataset given a single predictor
#
# Arguments
# see 'do.all.fits.sset' and 'do.all.fits'
#
# Value
# a list of lists of models (l2m) where the outer list corresponds to subsets
# (i.e. levels of the selected predictor) and the inner list to genes (see
# 'do.all.fits.sset')
l.do.all.fits.sset <- function(sel.pred = "Institution", Z = Y, G = E,
                               preds = e.vars[! e.vars %in% sel.pred],
                               s.mod = "logi2.S") {
    l <- levels(G[[sel.pred]])
    names(l) <- l
    lapply(l, function(lev)
           do.all.fits.sset(sel.pred, lev, Z, G, preds, s.mod))
}

# Fit a single model to a *l*ist of subsets for each dataset given *multiple* predictor
#
# Arguments
# see 'do.all.fits.sset' and 'do.all.fits'
#
# Value
# a list of lists of lists of models (l3m), where the outermost list
# corresponds to selected predictors and the lower level nestings (l2m, l1m)
# have structures that are described in 'l.do.all.fits.sset' and
# 'do.all.fits.sset'
l.l.do.all.fits.sset <- function(sel.preds = c("Institution", "Gender"),
                                 Z = Y, G = E, e.v = e.vars, s.mod = "logi.S") {
    names(sel.preds) <- sel.preds
    M <- lapply(sel.preds,
                function(s.p)
                    l.do.all.fits.sset(s.p, Z, G, preds = e.v[! e.v %in% s.p], s.mod))
    M$Marginal <-
        list(All = do.all.fits(Z = Z, G = G, preds = e.v, min.obs = 0, sel.models = s.mod)[[s.mod]])
    return(M)
}
