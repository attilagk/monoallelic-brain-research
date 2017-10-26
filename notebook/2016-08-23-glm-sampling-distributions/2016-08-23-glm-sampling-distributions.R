# create a long data frame for model predictions and probability density
#
#
# Arguments
#
# lll.M: a list of list of list of model objects
# M.multi: simple or multiple regression
# M.family: one of logi.S, logi2.S, wnlm.S, wnlm.R, unlm.S, unlm.R
# gene: a single selected gene
# xlim: the x and y limits of the density and predicted values
# l.ylims: a list of ylims whose components correspond to different M.families
# n.pnts: the x-resolution (in number of ponts) and the minimum y-resolution (see details)
#
#
# Value: a thin, long data frame with the following components (columns)
#
# x: n.pnts evenly spaced values within xlim of predictor x
# y: minimum n.pnts evenly spaced values within ylim of redicted response Y
# mu: the predicted response (i.e. conditional mean of Y given x)
# multi: a factor with levels 'simple < multiple' indicating the number of predictors
# family: a factor with levels 'logi.S < logi2.S < wnlm.S < wnlm.R < unlm.S < unlm.R'
# gene: a factor whose levels correspond to all genes in lll.M
# density: the probability density of the model, scaled to allow comparison across models
# id: not important; created as a result of calling the reshape function
# obs.x, obs.y: the observed data (x: Age and y: read count ratio S or its rank-transformed version R)
#
#
# Details
#
# n.pnts determines the x-resolution for all model families and the
# y-resolution for all the normal linear families (i.e those matching
# [wu]nlm.[RS]).  As for the logistic model families (logi.S and logi2.S) the
# y-resolution is n.pnts only if n.pnts > n and n othewise, where n is the average total read
# count for the model and plays a role of the denominator of the binomial
# error distribution of logistic models.  Note that a *scalar* n means that the
# weight on each predicted/theoretical/simulated data point is the same.  But,
# of course, all models (except for unlm.[RS]) have been fitted with the observed
# *vector* of total read counts.
#
# At the moment the grid.predictions can only deal with simple regression
# models; the processing of 'multi' models must be implemented if necessary.
grid.predictions <- function(lll.M = M, M.multi = "simple", M.family = "wnlm.R", gene = "GRB10",
                             l.ylim = list(nlm.Q = c(-1, 9), nlm.R = c(-20, 120), nlm.S = c(0.4, 1.1), logi.S = c(0.4, 1), logi2.S = c(0.5, 1)),
                             xlim = c(-50, 200), n.pnts = 201, ...) {
    beta <- coef(M <- lll.M[[M.multi]][[M.family]][[gene]])
    # observed data
    if(M.family == "logi2.S") {
        obs <- lll.M[[M.multi]]$wnlm.S[[gene]]$model
        weights <- lll.M[[M.multi]]$wnlm.S[[gene]]$weights
    } else {
        obs <- M$model
        if(is.matrix(obs$Y)) {
            weights <- apply(obs$Y, 1, sum)
            obs$Y <- obs$Y[ , 1] / weights
        } else {
            weights <- M$weights
        }
    }
    n <- as.integer(mean(weights, na.rm = TRUE))
    # functions to calculate density
    pdf.nlm <- function(mu) {
        dnorm(y, mean = mu, sd = 1, log = TRUE)
    }
    pdf.logi <- function(mu, logi2 = FALSE) {
        y <- y # assign y from enclosing environment to a local y
        if(logi2)
            y <- 2 * (y - 0.5)
        dbinom(as.integer(y * n), size = n, prob = mu, log = TRUE)
    }
    # select function according to M.family
    pdfun <- switch(M.family,
                    wnlm.S = pdf.nlm, unlm.S = pdf.nlm,
                    wnlm.Q = pdf.nlm, unlm.Q = pdf.nlm,
                    wnlm.R = pdf.nlm, unlm.R = pdf.nlm,
                    logi.S = pdf.logi, logi2.S = function(x) pdf.logi(x, logi2 = TRUE))
    # set resolution
    x <- seq(xlim[1], xlim[2], length.out = n.pnts)
    if(grepl("nlm.Q", M.family))
        ylim <- l.ylim$nlm.Q # scale to percentiles
    if(grepl("nlm.R", M.family))
        ylim <- l.ylim$nlm.R # scale to percentiles
    if(grepl("nlm.S", M.family))
        ylim <- l.ylim$nlm.S
    if("logi.S" == M.family)
        ylim <- l.ylim$logi.S
    if("logi2.S" == M.family)
        ylim <- l.ylim$logi2.S
    y <- seq(ylim[1], ylim[2], length.out = max(n.pnts, n, na.rm = TRUE))
    # prediction: conditinal mean mu given x
    mu <- predict(M, newdata = data.frame(Age = x), type = "response")
    # this data frame will be returned after several modifications
    df <- as.data.frame(t(sapply(mu, pdfun)))
    #return(df)
    names(df) <- y.names <- paste("y", seq_along(y), sep = ".")
    df$x <- x
    if(M.family == "logi2.S")
        mu <- mu / 2 + 0.5 # scale down vertically by a factor of 2 and shift upwards
    df$mu <- mu
    df$multi <- factor(M.multi, levels = c("simple", "multiple"), ordered = TRUE)
    df$family <- factor(M.family, levels = c("logi.S", "logi2.S", "wnlm.S", "wnlm.Q", "wnlm.R", "unlm.S", "unlm.Q", "unlm.R"), ordered = TRUE)
    df$gene <- factor(gene, levels = names(lll.M[[M.multi]][[M.family]]), ordered = TRUE)
    df$n <- n
    # from wide format to long format
    df <- reshape(df, direction = "long", varying = y.names, v.names = "density", timevar = "y", times = y)
    df$density[abs(df$density) == Inf] <- NA
    # crucial scaling step
    df$density <- scale(df$density, scale = TRUE)
    # add observations
    df$obs.x <- NA
    df$obs.y <- NA
    df$obs.x[seq_along(obs$Age)] <- obs$Age
    df$obs.y[seq_along(obs$Y)] <- obs$Y
    return(df)
}

# wrapper around grid.predictions to apply that function to multiple model
# families for a given gene
grid.predictions.1gene <- function(gene = "GRB10", lll.M = M, M.multi = "simple",
                                   M.families = c(logi.S = "logi.S", logi2.S = "logi2.S", wnlm.S = "wnlm.S",
                                                  wnlm.Q = "wnlm.Q", wnlm.R = "wnlm.R", unlm.S = "unlm.S",
                                                  unlm.Q = "unlm.Q", unlm.R = "unlm.R")[1:4]) {
    do.call(rbind,
            lapply(M.families,
                   function(x) grid.predictions(lll.M = M, M.multi = "simple", M.family = x, gene = gene)))
}

# A levelplot showing, scaled probability densities and predictions
#
#
# Arguments
# df: a thin, long data frame obtained with grid.predictions
# fm: a formula
# show.model: whether to show model desnsity besides the data
#
# Details
# A rectangle with dotted line shows demarcates the theoretically possible
# region (non-negative age and S between 0.5 and 1 or S between 0 and 100).
plot.predictions <- function(l.df, gene, fm = formula(density ~ x * y | family),
                             show.model = TRUE, ...) {
    df <- l.df[[gene]]
    levelplot(fm, data = df, M.family = df$family,
              panel = function(x, y, subscripts, M.family, ...) {
                  if(show.model) {
                      panel.levelplot(x, y, subscripts = subscripts, contour = TRUE, lwd = 0.5, col = "gray", ...)
                      panel.xyplot(x = x[subscripts][nx <- seq_along(unique(x))], y = df$mu[subscripts][nx],
                                   col = "black", type = "l", lwd = 2, ...)
                  }
                  panel.xyplot(x = df$obs.x[subscripts], y = df$obs.y[subscripts],
                               pch = 21, col = "darkgreen", fill = "green", alpha = 0.5, cex = 0.5, ...)
                  is.Q <- grepl("nlm.Q", M.family[subscripts][1])
                  panel.rect(xleft = 0, ybottom = ifelse(is.Q, 0, 0.5),
                             xright = 1e5, ytop = ifelse(is.Q, 100, 1), lty = "dotted")
              },
              par.settings = list(panel.background = list(col = "gray"), regions = list(col = rev(trellis.par.get("regions")$col))),
              scales = list(y = list(relation = "free")), colorkey = FALSE,
              strip = strip.custom(strip.levels = c(show.model, FALSE)),
              ylim = list(c(0.44, 1.06), c(0.44, 1.06), c(0.44, 1.06), c(0, 8)),
              xlab = "age", ylab = "S: read count ratio;  Q: transformed S",
              cut = 30,
              main = if(show.model) paste("fixed simple regression model,", gene) else gene,
              ...)
}
