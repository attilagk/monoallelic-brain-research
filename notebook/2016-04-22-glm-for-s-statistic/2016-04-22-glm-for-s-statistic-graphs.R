########################## FUNCTIONS #########################################

# Fits simple models S ~ x2 (ignoring all other explanatory variables). Returns a list of fitted models
fit.S.x1 <- function(g = genes.or.gsets[1], dt = d, x1 = "Age.of.Death") {
    md.l <- list()
    md.l$nlm.R <- glm(mk.form(paste0("R_", g), x1), data=dt)
    md.l$nlm2.R <- lm(mk.form(paste0("R_", g), x1), data=dt)
    md.l$nlm.S <- glm(mk.form(paste0("S_", g), x1), data=dt)
    md.l$logi.S <- glm(mk.form(paste0("C_", g), x1), data=dt, family=binomial)
    md.l$logi2.S <- glm(mk.form(paste0("C2_", g), x1), data=dt, family=binomial)
    return(md.l)
}
# Plots simple models S ~ x1 in model list 'md.l', supposed to contain "nlm.S", "logi.S", "logi2.S"
# An offset of the linear predictor can be specified by 'offset'.  The RNA
# quality variables RIN and RIN.2 are corrected for by adding fun.RIN(RIN) and
# fun.RIN(RIN.2) to the offset; by default 'fun.RIN' is the average.
# Data are also plotted
plot.S.x1 <- function(md.l, xlimits = c(0, 600), ylimits = c(0, 1),
                      x1 = "Age.of.Death", offset = c(0,0,0), add = FALSE,
                      fun.RIN = mean, ...) {
    # ensure that nlm.S, logi.S, logi2.S are contained in md.l
    stopifnot(length(intersect(names(md.l), md.types <- c("nlm.S", "logi.S", "logi2.S"))) == 3)
    names(offset) <- md.types
    # offset according to RIN
    rin <- "DLPFC_RNA_isolation..RIN"
    rin2 <- "DLPFC_RNA_isolation..RIN.2"
    if (length(grep(rin, names(coef(md.l[[1]]))))) {
        off <- sapply( md.l, function(m) coef(m)[[rin]] * fun.RIN(m$data[[rin]]))
        off2 <- sapply( md.l, function(m) coef(m)[[rin2]] * fun.RIN(m$data[[rin2]]))
        # add to offset
        offset <- offset + off + off2
    }
    # create a list of link functions
    lf <- list()
    lf$nlm.S <-
        function(x) md.l[["nlm.S"]]$coefficients[1] + x * md.l[["nlm.S"]]$coefficients[x1] + offset["nlm.S"]
    lf$logi.S <-
        function(x) 1 / (1 + exp( - md.l[["logi.S"]]$coefficients[1] - x * md.l[["logi.S"]]$coefficients[x1] - offset["logi.S"]) )
    lf$logi2.S <-
        function(x) 0.5 + 0.5 / (1 + exp( - md.l[["logi2.S"]]$coefficients[1] - x * md.l[["logi2.S"]]$coefficients[x1] - offset["logi2.S"]) )
    if (! add) {
        # plot data
        fm <- mk.form(as.character(formula( md.l[[ md.types[1] ]] ))[2], x1)
        plot(fm, data=md.l[[1]]$data, pch=".", col="black", xlim=xlimits,
             ylim=ylimits, cex=1.5)
    }
    # 'baseline' guide lines
    if (all(ylimits == c(0, 1))) {
        abline(h=0.5, lty=3); abline(h=0, lty=3); abline(h=1, lty=3)
    }
    # plot model lines
    sapply(seq_along(lf), function(i)
           plot(lf[[i]], from=xlimits[1], to=xlimits[2], add=TRUE, col=c(2,4,3)[i], ...))
    return(lf)
}

# Plots response stat Y against explanatory variable x1.
# When 'ar' is TRUE an array of plots is produced presenting genes on separate
# plots; when 'ar' is FALSE data for all genes are shown in a single plot with
# color code
plot.smooth.Y.x1 <- function(Y = "S", ar = TRUE, x1 = "Age.of.Death") {
    bar <- function(g) {
        if(ar) {
            d.x <- d$Age.of.Death
            d.y <- d[[ paste0(Y, "_", g) ]]
        } else {
            d.x <- d.p$pool16$Age.of.Death
            d.y <- d.p$pool16[[ paste0(Y, "_pool16") ]]
        }
        plot(d.x, d.y,
             type="n", xlab = ifelse(ar, "", x1),
             ylab = ifelse(ar, "", Y),
             main = ifelse(ar, g, ""),
             mar = ifelse(ar, c(0,0,0,0), c(5,4,4,2)),
             ylim = c(quantile(d.y, 1e-2, na.rm=TRUE), max(d.y, na.rm=TRUE))
             )
    }
    foo <- function(i) {
        g <- genes[i]
        if(ar) bar(g)
        S.sm <- lowess(age <- d$Age.of.Death,
                       ifelse(is.na(s.g <- d[[ paste0(Y, "_", g) ]]),
                              mean(s.g, na.rm=TRUE), s.g), f = 1/2)
        points(age, s.g, col = ifelse(ar, 1, i + 1), pch=".", cex=1.5)
        lines(S.sm, col = ifelse(ar, 2, i + 1))
    }
    if(ar) {
        par(mfrow=c(4,4), mar=c(2,2,3,1))
        par(cex=1, cex.lab=1)
    }
    else {
        bar("")
    }
    sapply(seq_along(genes), foo)
    return(invisible())
}

# compare model 1 to model 2; ml1 and ml2 are lists of models fitted to
# different data; mf1 and mf2 are strings of model families ("nlm.R", "logi.S", "logi2.S")
# fun compares the models's AIC values, default is scaled difference
m2m <- function(ml1, ml2, mf1, mf2, fun = function(x, y) (x - y) / abs(x)) {
    sapply(genes.or.gsets, function(r) fun(ml1[[r]][[mf1]]$aic, ml2[[r]][[mf2]]$aic))
}

plot.simple.cmp.multiple <- function(g = "PEG3") {
    par(mfrow=c(1,2))
    plot.S.x1(fit.S.x1(g), xlimits = c(0, 300), ylimits = c(0.7, 1))
    title(main=expression(paste(beta[0], ", ", beta[age], " from simple regr.")))
    legend("left", c("nlm.S", "logi.S", "logi2.S"), col=c(2,4,3), lty=1)
    plot.S.x1(m[[g]], xlimits = c(0, 300), ylimits = c(0.7, 1))
    plot.S.x1(m[[g]], add=TRUE, fun.RIN=min, lty=2, xlimits = c(0, 300), ylimits = c(0.7, 1))
    title(main=expression(paste(beta[0], ", ", beta[age], " from multiple. regr.")))
    legend("left", c("avg RIN", "nlm.S", "logi.S", "logi2.S"), col=c(NA, 2,4,3), lty=1)
    legend("bottomleft", c("min RIN", "nlm.S", "logi.S", "logi2.S"), col=c(NA, 2,4,3), lty=2)
}

# Gets property 'prop' of 'coef' from both md.l1 and md.l2 model lists
# The property 'prop' is an index for the vector c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
coefs.from.md.lists <- function(md.l1, md.l2,
                                prop = 1,
                                coef = "Age.of.Death") {
    md.types <- c("nlm.R", "nlm2.R", "nlm.S", "logi.S", "logi2.S")
    names(md.types) <- md.types
    lapply(md.types,
           function(ty) {
               m1 <- sapply(genes.or.gsets,
                            function(g) summary(md.l1[[g]][[ty]])$coefficients[coef, prop])
               m2 <- sapply(genes.or.gsets,
                            function(g) summary(md.l2[[g]][[ty]])$coefficients[coef, prop])
               cbind(m1, m2)
           })
}

plot.beta.pval <- function(mod = "nlm.S", bts = betas, pvs = pvals,
                           xlim.b=c(-2e-3,2e-3), xlim.logp=c(0,20), ...) {
    par(mar = c(5, 6, 4, 2), mfrow = c(1, 2))
    barplot(t(bts[[ mod ]]), beside=TRUE, horiz=TRUE, las=1, xlim=xlim.b,
            legend.text=c("simple", "multiple"),
            xlab=expression(hat(beta)[age]), main="estim. effect of age",
            col=gray(c(7,2)/8), ...)
    barplot(- log10(t(pvs[[ mod ]])), beside=TRUE, horiz=TRUE, las=1,
            xlim=xlim.logp, xlab=expression(-log[10](p)), main="significance",
            col=gray(c(7,2)/8))
}

# Extracts a matrix of effects under model type 'md.ty' for each gene or
# aggregated set
fx.summary <- function(md.l, md.ty="nlm.S", coefs = names(coef(m[[1]][[md.ty]])) ) {
    f <- function(cf)
        sapply(md.l, function(x) effects(x[[ md.ty ]])[ cf ])
    sapply(coefs, f)
}

# Extracts a matrix of deviance from ANOVA under model type 'md.ty' for each gene or
# aggregated set
anova.summary <- function(md.l, md.ty="nlm.S", coefs = row.names(anova(m[[1]][[md.ty]]))) {
    f <- function(cf)
        sapply(md.l, function(x) anova(x[[ md.ty ]])[ cf, 'Deviance' ])
    sapply(coefs, f)
}

# Prettify coefficient names by removing long substrings
nice.cf.names <- function(coef.n = names(coef(m[[1]][[1]]))) {
    sub("DLPFC_RNA_(report..Clustered.Library|isolation)\\.+", "RNA.", coef.n)
}

# Produces a boxplot of the matrix of effects.
# Also works with the matrix of deviances from ANOVA.
fx.summary.boxplot <- function(eff, names=nice.cf.names()[-1], ...) {
    par(mar = c(5, 8, 4, 2), mfrow = c(1, 1))
    boxplot(eff, horizontal=TRUE, las=1, add=FALSE, names=names, ...)
    grid()
    abline(v=0, col="red", lty="solid")
}

########################## CALCULATE A FEW OBJECTS ###########################

