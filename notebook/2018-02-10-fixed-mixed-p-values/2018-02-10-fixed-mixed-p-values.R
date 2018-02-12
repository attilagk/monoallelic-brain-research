# based on get.both.p.vals of 2018-02-10-fixed-mixed-p-values.R
extract.p.val <- function(mtype = "wnlm.Q", gene = "MEST", coef = "Age", M) {
    errfun <- function(e) NA
    data.frame(Model = mtype, Gene = gene, Coefficient = coef,
               p.val = tryCatch(summary(M[[c(mtype, gene)]])$coefficients[coef , 4],
                                  error = errfun))
}


my.dotplot <- function(do.mixed = TRUE, do.unlm.Q = ! do.mixed, df = d, param = pars,
                       col = list(wnlm.Q = "red", unlm.Q = "purple", b = "green4", b_g = "blue" ),
                       ...) {
    sset <- if(do.mixed) TRUE else df$Effects == "fixed"
    sset <- if(do.unlm.Q) sset else sset & df$Model != "unlm.Q"
    dotplot(Term ~ log10(p) | Predictor, groups = Model, data = df, layout = c(3, 1), subset = sset,
            xlab = expression(paste(log[10], "p")),
            auto.key = list(text =
                            if(do.mixed)
                                if(do.unlm.Q)
                                    rep("", 3)
                                else
                                    c("mixed: gene specific", "mixed: gene nonspecific", "earlier fixed effects model")
                            else
                                if(do.unlm.Q)
                                    c("", "unlm.Q", "wnlm.Q")
                                else
                                    c("", "", "earlier fixed effects model") ,
                            col = c(col$b_g, if(do.unlm.Q & ! do.mixed) col$unlm.Q else col$b, col$wnlm.Q),
                            points = FALSE),
            scales = list(y = list(limits = seq(from = 0, to = 1 + length(levels(df$Term))), at = seq_along(levels(df$Term)),
                                   labels = c(paste0("beta_", param$fixed), if(do.mixed) param$mixed else rep("", 3))
                                   ),
                          x = list(relation = "free", limits = list(c(-6, 1), c(-6, 1), c(-18, 3)))
            ),
            between = list(x = 0.5),
            par.settings = list(superpose.symbol = list(pch = 16, col = list(col$wnlm.Q, col$unlm.Q, c(col$b_g, col$b))),
                                dot.symbol = list(col = c("red", "blue", "green4")),
                                ...))
}
