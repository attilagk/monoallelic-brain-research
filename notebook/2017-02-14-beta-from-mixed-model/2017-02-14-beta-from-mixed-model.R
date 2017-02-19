# high level plotting function
my.plot <- function(x = "unlm.Q", y = "fixed.1", group.by = "gene", lbl.type = "coefficient", dt = cf, ...) {
    my.panel <- function(..., lbl) {
        panel.abline(a = 0, b = 1, lty = "dotted", col = "gray")
        lbl <-
            if(lbl.type == "gene") dt$gene
            else dt[[1]]
        panel.text(..., pch = lbl, cex = 0.7)
    }
    fm <- formula(paste(y, "~", x, "|", group.by))
    xyplot(fm, data = dt, panel = my.panel,
           scales = switch(group.by, gene = "same", coefficient = "free", predictor = "free"),
           layout = switch(group.by, gene = c(5, 6), coefficient = c(4, 6), predictor = c(4, 4)),
           key = my.key[[lbl.type]], ...)
}

# Expand "predictors"-length vector x to "coefficients"-length.
# When the components of vector x are named according to predictors, expand
# (repeat) each component to match the number of corresponding coefficients
expand.x.preds2coefs <- function(x, m) {
    helper <- function(e.v) {
        cf <- predictor2coefs(m, e.v)
        xx <- rep_len(x[e.v], length(cf))
        names(xx) <- cf
        return(xx)
    }
    unlist(lapply(names(x), helper))
}

tval.vp.1gene <- function(gene, m.type, vp, llm = M, e.v = e.vars) {
    lm <- llm[[m.type]]
    data.frame(t.value = tv <- summary(lm[[gene]])$coefficients[ -1, "t value"],
               var.part = expand.x.preds2coefs(unlist(vp[[m.type]][gene, e.v]), lm[[gene]]),
               predictor = factor(expand.x.preds2coefs(e.v, lm[[gene]]), ordered = TRUE, levels = e.v),
               gene = gene, coefficient = factor(names(tv), ordered = TRUE, levels = names(tv)))
}

tval.vp <- function(m.type = "fixed.1", llm = M) {
    df <- Reduce(rbind, lapply(names(llm[[ m.type ]]), tval.vp.1gene, m.type, vp, llm))
    df$model <- m.type
    return(df)
}

tval.vp.plot <- function(df, ...) {
    xyplot(abs(t.value) ~ var.part | gene, data = df,
           panel = function(..., pch, col) {
               pch <- as.character(df$coefficient)
               col <- rainbow(length(levels(df$predictor)))[df$predictor]
               panel.text(..., pch = pch, cex = 0.7, col = col)
           },
           scales = "free",
           key = my.key$coefficient.1,
           ...)
}

rnb <- rainbow(length(e.vars))
names(rnb) <- e.vars

# annotation keys
my.key <-
    list(coefficient =
             list(columns = 5, text = list(paste0("(", seq_along(levels(cf$coefficient)), ") ", levels(cf$coefficient)), cex = 0.7)),
         coefficient.1 =
             list(columns = 4, text = list(paste0("(", seq_along(levels(cf$coefficient[-1])), ") ", levels(cf$coefficient)[-1]), cex = 0.7, col = expand.x.preds2coefs(rnb, M[[1]][[1]]))),
         predictor =
             list(columns = 4, text = list(paste0("(", seq_along(e.vars), ") ", e.vars), cex = 0.7)),
         gene =
             list(columns = 5, text = list(paste0("(", seq_along(levels(cf$gene)), ") ", levels(cf$gene)), cex = 0.7)))

