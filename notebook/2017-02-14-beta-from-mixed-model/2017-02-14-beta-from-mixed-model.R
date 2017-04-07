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
           scales = switch(group.by, gene = "same", coefficient = "free", coefficient.1 = "free", predictor = "free"),
           layout = switch(group.by, gene = c(5, 6), coefficient = c(4, 6), coefficient.1 = c(4, 6), predictor = c(4, 4)),
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

tval.vp.1gene <- function(gene, m.type, vp, llm = M, e.v = e.vars, ref.m = M$unlm.Q[[1]]) {
    lm <- llm[[m.type]]
    m <- lm[[gene]]
    if(is.null(ref.m)) ref.m <- m
    coefs <- rownames(summary(ref.m)$coefficients)[-1]
    sm <- summary(m)$coefficients[ -1, ]
    t.value <- sapply(coefs, function(x) ifelse(x %in% rownames(sm), sm[x, "t value"], NA))
    df <- df.residual(m)
    var.part <- expand.x.preds2coefs(unlist(vp[[m.type]][gene, e.v]), ref.m)
    predictor <- factor(expand.x.preds2coefs(e.v, ref.m), ordered = TRUE, levels = e.v)
    coefficient <- factor(coefs, ordered = TRUE, levels = coefs)
    data.frame(gene = gene, predictor = predictor, var.part = var.part, coefficient = coefficient, t.value = t.value, df = df)
}

tval.vp <- function(m.type = "fixed.1", llm = M) {
    df <- Reduce(rbind, lapply(names(llm[[ m.type ]]), tval.vp.1gene, m.type, vp, llm))
    df$model <- m.type
    return(df)
}

tval.vp.plot <- function(df, lm = M$fixed.1, ...) {
    xyplot(abs(t.value) ~ var.part | gene, data = df,
           panel = function(..., pch, col, d.f, subscripts) {
               pch <- as.character(df$coefficient)
               col <- rainbow(length(levels(df$predictor)))[df$predictor]
               d.f <- df$df[subscripts]
               panel.abline(h = qt(c(5e-2, 1e-2) / 2, df = d.f, lower.tail = FALSE), lty = 2, lwd = 0.5, col = c("white", "black"))
               panel.text(..., pch = pch, cex = 0.7, col = col)
           },
           scales = "free", layout = c(4, 4), par.settings = list(background = list(col = "#909090")),
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
             list(columns = 4, text = list(c(paste0("(", seq_along(levels(cf$coefficient)[-1]), ") ", levels(cf$coefficient)[-1]), "signif.level 0.05", "signif.level 0.01"), cex = 0.7, col = c(expand.x.preds2coefs(rnb, M[[1]][[1]]), "white", "black"))),
         coefficient.2 =
             list(columns = 3, text = list(c(paste0("(", seq_along(levels(cf$coefficient)[-1]), ") ", levels(cf$coefficient)[-1]), "signif.level 0.05", "signif.level 0.01"), cex = 0.7, col = c(expand.x.preds2coefs(rnb, M[[1]][[1]]), "white", "black"))),
         predictor =
             list(columns = 4, text = list(paste0("(", seq_along(e.vars), ") ", e.vars), cex = 0.7)),
         gene =
             list(columns = 5, text = list(paste0("(", seq_along(levels(cf$gene)), ") ", levels(cf$gene)), cex = 0.7)))

