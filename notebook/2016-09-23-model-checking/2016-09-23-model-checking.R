get.diagnostics <- function(l.M, mtype = "logi.S") {
    helper <- function(m, gene = "PEG3", mtype = "logi.S") {
        df <- data.frame(gene = gene, model.type = mtype,
                         fitted = fitted(m),
                         leverage = hatvalues(m),
                         cooks.dist = cooks.distance(m),
                         res.std.dev = rstandard(m, type = "deviance"),
                         res.std.pear = rstandard(m, type = "pearson"))
        df$res.combined <- with(df, res.std.dev + log(res.std.pear / res.std.dev) / res.std.dev)
        return(df)
    }
    do.call(rbind,
            lapply(names(l.M), function(g) helper(l.M[[g]], gene = g, mtype = mtype)))
}

get.influence <- function(m, cases) {
    x <- rep(NA, length(cases))
    names(x) <- cases
    y <- cooks.distance(m)
    x[names(y)] <- y
    return(x)
}

myqqnorm2 <- function(genes = gene.ids, models = c("wnlm.Q", "unlm.Q"), dt = diagnostics, pch = if(length(genes) > 3) c("o", "+") else 16, ...) {
    col <- if(length(models) == 2) c("purple", "red")
    qqmath(~ res.std.dev | gene, data = dt, groups = model.type, subset = dt$gene %in% genes & dt$model.type %in% models,
           pch = pch, par.settings = list(superpose.symbol = list(col = col)),
           abline = c(0, 1), auto.key = list(text = rev(models), col = col, points = FALSE),
           ...)
}
