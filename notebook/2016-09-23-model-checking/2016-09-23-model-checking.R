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
