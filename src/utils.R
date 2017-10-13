
CI.p <- function(p.hat, n, conf.lev = 0.95) {
    alpha <- (1 - conf.lev) / 2
    zz <- lapply(list(lower = alpha, upper = 1 - alpha), qnorm)
    std <- sqrt(p.hat * (1 - p.hat) / n)
    lapply(zz, function(z) p.hat + z * std)
}

pad.fraction.imprint <- function(gene = "PWAR6", imp.stat = genes.fig1, fr = fraction) {
    padding <- as.integer(imp.stat[gene, "imprinting.status"]) - 1
    c(rep(rep(0, 6), padding), fr[[gene]], rep(rep(0, 6), 2 - padding))
}

padded.frac <- function(fr = frac$min.reads.15[51:1], imp.stat = gene.summary) {
    genes <- names(fr)
    names(genes) <- genes
    df <- data.frame(lapply(genes,
                            pad.fraction.imprint, imp.stat = imp.stat, fr = fr),
                     check.names = FALSE)
    t(as.matrix(df))
}



palette.ifat <- function(is.bw = FALSE, cols = if(is.bw) rep("black", 3) else c("blue", "darkgreen", "red"), detailed = FALSE, col2 = "tan") {
    color.fun <- function(col, detailed) {
        if(detailed)
            c(col, colorRampPalette(c("grey50", "grey90"))(4), col2)
        else
            c(colorRampPalette(c(col, "lightgray"))(3), "lightgray", "lightgray", col2)
    }
    is.bw
    c(sapply(cols, color.fun, detailed = detailed))
}

filter.min.obs <- function(min.obs = 10, X) {
    X[sapply(X, function(y) sum(! is.na(y)) >= min.obs)]
}

# filter X based on N
filter.min.read <- function(min.read = 15, X, N, min.obs = 0, ...) {
    genes <- names(X)
    names(genes) <- genes
    df <- data.frame(lapply(genes,
                function(g, Y = X) {
                    Y[[g]][N[[g]] < min.read] <- NA
                    return(Y[[g]])
                }, Y = X), check.names = FALSE)
    filter.min.obs(min.obs = min.obs, X = df)
}

# make a list of ECDFs ordered according to each ECDF component evaluated at sort.s
#
# Arguments
# S: a list of data frames of S[i,g], with rows i for individuals and columns g for genes
# sort.s: the s at which each ECDF F_g must be evaluated for sorting, i.e. F_g(s)
#
# Value
# a list of ECDFs where the inner lists (whose components are ECDFs
# F_g : for all genes g) are sorted according to F_g(s)
sorted.ecdfs <- function(S, sort.s = 0.9) {
    l.E <- lapply(S, ecdf)
    l.E[order(sapply(l.E, function(f) f(sort.s)))]
}

cum.frac.ecdf <- function(l.E, frac = 10:6 / 10) {
    cum.frac <- data.frame(lapply(l.E, function(f) f(10:6 / 10)), check.names = FALSE)
    row.names(cum.frac) <- as.character(10:6 / 10)
    cum.frac
}

cum.frac.andys.test <- function(genes, S, N, ucl.fun = CI.p, max.ucl = 0.7, max.s = 0.6) {
    helper <- function(g)
        sum(S[[g]] <= max.s & ucl.fun(S[[g]], N[[g]])$upper < max.ucl, na.rm = TRUE) / sum(! is.na(S[[g]]))
    names(genes) <- genes
    df <- data.frame(lapply(genes, helper), check.names = FALSE)
    row.names(df) <- "andys.test"
    return(df)
}

do.fractions <- function(l.E, S, N, frac = 10:6 / 10, ucl.fun = CI.p, max.ucl = 0.7, max.s = 0.6) {
    cum.frac <- rbind(cum.frac.ecdf(l.E, frac),
                      cum.frac.andys.test(names(l.E), S, N, ucl.fun, max.ucl, max.s))
    frac <- data.frame(lapply(cum.frac, function(y) - diff(c(y, 0))), check.names = FALSE)
    row.names(frac) <- row.names(cum.frac)
    frac
}
