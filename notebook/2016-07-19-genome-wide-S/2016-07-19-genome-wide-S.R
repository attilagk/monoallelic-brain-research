emp.distr.S <-
    function(Y,
             ss = seq(0.5, 1, length.out = 101),
             with.density = FALSE,
             ...) {
        E <- list()
        E$ecdf <- lapply(Y, ecdf)
        E$ecdf.val <- data.frame(lapply(E$ecdf, function(f) f(ss)))
        E$ss <- ss
        if(with.density) {
            E$density <-
                data.frame(lapply(Y, function(x)
                                  density(x, na.rm = TRUE, n = length(ss), from = 0.5, to = 1, ...)$y))
        }
        return(E)
}

ecdf.at.some.s <-
    function(ED, sel.g, type = "p", pch = 21, xlab = "F(s)", par.settings = list(dot.line = list(lty = 0)), ...){
        dotplot(~
                ED$ecdf.val[101, rev(sel.g)] +
                ED$ecdf.val[81, rev(sel.g)] +
                ED$ecdf.val[61, rev(sel.g)] +
                ED$ecdf.val[41, rev(sel.g)],
            type = type, pch = pch, xlab = xlab, par.settings = par.settings,
            auto.key = list(text = paste("s =", c(0.7, 0.8, 0.9, 1)), columns = 4),
            main = "Empirical CDF F at the following s values:",
            #scales = list(y = list(draw = FALSE)),
            ...)
    }

# complex figure on the distribution of S_ig

# reshape Y given selected genes
reshape.Y <- function(Y, sel.g = c("PEG10", "ZNF331", "AFAP1")){
    Y.l <- reshape(Y[sel.g], direction = "long", varying = list(sel.g), v.names = "s",
                   times = sel.g, timevar = "gene")
    Y.l$gene <- factor(Y.l$gene, levels = sel.g, ordered = TRUE)
    return(Y.l)
}

# separate plot strategy

# new density plot
my.densityplot <- function(Y.long, ...) {
    panel.my.densityplot <- function(...) {
        panel.densityplot(...)
        panel.grid(h = 0, v = -1)
    }

    densityplot(~ s, data = Y.long, groups = gene, plot.points = FALSE, type = c("p"),
                panel = panel.my.densityplot,
                scales = list(x = list(draw = FALSE), y = list(rot = 90, relation = "free")),
                xlim = c(0.5, 1),
                par.settings = list(superpose.line = standard.theme(color = FALSE)$superpose.line),
                auto.key = list(corner = c(0.1, 0.9), lines = TRUE, points = FALSE),
                xlab = NULL, ylab = "density")
}

# new ECDF plot
my.ecdfplot <- function(Y.long, eval.at = c(0.9), ...) {
    panel.my.ecdfplot <- function(x, groups, subscripts, eval.at = eval.at, ...){
        panel.superpose(x, groups = groups, subscripts = subscripts,
                        panel.groups = panel.ecdfplot, type = c("s"), ...)
        #panel.abline(v = eval.at, lty = 2, col = "blue")
        panel.grid(h = 0, v = -1)
        g.lev <- levels(groups)
        ss <- rep(c(eval.at,
                    rep(NA, length(x) / length(g.lev) - length(eval.at))),
                  length(g.lev))
        yy <- sapply(levels(groups), function(g) ecdf(x[groups == g])(ss))
        panel.superpose(x = ss, y = yy, groups = groups, subscripts = subscripts,
                        panel.groups = panel.xyplot, col = "blue", ...)
    }

    ecdfplot(~ s, data = Y.long, groups = gene, eval.at = c(0.9),
             scales = list(x = list(draw = FALSE), y = list(rot = 90, relation = "free")),
             panel = panel.my.ecdfplot,
             par.settings = list(superpose.line = standard.theme(color = FALSE)$superpose.line),
             xlim = c(0.5, 1),
             ylab = "ECDF, F(s)", xlab = NULL)

}

# ECDF genome-wide
my.levelplot <- function(ED.long, pct.top.g = 2, n.all.g = length(ok.genes), ...) {

    # data manipulation
    n.top.g <- ceiling(n.all.g * pct.top.g / 100)
    ED.long$rank.segment <- factor(rep(1, nrow(ED.long)), levels = c("top", "bottom"), ordered = TRUE)
    ED.long$rank.segment[ with(ED.long, gene %in% levels(gene)[seq_len(n.top.g)]) ] <- "top"
    ED.long$rank.segment[ with(ED.long, ! gene %in% levels(gene)[seq_len(n.top.g)]) ] <- "bottom"

    panel.my.levelplot <- function(...) {
        panel.levelplot(...)
        panel.grid(h = 0, v = -1)
    }

    lp <- levelplot(ECDF ~ s * gene | rank.segment, data = ED.long, legend = NULL, colorkey = FALSE,
                    panel = panel.my.levelplot,
                    scales =
                        list(y =
                             list(rot = 90, relation = "free",
                                  limits = list(c(n.top.g, 1),
                                                c(n.all.g, n.top.g + 1)))),
                    xlab = "imbalance score, s", ylab = "gene rank",
                    layout = c(1, 2))

    dimnames(lp)$rank.segment <- c(paste("top", pct.top.g, "% of genes"),
                                   paste("bottom", 100 - pct.top.g, "% of genes"))
    return(lp)
}



plot.all <- function(plots) {
    print(plots$density, position = c(0.0, 0.85, 1.0, 1.0), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$ecdf, position = c(0.0, 0.7, 1.0, 0.85), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$level, position = c(0.0, 0.0, 1.0, 0.7), more = FALSE)
}

