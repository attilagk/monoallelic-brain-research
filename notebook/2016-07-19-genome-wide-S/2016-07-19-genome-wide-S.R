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
        panel.my.dotplot <- function(..., sel.g) {
            panel.grid(h = -1, v = 0)
            #panel.grid(h = length(sel.g) - 0, v = 0)
            panel.dotplot(...)
        }
        dotplot(~
                ED$ecdf.val[101, rev(sel.g)] +
                ED$ecdf.val[81, rev(sel.g)] +
                ED$ecdf.val[61, rev(sel.g)] +
                ED$ecdf.val[41, rev(sel.g)],
            panel = panel.my.dotplot, sel.g = sel.g,
            type = type, pch = pch, xlab = xlab, par.settings = par.settings,
            auto.key = list(text = paste("s =", rev(c(0.7, 0.8, 0.9, 1))), columns = 4),
            main = "ECDF F at the following s values:",
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
    sel.g <- levels(Y.long$gene)
    panel.my.densityplot <- function(...) {
        panel.densityplot(...)
        panel.grid(h = 0, v = -1)
    }

    densityplot(~ s, data = Y.long, groups = gene, plot.points = FALSE, type = c("p"),
                panel = panel.my.densityplot,
                scales = list(x = list(draw = FALSE), y = list(rot = 90, relation = "free")),
                xlim = c(0.5, 1),
                par.settings = list(superpose.line = standard.theme(color = FALSE)$superpose.line),
                auto.key = list(text = paste0(letters[seq_along(sel.g)], ": ", sel.g),
                                corner = c(0.1, 0.9), lines = TRUE, points = FALSE),
                xlab = NULL, ylab = "density est.")
}

# new ECDF plot
my.ecdfplot <- function(Y.long, eval.at = c(0.9), ...) {
    df <- expand.grid(x = c(0.5, 0.55), y = 0:100 / 100)
    df$z <- rep(0:100 / 100, each = 2)
    df$z[seq(2, 202, by = 2)] <- NA
    panel.my.ecdfplot <- function(x, groups, subscripts, eval.at = eval.at, df, ...){
        with(df, panel.levelplot(x = x, y = y, z = z, at = 0:100 / 100, subscripts = 0:202, ...))
        panel.superpose(x, groups = groups, subscripts = subscripts,
                        panel.groups = panel.ecdfplot, type = c("s"), ...)
        #panel.abline(v = eval.at, lty = 2, col = "blue")
        panel.grid(h = 0, v = -1)
        g.lev <- levels(groups)
        ss <- rep(c(eval.at,
                    rep(NA, length(x) / length(g.lev) - length(eval.at))),
                  length(g.lev))
        yy <- sapply(levels(groups), function(g) ecdf(x[groups == g])(ss))
        gg <- sapply(levels(groups), function(g) ecdf(x[groups == g])(ss))
        #panel.superpose(x = ss, y = yy, groups = groups, subscripts = subscripts,
        #                panel.groups = panel.xyplot, col = "blue", pch = 21, ...)
        lapply(levels(groups),
               function(g) panel.xyplot(x = eval.at, y = ecdf(x[groups == g])(eval.at),
                                        col = "blue", ...))
        lapply(seq_along(levels(groups)),
               function(i)
                   panel.xyplot(x = eval.at - 0.02, y = ecdf(x[groups == levels(groups)[i]])(eval.at),
                                pch = letters[i], col = "black", cex = 1.5, ...))
    }

    ecdfplot(~ s, data = Y.long, groups = gene, eval.at = c(0.9),
             df = df,
             scales = list(x = list(draw = FALSE), y = list(rot = 90, relation = "free")),
             panel = panel.my.ecdfplot,
             par.settings = list(superpose.line = standard.theme(color = FALSE)$superpose.line),
             xlim = c(0.5, 1),
             ylab = "ECDF", xlab = NULL)

}

panel.my.levelplot <- function(...) {
    panel.levelplot(...)
    panel.grid(h = 0, v = -1)
}

my.levelplot.sel.g <- function(ED.long, sel.g = c("PEG10", "ZNF331", "AFAP1"), ...) {
    ED.long$gene <- with(ED.long, factor(gene, levels = rev(levels(gene)), ordered = TRUE))
    levelplot(ECDF ~ s * gene, data = ED.long[seq_len(nrow(ED.long)), ], legend = NULL, colorkey = FALSE,
                    panel = panel.my.levelplot, subset = gene %in% rev(sel.g),
                    xlab = NULL,
                    scales =
                        list(y =
                             list(rot = 90, relation = "free",
                                  labels = letters[rev(seq_along(sel.g))]
                                                ),
                             x = list(draw = FALSE)),
                    ...)
}

# ECDF genome-wide
my.levelplot <- function(ED.long, pct.top.g = 2, n.all.g = length(ok.genes), ...) {
    # data manipulation
    n.top.g <- ceiling(n.all.g * pct.top.g / 100)
    ED.long$rank.segment <- factor(rep(1, nrow(ED.long)), levels = c("top", "bottom"), ordered = TRUE)
    ED.long$rank.segment[ with(ED.long, gene %in% levels(gene)[seq_len(n.top.g)]) ] <- "top"
    ED.long$rank.segment[ with(ED.long, ! gene %in% levels(gene)[seq_len(n.top.g)]) ] <- "bottom"

    lp <- levelplot(ECDF ~ s * gene | rank.segment, data = ED.long, legend = NULL, colorkey = FALSE,
                    panel = panel.my.levelplot,
                    scales =
                        list(y =
                             list(rot = 91, relation = "free",
                                  limits = list(c(n.top.g, 1),
                                                c(n.all.g, n.top.g + 1)))),
                    xlab = "imbalance score, s", ylab = "gene rank",
                    layout = c(1, 2))

    dimnames(lp)$rank.segment <- c(paste("top", pct.top.g, "% of genes"),
                                   paste("bottom", 100 - pct.top.g, "% of genes"))
    return(lp)
}

rankplot <-
    function(ED, gene.order, pct.top.g = 2, sel.g, ...) {
        n.top.g <- ceiling(length(gene.order) * pct.top.g / 100)
        df <- data.frame(ECDF = sapply(ED$ecdf[gene.order], function(f) f(0.9)))
        ordered.genes <- names(ED$ecdf)[gene.order]
        df$gene <- factor(ordered.genes, levels = c(rev(ordered.genes)), ordered = TRUE)
        df$pch <- rep("", nrow(df))
        df$pch[with(df, gene %in% sel.g)] <- letters[seq_along(sel.g)]
        df$rank.segment <-
            factor(c(rep("top", n.top.g),
                     rep("bottom", length(gene.order) - n.top.g)),
                   levels = c("top", "bottom"), ordered = TRUE)

        panel.rankplot <- function(x, y, groups, subscripts, ...) {
            panel.superpose(x, y, groups = groups, subscripts = subscripts, col = "blue", ...)
            dx <- ifelse(x < 0.5, 0.1, - 0.1)
            panel.xyplot(x + dx, y, pch = as.character(groups[subscripts]), col = "black", cex = 1.5, ...)
        }

        lp <-
            with(df,
                 xyplot(rev(seq_along(ECDF)) ~ ECDF | rank.segment, groups = pch, layout = c(1, 2),
                        panel = panel.rankplot,
                        xlim = c(0, 1),
                        scales = list(y = list(relation = "free", draw = FALSE)),
                        #col = "blue",
                        ylab = NULL, xlab = "ECDF at s = 0.9"))
        dimnames(lp)$rank.segment <- c(paste("top", pct.top.g, "%"),
                                       paste("bottom", 100 - pct.top.g, "%"))
        return(lp)
    }

plot.all <- function(plots) {
    print(plots$density, position = c(0.0, 0.85, 0.8, 1.0), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$ecdf, position = c(0.0, 0.7, 0.8, 0.85), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$level.sel.g, position = c(0.0, 0.55, 0.8, 0.7), more = TRUE)
    print(plots$rank, position = c(0.8, 0.0, 1.0, 0.55), more = TRUE)
    print(plots$level, position = c(0.0, 0.0, 0.8, 0.55), more = FALSE)
}


panel.foo <- function(x, y, groups, my.col, ...) {
    panel.barchart(x, y, groups = groups, col = my.col, ...)
}
(ecdf.bar <- barchart(gene ~ ds, data = df.long, groups = obs.rank, stack = TRUE, panel = panel.foo, my.col = df.long$col))
