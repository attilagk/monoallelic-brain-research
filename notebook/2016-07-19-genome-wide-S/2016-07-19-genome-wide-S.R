# Empirical distribution of S
#
# Arguments
# Y: a data frame
#
# Value
#
# Details
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

#
# Arguments
#
# Value
#
# Details
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
                #par.settings = list(superpose.line = superpose.line),
                auto.key = list(text = paste0(letters[seq_along(sel.g)], ": ", sel.g),
                                corner = c(0.1, 0.9), lines = TRUE, points = FALSE),
                ...)
}

# new ECDF plot
my.ecdfplot <- function(Y.long, eval.at = c(0.9), survival = TRUE, ...) {
    df <- expand.grid(x = c(0.5, 0.55), y = 0:100 / 100)
    df$z <- rep(0:100 / 100, each = 2)
    df$z[seq(2, 202, by = 2)] <- NA
    panel.my.ecdfplot <- function(x, groups, subscripts, eval.at = eval.at, df, ...){
        with(df, panel.levelplot(x = x, y = y, z = z, at = 0:100 / 100, subscripts = 0:202, ...))
        panel.superpose(x, groups = groups, subscripts = subscripts,
                        panel.groups = panel.ecdfplot, type = c("s"), ...)
        panel.grid(h = 0, v = -1)
        g.lev <- levels(groups)
        ss <- rep(c(eval.at,
                    rep(NA, length(x) / length(g.lev) - length(eval.at))),
                  length(g.lev))
        yy <- sapply(levels(groups), function(g) ecdf(x[groups == g])(ss))
        gg <- sapply(levels(groups), function(g) ecdf(x[groups == g])(ss))
        panel.superpose(x = ss, y = yy, groups = groups, subscripts = subscripts,
                        panel.groups = panel.xyplot, col = "blue", pch = 21, ...)
        lapply(levels(groups),
               function(g) panel.xyplot(x = eval.at, y = ecdf(x[groups == g])(eval.at),
                                        ...))
        lapply(seq_along(levels(groups)),
               function(i)
                   panel.xyplot(x = eval.at - 0.00, y = ecdf(x[groups == levels(groups)[i]])(eval.at),
                                pch = letters[i], col = "black", cex = 1.5, ...))
    }

    ep <- ecdfplot(~ s, data = Y.long, groups = gene, eval.at = c(0.9),
                   df = df,
                   scales = list(x = list(draw = FALSE),
                                 y = list(rot = 90, relation = "free", at = yloc <- seq(0, 1, length.out = 5), label = yloc)),
                   panel = panel.my.ecdfplot.b,
                   xlim = c(0.5, 1),
                   ylab = "ECDF", xlab = NULL, ...)
    if(survival)
        ep <- update(ep, ylab = "1 - ECDF",
                     scales = list(x = list(limits = c(0.5, 1)), y = list(limits = c(1, 0), labels = rev(yloc))))
    return(ep)
}


# even newer ECDF plot; this colors curves with a ramp according to the
# cumulative probability (y axis)
my.ecdfplot2 <- function(Y.long, invert.col.grad = TRUE, ramp.cols = c("red", "blue"), eval.at = c(0.9), ...) {
    genes <- levels(Y.long$gene)
    names(genes) <- genes
    helper <- function(i) {
        g <- genes[i]
        # jump locations (knots) and y values (probabilities) of the ECDF for S distribution
        k <- knots(ecdf.g <- with(Y.long, ecdf(s[gene == g])))
        k <- c(0.5, k)
        y <- sapply(k, ecdf.g)
        # coloring vector
        color <- colorRampPalette(ramp.cols)(length(k))
        names(color) <- color
        df <- data.frame(y = y, k = k, color = color, stringsAsFactors = FALSE)
        # small 2-point segments; one for each group
        d <- data.frame(lapply(df, rep, each = 2), stringsAsFactors = FALSE)
        # without setting row names R warns
        row.names(d) <- seq.int(nr <- nrow(d))
        # shift k values relative to y values; necessary for
        # correct display of the 2-point segments
        d$k <- c(d$k[-1], 1)
        d <- cbind(d, data.frame(gene = g), stringsAsFactors = FALSE)
        #d <- cbind(d[c(seq_len(nr - 1) + 1, nr), c("y", "k")], d["color"], data.frame(gene = g), stringsAsFactors = FALSE)
        invert.fun <- if(invert.col.grad) rev else I
        d$group <- paste(invert.fun(d$color), d$gene, sep = ".")
        d$eval.at.y <- c(ecdf.g(eval.at), rep(NA, nr - 1))
        d$eval.at.x <- c(eval.at, rep(NA, nr - 1))
        d$letter <- c(letters[i], rep(NA, nr - 1))
        row.names(d) <- seq.int(nr)
        return(d)
    }
    #return(lapply(genes.ix, helper))
    df <- data.frame(do.call(rbind, lapply(seq_along(genes), helper)))
    #return(df)
    panel.my.ecdfplot2 <- function(x, y, groups, subscripts, mycol, eval.at.x, eval.at.y, letter, ...) {
        panel.grid(h = 0, v = -1)
        panel.superpose(x = x, y = y, groups = groups, subscripts = subscripts, col = mycol, type = "l", lty = "solid", ...)
        panel.xyplot(x = eval.at.x, y = 1 - eval.at.y, type = "p", ...)
        panel.xyplot(x = eval.at.x, y = 1 - eval.at.y, type = "p", pch = letter, col = "black", cex = 1.5, ...)
    }
    scl <- list(x = list(draw = FALSE),
                y = list(rot = 90, relation = "free", at = yloc <- seq(0, 1, length.out = 5), label = yloc))
    xyplot(1 - y ~ k, data = df, groups = df$group, panel = panel.my.ecdfplot2, mycol = df$color,
           eval.at.x = df$eval.at.x, eval.at.y = df$eval.at.y, letter = df$letter,
           xlim = c(0.5, 1), ylab = "1 - ECDF", xlab = NULL, scales = scl, ...)
}

panel.my.levelplot <- function(..., grid.h, grid.v) {
    panel.levelplot(...)
    panel.grid(h = grid.h, v = grid.v)
}

my.levelplot.sel.g <- function(ED.long, sel.g = c("PEG10", "ZNF331", "AFAP1"), letter.labels = FALSE, ...) {
    ED.long$gene <- with(ED.long, factor(gene, levels = rev(levels(gene)), ordered = TRUE))
    y.rot <- ifelse(letter.labels, 90, 0)
    y.labels <- ifelse(rep(letter.labels, length(sel.g)), letters[rev(seq_along(sel.g))], rev(sel.g))
    levelplot(ECDF ~ s * gene, data = ED.long, legend = NULL, colorkey = FALSE,
                    panel = panel.my.levelplot, subset = gene %in% rev(sel.g),
                    xlab = NULL,
                    grid.h = 0, grid.v = -1,
                    scales =
                        list(x = list(draw = FALSE),
                             y = list(rot = y.rot, labels = y.labels, relation = "free")),
                    ...)
}

# ECDF genome-wide
my.levelplot <- function(ED.long, pct.top.g = 2, n.all.g = length(levels(ED.long$gene)), top.on.top = TRUE, ...) {
    # data manipulation
    n.top.g <- ceiling(n.all.g * pct.top.g / 100)
    ED.long$rank.segment <- factor(rep(1, nrow(ED.long)), levels = c("top", "bottom"), ordered = TRUE)
    ED.long$rank.segment[ with(ED.long, gene %in% levels(gene)[seq_len(n.top.g)]) ] <- "top"
    ED.long$rank.segment[ with(ED.long, ! gene %in% levels(gene)[seq_len(n.top.g)]) ] <- "bottom"
    if(top.on.top) ED.long <- top.on.top(ED.long)
    # create trellis object
    n.grid.v <- 3
    lp <- levelplot(ECDF ~ s * gene | rank.segment, data = ED.long, legend = NULL, colorkey = FALSE,
                    panel = panel.my.levelplot,
                    grid.h = - n.grid.v, grid.v = -1,
                    scales =
                        list(y =
                             list(rot = 90, relation = "free",
                                  at = list(pretty(seq(n.top.g, 1), n = n.grid.v),
                                            pretty(seq(n.all.g, 1), n = n.grid.v)),
                                  limits = list(c(n.top.g, 1),
                                                c(n.all.g, n.top.g + 1)))),
                    layout = c(1, 2),
                    ...)
    # edit dimnames
    if(top.on.top)
        bottom <- "all genes"
    else
        bottom <- paste("bottom", 100 - pct.top.g, "% of", "genes")
    dimnames(lp)$rank.segment <- c(paste("top", pct.top.g, "% of genes"), bottom)
    return(lp)
}

rankplot.2 <- function(frac, pct.top.g = 2, sel.g = c("PEG10", "ZNF331", "AFAP1"), top.on.top = TRUE, ...) {
    n.top.g <- ceiling(length(frac) * pct.top.g / 100)
    df <- data.frame(score = unlist(frac["1", ]))
    df$gene <- factor(names(frac), levels = names(frac), ordered = TRUE)
    df$pch <- rep("", nrow(df))
    df$pch[with(df, gene %in% sel.g)] <- letters[seq_along(sel.g)]
    df$rank.segment <-
        factor(c(rep("top", n.top.g),
                 rep("bottom", length(frac) - n.top.g)),
               levels = c("top", "bottom"), ordered = TRUE)
    df$score.sel.g <- df$score
    df$score.sel.g[! with(df, gene %in% sel.g)] <- NA
    if(top.on.top) df <- top.on.top(df)
    my.prepanel <- function(x, y, ...) {
        list(ylim = rev(range(y)))
    }

    panel.rankplot <- function(x, y, groups, subscripts, gene, ...) {
        panel.grid(h = -3, v = 0)
        panel.superpose(x, y, groups = groups, subscripts = subscripts, ...)
        panel.xyplot(x, y, pch = as.character(gene[subscripts]), col = "black", cex = 1.5, ...)
    }

    rp <- xyplot(seq_along(score) ~ score + score.sel.g | rank.segment, data = df,
                 x.sel.g = unlist(frac["1", sel.g]),
                 y.sel.g = match(sel.g, names(frac)),
                 gene = df$pch,
                 layout = c(1, 2), prepanel = my.prepanel, panel = panel.rankplot,
                 scales = list(y = list(relation = "free", draw = FALSE)),
                 ...)
    # edit dimnames
    if(top.on.top)
        bottom <- "all genes"
    else
        bottom <- paste("bottom", 100 - pct.top.g, "%")
    dimnames(rp)$rank.segment <- c(paste("top", pct.top.g, "%"), bottom)
    return(rp)
}

rankplot <-
    function(cum.freq, sorted.genes = names(cum.freq),
             plot.andys.test = FALSE, pct.top.g = 2, sel.g = c("PEG10", "ZNF331", "AFAP1"),
             top.on.top = TRUE,
             ...) {
        n.top.g <- ceiling(length(sorted.genes) * pct.top.g / 100)
        df <- data.frame(ECDF = unlist(cum.freq["0.9", ]))
        df$andys.test = unlist(cum.freq["andys.test", ])
        df$gene <- factor(sorted.genes, levels = c(rev(sorted.genes)), ordered = TRUE)
        df$pch <- rep("", nrow(df))
        df$pch[with(df, gene %in% sel.g)] <- letters[seq_along(sel.g)]
        df$rank.segment <-
            factor(c(rep("top", n.top.g),
                     rep("bottom", length(sorted.genes) - n.top.g)),
                   levels = c("top", "bottom"), ordered = TRUE)
        df$ECDF.sel.g <- df$ECDF
        df$ECDF.sel.g[! with(df, gene %in% sel.g)] <- NA
        if(top.on.top) df <- top.on.top(df)

        panel.rankplot <- function(x, y, groups, subscripts, ...) {
            panel.superpose(x, y, groups = groups, subscripts = subscripts, col = "blue", ...)
            dx <- ifelse(x < 0.5, 0.1, - 0.1)
            panel.xyplot(x + dx, y, pch = as.character(groups[subscripts]), col = "black", cex = 1.5, ...)
        }

        if(plot.andys.test)
            fm <- formula(rev(seq_along(ECDF)) ~ andys.test + ECDF + ECDF.sel.g | rank.segment)
        else
            fm <- formula(rev(seq_along(ECDF)) ~ ECDF + ECDF.sel.g | rank.segment)
        lp <-
            xyplot(fm, layout = c(1, 2),
                   data = df, xlim = c(0, 1), scales = list(y = list(relation = "free", draw = FALSE)),
                   ylab = NULL, xlab = "ECDF at s = 0.9",
                   ...)
        # edit dimnames
        if(top.on.top)
            bottom <- "all genes"
        else
            bottom <- paste("bottom", 100 - pct.top.g, "%")
        dimnames(lp)$rank.segment <- c(paste("top", pct.top.g, "%"), bottom)
        return(lp)
    }

# place the top portion of a data frame on the entire data frame
#
top.on.top <- function(df, tb.var = "rank.segment", tb.str = c("top", "bottom", "all")) {
    top <- df[df[[tb.var]] == tb.str[1], ]
    levels(df[[tb.var]]) <- tb.str
    df[[tb.var]] <- tb.str[3] # the entire original data frame becomes the new "bottom" labeled as all
    return(rbind(top, df))
}

plot.all <- function(plots) {
    print(plots$density, position = c(0.0, 0.85, 0.8, 1.0), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$ecdf, position = c(0.0, 0.7, 0.8, 0.85), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$level.sel.g, position = c(0.0, 0.55, 0.8, 0.7), more = TRUE)
    print(plots$rank, position = c(0.8, 0.0, 1.0, 0.55), more = TRUE)
    print(plots$level, position = c(0.0, 0.0, 0.8, 0.55), more = FALSE)
}


plot.all.b <- function(plots) {
    print(plots$density, position = c(0.0, 0.85, 0.8, 1.0), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$ecdf, position = c(0.0, 0.7, 0.8, 0.85), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$level.sel.g, position = c(0.0, 0.58, 0.8, 0.73), more = TRUE)
    print(plots$rank, position = c(0.8, 0.0, 1.0, 0.65), more = TRUE)
    print(plots$level, position = c(0.0, 0.0, 0.8, 0.65), more = FALSE)
}


plot.all.c <- function(plots) {
    print(plots$ecdf, position = c(0.0, 0.75, 0.8, 1.0), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$level.sel.g, position = c(0.0, 0.6, 0.8, 0.75), more = TRUE)
    print(plots$rank, position = c(0.8, 0.0, 1.0, 0.65), more = TRUE)
    print(plots$level, position = c(0.0, 0.0, 0.8, 0.65), more = FALSE)
}


# Prepare data Y to ECDF bar plot
ecdf.bar.prep <- function(Y) {
    E <- lapply(Y, ecdf)
    n.tot <- nrow(Y) # sample size including NAs
    # knots: the sets consisting of 0.5 and the jump points (discontinuities) of the ECDF
    # the m-length vector of knots must be padded to length n.tot + 1
    K <- data.frame(lapply(E,
                           function(f) c(k <- c(0.5, knots(f)), rep(1, n.tot - length(k) + 1))))
    # the lengths of continuous intervals of the ECDF of S within the interval [0.5, 1]
    # a vector of length n.tot
    L <- data.frame(lapply(K, diff))
    sel.g <- names(Y)
    # ECDF values at K except for the last value (which must be one); also of length n.tot
    Fs <- data.frame(lapply(sel.g, function(g) E[[g]](K[[g]][seq_len(n.tot)])))
    # color coding
    #col <- data.frame(lapply(Fs, function(x) trellis.par.get("regions")$col[cut(x, 100, labels = FALSE)]))
    col <- data.frame(lapply(Fs, function(x) trellis.par.get("regions")$col[cut(x, 100, labels = FALSE)]),
                      stringsAsFactors = FALSE)
    names(K) <- paste("K", sel.g, sep = ".")
    names(L) <- paste("L", sel.g, sep = ".")
    names(Fs) <- paste("Fs", sel.g, sep = ".")
    names(col) <- paste("col", sel.g, sep = ".")
    # rank of observation (including NAs)
    R <- data.frame(R = factor(seq_len(n.tot)))
    Z <- reshape(cbind(K[seq_len(n.tot), ], L, Fs, col, R), direction = "long", v.names = c("K", "L", "Fs", "col"),
                 varying = list(paste0("K.", sel.g), paste0("L.", sel.g), paste0("Fs.", sel.g), paste0("col.", sel.g)),
                 timevar = "gene", times = sel.g, idvar = "R", ids = R)
    Z$gene <- factor(Z$gene, levels = sel.g[order(sapply(E, function(f) f(0.9)))], ordered = TRUE)
    return(Z)
}

ecdf.barchart <- function(Y.long = ecdf.bar.prep(Y[ok.genes[gene.order[1:100]]]), ...) {
    panel.ecdf.barchart <- function(x, y, subscripts, my.col, ...) {
        panel.barchart(x, y, subscripts = subscripts, col = my.col[subscripts], ...)
        panel.grid(h = 0, v = -1)
    }

    with(Y.long,
         barchart(~ L | gene, groups = R, panel = panel.ecdf.barchart, my.col = col,
                  stack = TRUE, layout = c(1, length(levels(gene))), as.table = TRUE, box.width = 1.5,
                  strip = FALSE,
                  par.settings =
                      list(axis.line = list(col = "transparent"),
                           superpose.polygon = list(lty = 0)
                           ),
                  between = list(y = 0)), ...)
}
