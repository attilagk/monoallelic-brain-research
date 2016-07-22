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

# separate plot strategy
# density for selected genes
plots$density <-
    densityplot(~ PEG10 + ZNF331, data = Y, plot.points = FALSE, type = c("p", "g"),
                scales = list(x = list(draw = FALSE), y = list(rot = 90, relation = "free")),
                xlim = c(0.5, 1),
                auto.key = list(corner = c(0.1, 0.9), lines = TRUE, points = FALSE),
                xlab = NULL, ylab = "density")

panel.my.ecdfplot <- function(x, groups, subscripts, ...){
    xx <- x[groups]
    yy <- ecdf(xx)(xx)
    #y1 <- ecdf(x[groups[subscripts]])(0.9)
    #panel.xyplot(xx, yy, ...)
    panel.densityplot(xx, ...)
}

# ECDF for selected genes
plots$ecdf <-
    ecdfplot(~ PEG10 + ZNF331, data = Y, type = c("s", "g"),
             panel = panel.ecdfplot,
             scales = list(x = list(draw = FALSE), y = list(rot = 90, relation = "free")),
             xlim = c(0.5, 1),
             ylab = "empirical CDF",
             xlab = NULL)

# ECDF genome-wide
# data manipulation
ED.long$rank.segment <- factor(rep(1, nrow(ED.long)), levels = c("top", "bottom"), ordered = TRUE)
ED.long$rank.segment[ with(ED.long, gene %in% levels(gene)[1:100]) ] <- "top"
ED.long$rank.segment[ with(ED.long, ! gene %in% levels(gene)[1:100]) ] <- "bottom"
#n.all.g <- length(ok.genes)
n.all.g <- 1000
pct.top.g <- 5 # percent 
n.top.g <- ceiling(n.all.g * pct.top.g / 100)
plots$level <- levelplot(ECDF ~ s * gene | rank.segment, data = ED.long, legend = NULL, colorkey = FALSE,
                         scales =
                             list(y =
                                  list(rot = 90, relation = "free",
                                       limits = list(c(n.top.g, 1),
                                                     c(n.all.g, n.top.g + 1)))),
                         xlab = "s, imbalance score", ylab = "gene rank",
                         layout = c(1, 2))
dimnames(plots$level)$rank.segment <- c(paste("top", pct.top.g, "% of genes"),
                                        paste("bottom", 100 - pct.top.g, "% of genes"))

if(FALSE) {
    print(plots$density, position = c(0.0, 0.85, 1.0, 1.0), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$ecdf, position = c(0.0, 0.7, 1.0, 0.85), panel.height = list(0.9, "npc"), more = TRUE)
    print(plots$level, position = c(0.0, 0.0, 1.0, 0.7), more = FALSE)
}

