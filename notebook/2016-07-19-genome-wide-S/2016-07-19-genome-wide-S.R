emp.distr.S <-
    function(Y,
             ss = seq(0.5, 1, length.out = 101),
             with.density = FALSE,
             ...) {
        E <- list()
        E$ecdf <- lapply(Y, ecdf)
        E$ecdf.val <- data.frame(lapply(E$ecdf, function(f) f(ss)))
        E$ecdf.val$ss <- ss
        if(with.density) {
            E$density <-
                data.frame(lapply(Y, function(x)
                                  density(x, na.rm = TRUE, n = length(ss), from = 0.5, to = 1, ...)$y))
            E$density$ss <- ss
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
            ...)
    }

# complex figure on the distribution of S_ig

#print(dp, split = c(1, 1, 1, 4), more = TRUE)
#print(ep, split = c(1, 2, 1, 4), more = TRUE)
#print(lp.100, split = c(1, 3, 1, 4), more = TRUE)
#print(lp, split = c(1, 4, 1, 4), more = FALSE)

#dp <- densityplot(~ PEG10 + ZNF331 + ZNF286A + MAPK9, data = Y, plot.points = FALSE, xlab = FALSE)
#ep <- with(ED$ecdf.val, xyplot(PEG10 + ZNF331 + ZNF286A + MAPK9 ~ ss, type = c("g", "s"), xlab = NULL, ylab = NULL))
lp <- levelplot(as.matrix(ED[[2]][rev(gene.order)]), aspect = "fill", par.settings = tp, scales = sc)
lp.100 <- levelplot(as.matrix(ED[[2]][rev(gene.order[1:100])]), aspect = "fill", par.settings = tp, scales = sc)
