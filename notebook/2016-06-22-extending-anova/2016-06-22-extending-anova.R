mtype.compare.plot <- function(mtypeA = "logi.S", mtypeB = "wnlm.Q", dt = do.call(cbind, Betas), do.key = TRUE, fit.OK = logi.S.OK$logi.S.fit.OK, ...) {
    mycol <- rainbow(n = length(levels(dt$wnlm.Q.Gene)), start = 0, end = 5 / 6)
    mypch <- c(23, 22)
    key <- 
        if(do.key)
            list(columns = 4, points = list(pch = ifelse(fit.OK, mypch[1], mypch[2]), fill = mycol, alpha = 0.5),
                 text = list(paste0("(", seq_along(levels(Betas[[1]]$Gene)), ") ", levels(Betas[[1]]$Gene)),
                             cex = 0.7, font = ifelse(fit.OK, 2, 1)))
        else
            FALSE
    fm <- formula(paste0(mtypeA, ".Estimate ~ ", mtypeB, ".Estimate | wnlm.Q.Coefficient"))
    xyplot(fm, data = dt,
           groups = dt$wnlm.Q.Gene,
           gene = as.character(as.integer(dt$wnlm.Q.Gene)),
           par.settings = list(superpose.symbol = list(fill = mycol), add.text = list(cex = 1.0)),
           scales = list(relation = "free", draw = FALSE),
           panel = function(x, y, subscripts, gene, fit.OK,...){
               panel.abline(h = 0, lty = 2)
               panel.abline(v = 0, lty = 2)
               col <- trellis.par.get("superpose.symbol")$fill
               panel.xyplot(x, y, pch = ifelse(fit.OK, mypch[1], mypch[2]), fill = col, col = col, cex = 1.5, alpha = 0.5)
               panel.text(x, y, labels = gene[subscripts], cex = 0.7, font = ifelse(fit.OK, 2, 1), ...)
           },
           main = expression(paste("Estimate ", hat(beta), "[...] under two models")),
           xlab = paste("under", mtypeB), ylab = paste("under", mtypeA),
           key = key,
           fit.OK = fit.OK,
           ...)[1:20]
}


reverse.genes <- function(df) {
    g <- df$Gene
    df$Gene <- factor(g, levels = rev(levels(g)))
    return(df)
}
