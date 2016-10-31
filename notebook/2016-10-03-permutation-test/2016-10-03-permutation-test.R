get.both.p.vals <- function(mtype = "wnlm.Q", gene = "MEST", coef = "Age", M, B) {
    errfun <- function(e) NA
    P <- B[ with(B, Model == mtype & Gene == gene & Coefficient == coef), c("Estimate", "Permutation")]
    data.frame(Model = mtype, Gene = gene, Coefficient = coef,
               Estimate = tryCatch(summary(M[[c(mtype, gene)]])$coefficients[coef , 1],
                                  error = errfun),
               p.val.t.dist = tryCatch(summary(M[[c(mtype, gene)]])$coefficients[coef , 4],
                                  error = errfun),
               p.val.perm = tryCatch(get.p.val(P[P$Permutation == "U", "Estimate"],
                                               P[P$Permutation != "U", "Estimate"]),
                                     error = errfun))
}


annotate.signif <- function(v) {
    sapply(v, function(x) {
               if(is.na(x)) return(NA)
               else {
                   if(x < 1e-3) return("***")
                   if(x < 1e-2) return("**")
                   if(x < 5e-2) return("*")
                   else return("")
               }})}


pvalplot.genes.as.panels <- function(bothpv = both.p.val, ...) {
    xyplot(- log10(p.val.t.dist) ~ - log10(p.val.perm) | Gene, data = bothpv, groups = Model,
           par.settings = list(add.text = list(cex = 0.8)),
           panel = function(...) {
               panel.abline(a = 0, b = 1, col = "gray", lty = "dotted")
               panel.xyplot(...)
           },
           auto.key = list(columns = 2),
           xlab = expression(paste(plain(-log)[10], "p  (permutations)")),
           ylab = expression(paste(plain(-log)[10], "p  (t-distribution)")),
           ...)
}
