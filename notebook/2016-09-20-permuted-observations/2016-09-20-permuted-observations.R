is.signif <- function(betas, sel.coef = levels(betas$Coefficient)) {
    betas <- betas[levels(betas$Coefficient) %in% sel.coef, ]
    betas$Lower.CL > 0 | betas$Upper.CL < 0 
}

my.segplot2 <- function(coef, type, gene.ix, main = type, ...) {
    segplot(ordered(Permutation, levels = rev(levels(Permutation))) ~ Lower.CL + Upper.CL | Gene, data = Betas,
            subset = Coefficient %in% coef & Model %in% type,
            level = factor(! Permutation == "U"),
            colorkey = FALSE,
            panel = function(x, y, ...) {
                panel.segplot(x, y, ...)
                panel.abline(v = 0, col = "blue")
            },
            xlab = eval(substitute(expression(paste(beta, "[ ", coef, " ]")), list(coef = coef))),
            ylab = "permutation", main = main,
            scales = list(draw = FALSE, x = list(relation = "free")))[gene.ix]
}

#before.after <- function(m1 = "logi.S", m2 = "wnlm.R",
#                         betas1 = list(U = Betas.Unpermuted, P = Betas.Permuted),
#                         betas2 = list(U = Betas.Unpermuted, P = Betas.Permuted),
#                         sel.coef = levels(betas1$U[[m1]]$Coefficient)
#                         ) {
#    list(U = sum(is.signif(betas1$U[[m1]], sel.coef) & is.signif(betas2$U[[m2]], sel.coef), na.rm = TRUE),
#         P = sum(is.signif(betas1$P[[m1]], sel.coef) & is.signif(betas2$P[[m2]], sel.coef), na.rm = TRUE))
#}

## extract p-values of every coefficient for a given gene, and model type
#foo <- function(gene, M, type = "logi.S", e.vars) {
#    helper <- function(v) {
#                      x <- summary(M[[c(v, type, gene)]])$coefficients
#                      x[ grep(v, row.names(x)), 4]
#    }
#    data.frame(pvalue = unlist(sapply(e.vars[! e.vars %in% "RIN2"], helper)),
#               coef = grep("Intercept", names(coef(M[[c(1, 1, 1)]])), value = TRUE, invert = TRUE),
#               gene = gene)
#}
#
## extract p-values of every gene for a given coefficient, and model type
#foo2 <- function(v, M, type = "logi.S", e.vars) {
#    helper <- function(gene) {
#                      x <- summary(M[[c(v, type, gene)]])$coefficients
#                      x[ grep(v, row.names(x)), 4]
#    }
#    data.frame(pvalue = unlist(sapply(e.vars[! e.vars %in% "RIN2"], helper)),
#               coef = grep("Intercept", names(coef(M[[c(1, 1, 1)]])), value = TRUE, invert = TRUE),
#               gene = gene)
#}
#
#helper <- function(gene, v, type) {
#    x <- summary(M[[c(v, type, gene)]])$coefficients
#    x[ grep(v, row.names(x)), 4]
#}
