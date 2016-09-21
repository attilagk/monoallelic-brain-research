get.estimate.CI.permut <- function(type = "logi.S", M, e.v, conf.lev = 0.99) {
    betas <- lapply(e.v,
                    function(e) {
                        x <- get.estimate.CI(M[[e]][[type]], conf.lev = conf.lev)
                        x <- x[ grep(e, row.names(x), value = TRUE), ]
                    })
    betas$RIN2 <- NULL #  regex "RIN" already matched both RIN and RIN2
    return(do.call(rbind, betas))
}

is.signif <- function(betas, sel.coef = levels(betas$Coefficient)) {
    betas <- betas[levels(betas$Coefficient) %in% sel.coef, ]
    betas$Lower.CL > 0 | betas$Upper.CL < 0 
}

#before.after <- function(m1 = "logi.S", m2 = "wnlm.R",
#                         betas1 = list(U = Betas.Unpermuted, P = Betas.Permuted),
#                         betas2 = list(U = Betas.Unpermuted, P = Betas.Permuted),
#                         sel.coef = levels(betas1$U[[m1]]$Coefficient)
#                         ) {
#    list(U = sum(is.signif(betas1$U[[m1]], sel.coef) & is.signif(betas2$U[[m2]], sel.coef), na.rm = TRUE),
#         P = sum(is.signif(betas1$P[[m1]], sel.coef) & is.signif(betas2$P[[m2]], sel.coef), na.rm = TRUE))
#}
