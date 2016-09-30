is.signif <- function(betas, sel.coef = levels(betas$Coefficient)) {
    betas <- betas[levels(betas$Coefficient) %in% sel.coef, ]
    betas$Lower.CL > 0 | betas$Upper.CL < 0 
}
