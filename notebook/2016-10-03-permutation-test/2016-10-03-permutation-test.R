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
