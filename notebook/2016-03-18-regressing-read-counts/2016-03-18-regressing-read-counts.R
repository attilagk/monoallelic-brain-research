logistic.upper <- function(x, b1, x0) {
    l <- 1 / (1 + exp(- (x - x0) * b1))
    l[l < 1/2] <- 1/2
    return(l)
}

# x is age and p_v is expected read fraction of alt. allele
plot(x, logistic.upper(x,-1 / 30, 150), type='l', xlab='x', ylab=expression(p[v]))
