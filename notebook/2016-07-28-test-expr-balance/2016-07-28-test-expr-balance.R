CI.p <- function(p.hat, n, conf.lev = 0.95) {
    alpha <- (1 - conf.lev) / 2
    zz <- lapply(list(lower = alpha, upper = 1 - alpha), qnorm)
    std <- sqrt(p.hat * (1 - p.hat) / n)
    lapply(zz, function(z) p.hat + z * std)
}
