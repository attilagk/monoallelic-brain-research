# To check consistency with model logLik values that R provides
my.logLik <- function(SS, n) {
    - n / 2 * (log(2 * pi) + log(SS) - log(n) + 1)
}

aic.shift <- function(c, n = 491)
    n * (1 + log(c) - c)


ci4logifun <- function(logifun = function(x) plogis(-x, scale=1), delim = 300, pr = 0.05) {
    cifun <- list(lower=NA, upper=NA)
    cifun$lower <- function(x) {
        q <- qbinom(logifun(x), size=delim, prob=0.05, lower.tail=TRUE)
        return(0.5 - q / delim)
        #return(logifun(x) - (0.5 - q / delim))
        #return((delim / 2 - q) / delim - 0.5 + logifun(x))
    }
    cifun$upper <- function(x) {
        q <- qbinom(logifun(x), size=delim, prob=0.05, lower.tail=FALSE)
        return((delim / 2 - q) / delim - 0.5 + logifun(x))
    }
    return(cifun)
}

# median(unlist(d[ paste0("N_", genes) ]), na.rm=TRUE) # should be 302

