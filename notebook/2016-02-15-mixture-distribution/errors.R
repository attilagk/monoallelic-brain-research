get.error.rates <- function(pi1, lambda, thrs) {
    # calculate probability of true positives, false negatives,...
    get.joint.probs <- function() {
        pfn <- (1 - exp(- lambda)) * (exp(- lambda * thrs) - exp(- lambda))
        list(
             ptp = pi1 * (1 - pfn),
             pfn = pi1 * pfn,
             pfp = thrs * (1 - pi1),
             ptn = (1 - thrs) * (1 - pi1)
             )
    }
    p <- get.joint.probs()
    precision <- p$ptp / (p$ptp + p$pfp)
    negpredval <- p$ptn / (p$ptn + p$pfn)
    list(
        ppv = precision,
        fdr = 1 - precision,
        npv = negpredval,
        fom = 1 - negpredval,   # false omission rate
        fpr = thrs,
        tpr = p$ptp / (p$ptp + p$pfn)   # sensitivity, recall
    )
}

# Creates a prob. density fun. for a mixture of two p-value densities: one uniform, one trimmed exponential.
# p1 is the prior probability of the exponential component, lamda is its rate constant.
dpval.maker <- function(pi1, lambda) {
    foo <- function(p)
        1 - pi1 + pi1 * lambda * (1 - exp(- lambda)) * exp(- lambda * p)
    class(foo) <- 'dpval'
    attributes(foo) <- list(pi1=pi1, lambda=lambda)
    return(foo)
}

# Samples n p-values from a mixture density parametrized by pi1 and lambda.
# See dpval.maker for creating pdf
pval.sampler <- function(dpval, n) {
    pi1 <- attributes(dpval)$pi1
    lambda <- attributes(dpval)$lambda
    # function for sampling m1 points from a truncated exponential distribution
    trunc.exp.sampler <- function(m1, s1) {
        l <- ((s <- rexp(m1, rate=lambda)) > 1)
        if(m2 <- length(which(l))) {
            trunc.exp.sampler(m2, c(s1, s[!l]))
        }
        else return(c(s1, s))
    }
    numexp <- rbinom(1, size=n, prob=pi1)
    s.unif <- runif(n=n-numexp)
    s.trunc.exp <- trunc.exp.sampler(numexp, c())
    c(s.trunc.exp, s.unif)
}
