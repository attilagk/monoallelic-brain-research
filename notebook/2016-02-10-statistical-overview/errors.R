get.error.rates <- function(pi1, thrs, lambda) {
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
