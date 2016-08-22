#const.term <- function(sigma2, y)
#    return(- length(y) / 2 * log10(2 * pi * sigma2))
#
#ll.nlm <- function(beta, sigma2, y, X, ll.const = 0) {
#    eps <- y - X %*% beta
#    ll.var <- - drop(t(eps) %*% eps) / (2 * sigma2)
#    return(ll.const + ll.var)
#}
#
#ll.args <- function(M, yname = "Y") {
#    w <- M$model[["(weights)"]]
#    w <- w / sum(w)
#    #W <- diag(w)
#    X <- diag(w) %*% model.matrix(M$formula, data = M$data)
#    y <- M$model[[ yname ]]
#    eps <- y - X %*% M$coefficients
#    sigma2 <- drop(t(eps) %*% eps) / length(eps) # MLE of residual error variance
#    list(sigma2 = sigma2,
#         y = y,
#         X = X,
#         ll.const = const.term(sigma2, y))
#}

ll.logi <- function(beta, y, n, X, do.neg = FALSE) {
    eta <- X %*% beta
    sgn <- ifelse(do.neg, -1, 1)
    sgn * sum(y * eta - n * log(1 + exp(eta)) + log(choose(n, y)))
}

args.logi <- function(m, do.neg = FALSE)
    list(beta = m$coefficients,
         y = m$model[["Y"]][ , 1],
         n = apply(m$model[["Y"]], 1, sum),
         X = model.matrix(m$formula, m$model),
         do.neg = do.neg)

ll.grid <- function(l.M = M$logi.S, gene = "PEG3", n.pnts = 101, v.name.A = "InstitutionPenn", v.name.B = "Age", CI.lev.A = 0.99, CI.lev.B = 0.99, ll.fun = ll.logi, args.fun = args.logi, do.neg = FALSE) {
    m <- l.M[[gene]]
    args <- args.fun(m, do.neg = do.neg)
    foo <- function(beta.A, beta.B) {
        beta <- args$beta
        beta[v.name.A] <- beta.A
        beta[v.name.B] <- beta.B
        do.call(ll.fun, c(beta = list(beta), args[-1]))
    }
    CI.A <- confint(m, level = CI.lev.A)
    CI.B <- confint(m, level = CI.lev.B)
    b.A <- seq(CI.A[v.name.A, 1], CI.A[v.name.A, 2], length.out = n.pnts)
    b.B <- seq(CI.B[v.name.B, 1], CI.B[v.name.B, 2], length.out = n.pnts)
    df <- expand.grid(b.A, b.B)
    names(df) <- c("beta.A", "beta.B")
    df$log.L <- sapply(as.data.frame(t(as.matrix(df))), function(x) foo(x[1], x[2]))
    df$rel.log.L <- df$log.L - logLik(m)
    df$beta.hat.A <- args$beta[v.name.A] 
    df$beta.hat.B <- args$beta[v.name.B] 
    df$gene <- gene
    df$v.name.A <- v.name.A
    df$v.name.B <- v.name.B
    return(df)
}

ll.surfaceplot <- function(fm = formula(- rel.log.L ~ beta.A * beta.B | gene), df = ll.grid(), ...) {
    levelplot(fm, data = df,
              scales = list(relation = "free"), colorkey = list(space = "top"),
              aspect = 1, layout = c(3, 1),
              b.hat.A = df$beta.hat.A, b.hat.B = df$beta.hat.B,
              panel = function(x, y, z, b.hat.A, b.hat.B, subscripts, ...) {
                  panel.levelplot(x = x, y = y, z = z, subscripts = subscripts, contour = TRUE, ...)
                  panel.xyplot(x = b.hat.A[subscripts[1]], y = b.hat.B[subscripts[1]], col = "black", pch = 16)
                  panel.text(x = b.hat.A[subscripts[1]], y = b.hat.B[subscripts[1]], pos = 4, labels = expression(hat(beta)))
              },
              ...)
}
