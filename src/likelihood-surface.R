# log-likelihood of logistic model
#
# Arguments
# beta: regression coefficients
# y: count of success (higher readcount), possibly vector
# n: denominator (total readcount), possibly vector
# X: design matrix (as if made by model.matrix)
# do.neg: take the negative of logL or not
#
# Value
# The log-likelihood
ll.logi <- function(beta, y, n, X, do.neg = FALSE) {
    eta <- X %*% beta
    sgn <- ifelse(do.neg, -1, 1)
    sgn * sum(y * eta - n * log(1 + exp(eta)) + log(choose(n, y)))
}

# prepare arguments from model 'm' to be passed to ll.logi
args.logi <- function(m, do.neg = FALSE)
    list(beta = m$coefficients,
         y = m$model[["Y"]][ , 1],
         n = apply(m$model[["Y"]], 1, sum),
         X = model.matrix(m$formula, m$model),
         do.neg = do.neg)


# log-likelihood of weighted or unweighted normal linear model
#
# Arguments
# beta: regression coefficients
# sigma2: error variance
# y: count of success (higher readcount), possibly vector
# X: design matrix (as if made by model.matrix)
# W: a vector of weights
# do.neg: take the negative of logL or not
#
# Value
# The log-likelihood
ll.wnlm <- function(beta, sigma2, y, W, X, do.neg = FALSE) {
    eta <- X %*% beta
    sgn <- ifelse(do.neg, -1, 1)
    - sgn * sum((y - eta)^2 * W) / sigma2 / 2
}

# prepare arguments from model 'm' to be passed to ll.wnlm
args.wnlm <- function(m, do.neg = FALSE)
    list(beta = m$coefficients,
         # S^2, the unbiased estimate of sigma^2; I am not sure if it is
         # correct in the general weighted case
         sigma2 = sum((m$residuals)^2) / m$df.residual,
         y = m$model[["Y"]],
         W = m$weights,
         X = model.matrix(m$formula, m$model),
         do.neg = do.neg)


# create logL surface
#
# Arguments
# l.M: a list of model object; each component reflects a gene
# gene: the model to be selected from 'l.M' based on this gene
# n.pnts: the 2D parameter domain of logL is represented by an n.pnts x n.pnts grid
# v.name.A, v.name.B: the name of the two selected coefficients
# CI.lev.A, CI.lev.B: define the boundaries of the rectangular 2D parameter domain
# ll.fun: the function to compute logL
# args.fun: the function to prepare arguments to ll.fun
# do.neg: take the negative of logL or not
#
# Value
# A long, thin data frame that represents the grid's x and y coordinate with
# its 'beta.A' and 'beta.B' components.  The data frame contains a 'log.L' and
# a 'rel.log.L' component to be used as z coordinate.  Other components are
# useful for annotation and conditioning in lattice plots.  See
# 'll.surfacepolot' for such plotting function.
ll.grid <- function(l.M = M$logi.S, gene = "PEG3", n.pnts = 101,
                    v.name.A = "InstitutionPenn", v.name.B = "Age", CI.lev.A = 0.99, CI.lev.B = 0.99,
                    ll.fun = ll.logi, args.fun = args.logi, do.neg = FALSE) {
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
    df$rel.log.L <- df$log.L - do.call(ll.fun, args)
    #df$rel.log.L <- df$log.L - logLik(m)
    df$beta.hat.A <- args$beta[v.name.A] 
    df$beta.hat.B <- args$beta[v.name.B] 
    df$gene <- gene
    df$v.name.A <- v.name.A
    df$v.name.B <- v.name.B
    return(df)
}

# create levelplot from a logL surface produced by 'll.grid'
#
# Arguments
# fm: a formula to be passed to levelplot as first argument
# df: the logL surface to be passed to levelplot as the 'data' argument
# bits: to specify color bitdepth
# zeroline: whether to draw it
ll.surfaceplot <- function(fm = formula(- rel.log.L ~ beta.A * beta.B | gene), df = ll.grid(),
                           bits = 8, hv.A = NA, hv.B = NA, ...) {
    my.panel.abline <- function(hv, direction = "h") {
        l.B <- length(hv)
        l <- list(lty = c(1, rep(seq_len(l.B / 2) + 1, each = 2)), lwd = c(1, rep(2, length(hv))),
                  col = c("black", rep("white", length(hv))))
        if(direction == "h") l$h <- hv
        else l$v <- hv
        do.call(panel.abline, l)
    }
    levelplot(fm, data = df,
              scales = list(x = list(relation = "free")), colorkey = list(space = "top"),
              aspect = 1,
              b.hat.A = df$beta.hat.A, b.hat.B = df$beta.hat.B,
              panel = function(x, y, z, b.hat.A, b.hat.B, subscripts, ...) {
                  panel.levelplot(x = x, y = y, z = z, subscripts = subscripts, contour = FALSE, lwd = 0.5, ...)
                  panel.xyplot(x = b.hat.A[subscripts[1]], y = b.hat.B[subscripts[1]], col = "black", pch = 16)
                  my.panel.abline(hv.A, "v")
                  my.panel.abline(hv.B, "h")
                  panel.text(x = b.hat.A[subscripts[1]], y = b.hat.B[subscripts[1]], pos = 4, labels = expression(hat(beta)))
              },
              col.regions = hsv(h = 2^bits:1 / 2^bits, s = 0.9, v = 0.9),
              strip = strip.custom(bg = "white"),
              cuts = 2^bits,
              ...)
}


# wireframe-type 3D plot of log-likelihood surface
#
# Arguments
# df: the logL surface to be passed to levelplot as the 'data' argument
ll.wireframe <- function(dt = dat$coefs.wnlm.Q, v.A = "Ancestry.2", ...) {
    wireframe(rel.log.L ~ beta.A * beta.B | v.name.A, data = dt,
              subset = v.name.A == v.A, shade = TRUE, strip = FALSE,
              shade.colors.palette = function(irr, ref, height, ...)
                  trellis.par.get("shade.colors")$palette(irr, ref, 1 - height, ...),
              scales = list(arrows = FALSE),
              xlab = eval(substitute(expression(paste(beta, "[ ", v.A, " ]")), list(v.A = v.A))),
              ylab = expression(paste(beta, "[ Age ]")),
              zlab = "rel log L",
              ...)
}
