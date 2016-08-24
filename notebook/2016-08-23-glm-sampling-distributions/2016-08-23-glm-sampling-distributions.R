logistic.f <- function(eta)
    1 / (1 + exp(-eta))

eta.f <- function(x, beta = c(0, -1))
    cbind(1, x) %*% beta

pdf.logi <- function(x, y, beta = c(0, -1), n = 50) {
    p <- logistic.f(eta.f(x, beta = beta))
    dbinom(round(y * n), size = n, prob = p, log = TRUE)
}

grid.density <- function(xlim = c(-10, 10), ylim = c(0, 1), n.pnts = 201, pdfun = pdf.logi, beta = c(0, -1), ...) {
    x <- seq(xlim[1], xlim[2], length.out = n.pnts)
    y <- seq(ylim[1], ylim[2], length.out = n.pnts)
    df <- expand.grid(x, y)
    names(df) <- c("x", "y")
    df$density <- sapply(as.data.frame(t(as.matrix(df))), function(x) pdfun(x[1], x[2], beta = beta, ...))
    return(df)
}

plot.density <- function(beta = c(0, -1), xlim = c(-30, 30), ...) {
    gd = grid.density(xlim = xlim, beta = beta, ...)
    levelplot(density ~ x * y, data = gd, beta = beta,
              panel = function(x, ...) {
                  panel.levelplot(x = x, ...)
                  panel.xyplot(x = x, y = logistic.f(eta.f(x, beta = beta)), type = "l", col = "black")
              },
              par.settings = list(regions = list(col = rev(trellis.par.get("regions")$col))))
}
