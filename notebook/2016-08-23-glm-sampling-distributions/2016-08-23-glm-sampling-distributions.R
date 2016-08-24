logistic.f <- function(eta)
    1 / (1 + exp(-eta))

logistic.f2 <- function(eta)
    0.5 / (1 + exp(-eta)) + 0.5

eta.f <- function(x, beta = c(0, -1))
    cbind(1, x) %*% beta

pdf.nlm <- function(x, y, beta = c(0, -1), n = 50, sigma = 0.1) {
    mu <- eta.f(x, beta = beta)
    dnorm(y, mean = mu, sd = sigma, log = TRUE)
}

pdf.logistic <- function(x, y, beta = c(0, -1), n = 50, logi.f) {
    p <- logi.f(eta.f(x, beta = beta))
    dbinom(round(y * n), size = n, prob = p, log = TRUE)
}

pdf.logi <- function(x, y, beta = c(0, -1), n = 50)
    pdf.logistic(x = x, y = y, beta = beta, n = n, logi.f = logistic.f)

pdf.logi2 <- function(x, y, beta = c(0, -1), n = 50)
    pdf.logistic(x = x, y = y, beta = beta, n = n, logi.f = logistic.f2)

grid.density <- function(pdfun = pdf.logi, beta = c(0, -1), n = 50, sigma = 0.1, xlim = c(-50, 200), ylim = c(0.4, 1.1), n.pnts = 201, ...) {
    x <- seq(xlim[1], xlim[2], length.out = n.pnts)
    y <- seq(ylim[1], ylim[2], length.out = max(n.pnts, n, na.rm = TRUE))
    df <- expand.grid(x, y)
    names(df) <- c("x", "y")
    df$density <- sapply(as.data.frame(t(as.matrix(df))), function(x) pdfun(x[1], x[2], beta = beta, n = n, ...))
    return(df)
}

plot.density <- function(mdl, type = "nlm", is.R = FALSE,
                         xlim = c(-20, 150), ylim = c(0.4, 1.1), n.pnts = 201,
                         ...) {
    beta <- coef(mdl)
    if(type == "nlm") {
        n <- 50
        S <- mdl$model$Y
        sigma <- sqrt(sum(mdl$residuals^2) / mdl$df.residual)
    } else {
        n <- as.integer(mean(N <- apply(mdl$model$Y, 1, sum), na.rm = TRUE))
        S <- mdl$model$Y[ , 1] / N
        sigma <- NULL
    }
    pdfun <- switch(type, nlm = pdf.nlm, logi = pdf.logi)
    #pdfun <- switch(type, nlm = pdf.nlm, logi = pdf.logi, logi2 = pdf.logi2)
    gd = grid.density(pdfun = pdfun, beta = beta, n = n, sigma = sigma, xlim = xlim, ylim = ylim, n.pnts = n.pnts, ...)
    levelplot(density ~ x * y, data = gd, beta = beta,
              S = S,
              age = mdl$model$Age,
              panel = function(x, beta = beta, S = S, age = age, ...) {
                  panel.levelplot(x = x, ...)
                  xx <- x[1:n.pnts]
                  yy <- switch(type,
                               nlm = predict(mdl, data.frame(Age = xx)),
                               logi = logistic.f(eta.f(xx, beta = beta)))
                               #logi2 = logistic.f2(eta.f(xx, beta = beta)))
                  panel.rect(xleft = 0, ybottom = ifelse(is.R, 0, 0.5), xright = 1e5, ytop = ifelse(is.R, 100, 1), lty = "dotted")
                  #panel.rect(xleft = 0, ybottom = 0.5, xright = 1e5, ytop = 1, lty = "dotted")
                  panel.xyplot(x = xx, y = yy, type = "l", col = "black")
                  panel.xyplot(x = age, y = S, pch = 21, col = "darkgreen", fill = "green", alpha = 0.5, cex = 0.5)
              },
              colorkey = FALSE,
              xlab = "", ylab = "",
              par.settings = list(regions = list(col = rev(trellis.par.get("regions")$col))))
}


plot.mdls.1gene <- function(gene = "GRB10", l.l.M = M$simple) {
    lp <- list()
    lp$logi.S <-
        update(plot.density(l.l.M$logi.S[[gene]], type = "logi", ylim = c(0.4, 1.1), is.R = FALSE),
               main = "logi.S")
    lp$logi.S.large <-
        update(plot.density(l.l.M$logi.S[[gene]], type = "logi", xlim = c(-500, 1500), ylim = c(-0.1, 1.1), is.R = FALSE),
               main = "logi.S")
    lp$wnlm.S <-
        update(plot.density(l.l.M$wnlm.S[[gene]], type = "nlm", ylim = c(0.4, 1.1), is.R = FALSE),
               main = "wnlm.S")
    lp$wnlm.R <-
        update(plot.density(l.l.M$wnlm.R[[gene]], type = "nlm", ylim = c(-9, 110), is.R = TRUE),
               main = "wnlm.R")
    lp$unlm.S <-
        update(plot.density(l.l.M$unlm.S[[gene]], type = "nlm", ylim = c(0.4, 1.1), is.R = FALSE),
               main = "unlm.S")
    lp$unlm.R <-
        update(plot.density(l.l.M$unlm.R[[gene]], type = "nlm", ylim = c(-9, 110), is.R = TRUE),
               main = "unlm.R")
    return(lp)
}
