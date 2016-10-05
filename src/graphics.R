my.segplot <- function(data = Betas$logi.S,
                       main = expression(paste("99 % CI for ", beta, " under logi.S")),
                       ...) {
    segplot(Gene ~ Lower.CL + Upper.CL | Coefficient,
            data = data,
            scales = list(cex = 0.5, x = list(relation = "free")),
            layout = c(6,4),
            par.settings = list(add.text = list(cex = 0.7)),
            between = list(x = 0.2),
            panel = function(x, y, ...) {
                panel.grid(h = -1, v = 0)
                panel.segplot(x, y, ...)
                panel.abline(v = 0, col = "red")
            },
            main = main, as.table = TRUE, center = Estimate, draw.bands = FALSE,
            #est = data$Estimate, prepanel = my.prepanel.segplot,
            ...)
}
# for manual setting of xlim
xl <- c(-1, 1)
my.xlim <- list(Age = 0.05 * xl,
                InstitutionPenn = 3 * xl,
                InstitutionPitt = 3 * xl,
                GenderMale = 1.5 * xl,
                PMI = 0.1 * xl,
                DxControl = 2 * xl,
                DxSCZ = 2 * xl,
                RIN = 1 * xl,
                #RIN2 = 1 * xl,
                RNA_batchA = 3 * xl,
                RNA_batchB = 3 * xl,
                RNA_batchC = 3 * xl,
                RNA_batchD = 3 * xl,
                RNA_batchE = 3 * xl,
                RNA_batchF = 3 * xl,
                RNA_batchG = 3 * xl,
                RNA_batchH = 3 * xl,
                Ancestry.1 = 20 * xl,
                Ancestry.2 = 20 * xl,
                Ancestry.3 = 20 * xl,
                Ancestry.4 = 20 * xl,
                Ancestry.5 = 20 * xl)


my.segplot2 <- function(coef, type, main = type, skip = FALSE, layout = c(6, 5), ...) {
    segplot(ordered(Permutation, levels = rev(levels(Permutation))) ~ Lower.CL + Upper.CL | Gene, data = Betas,
            subset = Coefficient %in% coef & Model %in% type,
            level = factor(! Permutation == "U"),
            colorkey = FALSE,
            par.settings = list(add.text = list(cex = 0.8)),
            panel = function(x, y, ...) {
                panel.segplot(x, y, ...)
                panel.abline(v = 0, col = "red")
            },
            xlab = eval(substitute(expression(paste(beta, "[ ", coef, " ]")), list(coef = coef))),
            ylab = "permutation", main = main,
            skip = skip,
            layout = layout,
            scales = list(draw = FALSE, x = list(relation = "free")))[c(1:30)[! skip]]
}

beta0densityplot <- function(coef = "Age", mtype = "wnlm.Q", data = Betas, ...) {
    densityplot(~ Estimate | Gene, data = data, subset = Coefficient == coef & Model == mtype,
                par.settings = list(add.text = list(cex = 0.8)),
                scales = list(relation = "free", draw = FALSE),
                U = data$Permutation == "U",
                panel = function(x, ..., U = U, subscripts) {
                    panel.densityplot(x, ..., plot.points = FALSE, ref = TRUE)
                    panel.abline(v = 0, col = "red")
                    beta.U <- x[U[subscripts]]
                    p.val <- beta.U
                    panel.abline(v = beta.U, col = trellis.par.get("plot.line")$col, lty = "dotted")
                    panel.text(x = beta.U, y = 0,
                               labels = format(get.p.val(beta.U, x[! U[subscripts]], do.norm = TRUE),
                                               digits = 2, scientific = FALSE))
                }, main = coef,
                ...)
}
