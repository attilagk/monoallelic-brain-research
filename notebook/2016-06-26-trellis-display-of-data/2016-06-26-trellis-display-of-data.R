my.update <- function(dp, ...) {
    update(dp, scales = list(x = list(draw = FALSE, relation = "free")),
           par.settings = list(add.text = list(cex = 0.8)),
           xlab = "read count ratio, S",
           ...)
}
