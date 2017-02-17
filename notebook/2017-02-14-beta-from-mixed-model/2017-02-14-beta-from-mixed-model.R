# annotation keys
my.key <-
    list(coefficient =
             list(columns = 5, text = list(paste0("(", seq_along(levels(cf$coefficient)), ") ", levels(cf$coefficient)), cex = 0.7)),
         predictor =
             list(columns = 4, text = list(paste0("(", seq_along(e.vars), ") ", e.vars), cex = 0.7)),
         gene =
             list(columns = 5, text = list(paste0("(", seq_along(levels(cf$gene)), ") ", levels(cf$gene)), cex = 0.7)))

# high level plotting function
my.plot <- function(x = "unlm.Q", y = "fixed.1", group.by = "gene", lbl.type = "coefficient", dt = cf, ...) {
    my.panel <- function(..., lbl) {
        panel.abline(a = 0, b = 1, lty = "dotted", col = "gray")
        lbl <-
            if(lbl.type == "gene") dt$gene
            else dt[[1]]
        panel.text(..., pch = lbl, cex = 0.7)
    }
    fm <- formula(paste(y, "~", x, "|", group.by))
    xyplot(fm, data = dt, panel = my.panel,
           scales = switch(group.by, gene = "same", coefficient = "free", predictor = "free"),
           layout = switch(group.by, gene = c(5, 6), coefficient = c(4, 6), predictor = c(4, 4)),
           key = my.key[[lbl.type]], ...)
}

reshape.x.y.y.hat <- function(m, e.v = e.vars) {
    XY.obs <- reshape(m$model, varying = e.v, v.names = "X", timevar = "predictor",
                      ids = row.names(df), times = e.v, direction = "long")
    Y.hat <- reshape(df <- data.frame(predict(m, type = "terms", se.fit = FALSE)),
                     varying = e.v, v.names = "Y.hat", timevar = "predictor",
                     ids = row.names(df), times = e.v, direction = "long")
    cbind(XY.obs[ c("id", "predictor", "X", "Y") ], Y.hat[ "Y.hat" ])
}
