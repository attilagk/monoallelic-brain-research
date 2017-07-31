get.coef <- function(batch = "Gene", model = M5, fun = ranef, reverse = TRUE, prettify = TRUE) {
    df <- fun(model)[[batch]]
    lv <-
        if (reverse) rev(rownames(df))
        else rownames(df)
    f <- factor(rownames(df), levels = lv, ordered = TRUE)
    df <- stack(df)
    df[["batch.levels"]] <- f
    df[["batch"]] <- batch
    if (prettify)
        df$ind <- factor(pretty.batch(df$ind, batch), levels = pretty.batch(levels(df$ind), batch))
    return(df)
}

pretty.batch <- function(batch.levels, batch.name = "Gene", interceptAs1 = TRUE) {
    x <- gsub("^.*\\(([^)]+)\\).*$", paste0("\\(\\1\\|", batch.name, "\\)"), batch.levels)
    if (interceptAs1) return(gsub("Intercept", "1", x))
    else return(x)
}

contrast.coef <- function(cf = get.coef("Gender:Gene"), e.var = "Gender", df = dat) {
    lev1 <- levels(df[[e.var]])[1]
    helper <- function(y) {
        x <- cf[grep(y, cf$batch.levels), ]
        x$batch.levels <- gsub(y, "", x$batch.levels)
        x$batch.levels <- gsub(":", "", x$batch.levels)
        x$batch <- gsub(e.var, "", x$batch)
        x$batch <- gsub(":", "", x$batch)
        return(x)
    }
    lapply(levels(df[[e.var]]), helper)
}

mydotplot <- function(x = get.coef(), ...) {
    n <- names(x)
    fo <- as.formula(paste0("`", n[3], "`", " ~ `", n[1], "` | `", n[2], "`"))
    dotplot(fo, data = x, scales = list(x = list(relation = "free")), ...)
}
