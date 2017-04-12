print.av <- function(av) {
    print.av.col <- function(cname) {
        if(cname %in% c("Delta.AIC", "Delta.BIC", "Chisq"))
            fm <- formatC(av[[cname]], format = "f", width = 8, digits = 1)
        if(cname %in% c("df"))
            fm <- formatC(av[[cname]], format = "d", width = 5)
        if(cname %in% c("p.Chi"))
            fm <- formatC(av[[cname]], format = "e", width = 10, digits = 1, drop0trailing = TRUE)
        return(fm)
    }
    fav <- data.frame(sapply(names(av), print.av.col, simplify = FALSE))
    row.names(fav) <- row.names(av)
    fav
}
