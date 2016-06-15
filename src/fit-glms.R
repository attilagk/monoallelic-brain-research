# predictors (a.k.a. explanatory variables)
predictors <- c("Age",
               "Institution",
               "Gender",
               "PMI",
               "Dx",
               "RIN", "RIN2",
               "RNA_lib_batch",
               "Ancestry_EV.1", "Ancestry_EV.2", "Ancestry_EV.3", "Ancestry_EV.4", "Ancestry_EV.5" )

do.fit <- function(y = Y[[1]]$S,
                   X = E,
                   preds = predictors,
                   ...) {
    X$Y <- y
    fm <- as.formula(paste("Y", "~", paste0(preds, collapse = " + ")))
    glm(fm, data = X, ...)
}
