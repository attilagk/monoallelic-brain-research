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
                RIN = 10 * xl,
                RIN2 = 1 * xl,
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

