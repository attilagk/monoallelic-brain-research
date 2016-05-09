# list of simple regression models
m.smpl <- lapply(genes.or.gsets[1:19], fit.S.x1, d)
for(g in genes.or.gsets[20:22])
    m.smpl[[g]] <- fit.S.x1(g, dt = d.p[[g]])
# get estimates and p-values for the age regression coefficient
betas <- coefs.from.md.lists(m.smpl, m, prop = 1)
pvals <- coefs.from.md.lists(m.smpl, m, prop = 4)

# test the effect of filtering
rel.aic <- - t(sapply(c("nlm.R","nlm.S", "logi.S", "logi2.S"), function(x) m2m(m, m.nof, x, x)))
n.i <- t(cbind(retained <-
    sapply(genes.or.gsets, function(g) sum(! is.na(m[[g]]$nlm.R$data[[
                                                   paste0("N_", g) ]]))),
             filtered <-
                 sapply(genes.or.gsets,
                        function(g) sum(! is.na(m.nof[[g]]$nlm.R$data[[ paste0("N_", g) ]])))
             ))

