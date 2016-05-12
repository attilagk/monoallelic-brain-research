# data files
dir <- "../../data/ifat/age-dependence/"
files <- list(S="pop_skew_3June15.txt",
              N="pop_cov_3June15.txt",
              X="DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv",
              X2="samples.csv")
files <- lapply(files, function(f) paste0(dir, f))

             # 8 genes analyzed by Ifat
genes <- c("PEG3", "INPP5F", "SNRPN", "PWAR6", "ZDBF2", "MEG3", "ZNF331", "GRB10",
             # 5 more genes analyzed by AGK 3/2/16
             "PEG10", "SNHG14", "NAP1L5", "KCNQ1OT1", "MEST",
             # 3 more genes present in data files
             "IGF2", "NLRP2", "UBE3A")

# explanatory variables
expl.var <- c("Age.of.Death",
               "Institution",
               "Gender",
               "PMI..in.hours.",
               "Dx",
               "DLPFC_RNA_isolation..RIN", "DLPFC_RNA_isolation..RIN.2",
               "DLPFC_RNA_report..Clustered.Library.Batch",
               "Ancestry.EV.1", "Ancestry.EV.2", "Ancestry.EV.3", "Ancestry.EV.4", "Ancestry.EV.5" )

# response variables
genes.or.gsets <- c(genes, "avg8", "avg13", "avg16", "pool8", "pool13", "pool16")
names(genes.or.gsets) <- genes.or.gsets
# response variables
genes.or.gsets <- c(genes, "avg8", "avg13", "avg16", "pool8", "pool13", "pool16")
names(genes.or.gsets) <- genes.or.gsets


# the data ready for fit
d <- data4fit(thrs=50, gns=genes)
# using no filter
d.nof <- data4fit(thrs=0)
# pool data across genes of sets 8, 13, 16
d.p <- lapply(list(pool8=8, pool13=13, pool16=16), function(n) pool.gset(genes[1:n], d))
d.p.nof <- lapply(list(pool8=8, pool13=13, pool16=16), function(n) pool.gset(genes[1:n], d.nof))


# fit models
m <- list()
m.nof <- list()
# omit pool8, pool13, pool16 from this loop
for(y in genes.or.gsets[ setdiff(genes.or.gsets, names(d.p)) ] ) {
    m[[y]] <- lapply(f, function(fun) fun(y, d))
    m.nof[[y]] <- lapply(f, function(fun) fun(y, d.nof))
}
# include only pool8, pool13, pool16 in this loop
for(y in names(d.p)) {
    m[[y]] <- lapply(f, function(fun) fun(y, d.p[[y]]))
    m.nof[[y]] <- lapply(f, function(fun) fun(y, d.p.nof[[y]]))
}
rm(list = c("y"))
