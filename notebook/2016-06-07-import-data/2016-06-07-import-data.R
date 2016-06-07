# table of predictors: RNA samples (rows) and predictors (columns)
predictors <- read.csv("~/projects/monoallelic-brain/data/predictors/predictors.csv",
                       row.names = "DLPFC_RNA_ID")
# mapping between RNAseqID (present in read count tables) and DLPFC_RNA_ID (present in the table of predictors)
rna_id <- read.csv("~/projects/monoallelic-brain/data/predictors/RNAseq_ID.DLPFC_RNA_ID.csv",
                   row.names = "RNAseq_ID")
# extract those RNA samples from the table of predictors whose RNAseq_ID is known
predictors <- predictors[ as.character(rna_id$DLPFC_RNA_ID), ]
row.names(predictors) <- row.names(rna_id)
rm(rna_id) # no more necessary


# explanatory variables
expl.var <- c("Age",
               "Institution",
               "Gender",
               "PMI",
               "Dx",
               "RIN", "RIN2",
               "RNA_lib_batch",
               "Ancestry_EV.1", "Ancestry_EV.2", "Ancestry_EV.3", "Ancestry_EV.4", "Ancestry_EV.5" )

