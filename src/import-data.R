# Import observed values of predictors to be used as a basis of design matrix
# for regression.
#
# Parameters
# f.predictors: file containing the table of predictors: RNA samples (rows) and predictors (columns)
# f.rna.ids: file with mapping between RNAseqID (present in read count tables)
# and DLPFC_RNA_ID (present in the table of predictors)
#
# Value
# a data frame whose rows are observations (RNA samples/individuals) and columns predictors
#
# Details
# The input files must be obtained using the shell script ~/projects/monoallelic-brain/src/predictors-for-R
get.predictors <- function(f.predictors = "~/projects/monoallelic-brain/data/predictors/predictors.csv",
                      f.rna.ids = "~/projects/monoallelic-brain/data/predictors/RNAseq_ID.DLPFC_RNA_ID.csv") {
    predictors <- read.csv(f.predictors, row.names = "DLPFC_RNA_ID")
    rna_id <- read.csv(f.rna.ids, row.names = "RNAseq_ID")
    # extract those RNA samples from the table of predictors whose RNAseq_ID is known
    predictors <- predictors[ as.character(rna_id$DLPFC_RNA_ID), ]
    row.names(predictors) <- row.names(rna_id)
    rm(rna_id) # no more necessary
    predictors
}

# Get observed readcounts for a set of selected genes
#
# Parameters
# gene.ids: a vector of IDs (symbols) of selected genes
# data.dir: the main directory containing gene-wise subdirectories and read count files
# count.thrs: lower filtering threshold for the total readcount (summing over alleles and SNPs)
# sel.obs: the RNA IDs of selected observations (normally those on predictors)
# aggregators: a list of functions to be used to aggregate read counts over sets of genes according to some method (i.e. averaging or pooling)
#
# Value
# a list of data frames, one data frame for each gene or "aggregate", whose
# rows are observations and columns are higher (H) or lower (L) readcounts and
# derived statistics such as S and R
#
# Details
# The input directory/file structure under data.dir must be generated with the
# shell scripts summary-html2csv and all-genes-readcounts-for-R, both in ~/projects/monoallelic-brain/src/
get.readcounts <- function(gene.ids,
                           data.dir = "~/projects/monoallelic-brain/data/readcount/genes",
                           count.thrs = 50,
                           sel.obs = row.names(get.predictors()),
                           aggregators = list(identity)) {
    import.rc <- function(gene) {
        fpath <- paste0(data.dir, "/", gene, "/", gene, ".csv")
        y <- read.csv(fpath, row.names = NULL, header = FALSE, skip = 1)
        #y <- read.csv(fpath, row.names = NULL)
        #y[ , c("L", "H") ] # retain only lower and higher read counts
        y[ , c(1, 7, 8) ] # retain only lower and higher read counts
    }
    get.sel.obs <- function(counts) {
        y <- counts[ sel.obs, ]
        row.names(y) <- sel.obs
        y
    }
    derive.stats <- function(y) {
        # transform 'S' statistics into LOI_R by ranking
        S.rank <- function(S) {
            y <- rank(S, ties.method = "min", na.last = "keep")
            as.numeric(round(y / max(y, na.rm = TRUE) * 100, 0))
        }
        # total read counts
        y$N <- y$L + y$H
        # S statistic
        y$S <- y$H / y$N
        y$R <- S.rank(y$S)
        y
    }
    #Y <- lapply(gene.ids, function(g) derive.stats(get.sel.obs(import.rc(g))))
    Y <- lapply(gene.ids, import.rc)
    names(Y) <- gene.ids
    filter.rc <- function(counts) {}
    Y
}


# explanatory variables
expl.var <- c("Age",
               "Institution",
               "Gender",
               "PMI",
               "Dx",
               "RIN", "RIN2",
               "RNA_lib_batch",
               "Ancestry_EV.1", "Ancestry_EV.2", "Ancestry_EV.3", "Ancestry_EV.4", "Ancestry_EV.5" )
