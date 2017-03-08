# Import observed values of predictors to be used as a basis of design matrix
# for regression.
#
# Arguments
# f.predictors: file containing the table of predictors: RNA samples (rows) and predictors (columns)
# f.rna.ids: file with mapping between RNAseqID (present in read count tables) and DLPFC_RNA_ID (present in the table of predictors)
# adj.levels: logical; whether to adjust levels of certain factors e.g. Dx
#
# Value
# a data frame whose rows are observations (RNA samples/individuals) and columns predictors
#
# Details
# The input files must be obtained using the shell script ~/projects/monoallelic-brain/src/predictors-for-R
get.predictors <- function(f.predictors = "~/projects/monoallelic-brain/data/predictors/predictors.csv",
                      f.rna.ids = "~/projects/monoallelic-brain/data/predictors/RNAseq_ID.DLPFC_RNA_ID.csv",
                      adj.levels = TRUE) {
    predictors <- read.csv(f.predictors, row.names = "DLPFC_RNA_ID")
    rna_id <- read.csv(f.rna.ids, row.names = "RNAseq_ID")
    # extract those RNA samples from the table of predictors whose RNAseq_ID is known
    predictors <- predictors[ as.character(rna_id$DLPFC_RNA_ID), ]
    row.names(predictors) <- row.names(rna_id)
    rm(rna_id) # no more necessary
    if(! adj.levels) return(predictors)
    adj.levels.predictors(predictors)
}

# Adjust levels of certain factors (e.g. set Control as first (baseline) level for Dx)
#
# Arguments
# X: the matrix of observed predictor values
#
# Value
# The adjusted X
adj.levels.predictors <- function(X) {
    # an alternative implementation would set contrasts instead of reordering factor levels
    X$Dx <- factor(X$Dx, levels = c("Control", "SCZ", "AFF"))
    X$RNA_batch <- factor(X$RNA_batch, levels = c(LETTERS[1:8], "0"))
    return(X)
}


# Get observed readcounts for a set of selected genes
#
# Arguments
# gene.ids: a vector of IDs (symbols) of selected genes
# data.dir: the main directory containing gene-wise subdirectories and read count files
# count.thrs: lower filtering threshold for the total readcount (summing over alleles and SNPs)
# sel.obs: the RNA IDs of selected observations (normally those on predictors)
# g.subsets: a named list of subsets of genes for aggregation over each subset using each of two methods: weighted and unweighted average
# rm.conflict: "a filter giving priority to even a single SNP showing biallelic expression" (see Ifat's old ms)
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
                           count.thrs = 0,
                           sel.obs = row.names(get.predictors()),
                           g.subsets = list(A.8 = gene.ids[1:8], A = gene.ids),
                           rm.conflict = FALSE,
                           sel.var = c("L", "H", "N", "S", "R", "Q")) {
    import.rc <- function(gene) {
        fpath <- paste0(data.dir, "/", gene, "/", gene, ".csv")
        # shell command to select 'Sample RNA ID', 'L' and 'H' columns using GNU's cut;
        cmd <- "cut --delimiter=, --fields=1,5,7,8"
        cmds <- paste(cmd, fpath)
        read.csv(pipe(cmds), row.names = 1)
    }
    conflict.rm <- function(counts) {
        if(! rm.conflict) return(counts)
        is.confl <- counts$Allelic.state == "1_conflict"
        counts[ is.confl, "L" ] <- NA
        counts[ is.confl, "H" ] <- NA
        return(counts)
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
        shift.log.transf <- function(H, L, pseudo.cnt = 1) {
            S <- H / (H + L + pseudo.cnt)
            - log(1 - S)
        }
        # total read counts
        y$N <- y$L + y$H
        # S statistic
        y$S <- y$H / y$N
        y$R <- S.rank(y$S)
        y$Q <- shift.log.transf(y$H, y$L, pseudo.cnt = 1)
        y[sel.var]
    }
    aggregator <- function(g.subset, Z, vars = c("S", "R"), fun = mean) {
        helper <- function(var) {
            M <- as.matrix(as.data.frame(lapply(Z[g.subset], function(y) y[[var]])))
            apply(M, 1, fun, na.rm = TRUE)
        }
        names(vars) <- vars
        res <- data.frame(lapply(vars, helper))
        row.names(res) <- sel.obs
        res
    }
    filter.rc <- function(y) {
        # observations pass filter if total read count is > threshold
        passed <- y$L + y$H > count.thrs
        # apply filter to each statistic, i.e. each column of the data frame
        z <- data.frame(lapply(y, function(x) ifelse(passed, x, NA)))
        row.names(z) <- sel.obs
        z
    }
    # import, remove "1_conflict"s and adjust set of observations
    Y <- lapply(gene.ids, function(g) get.sel.obs(conflict.rm(import.rc(g))))
    names(Y) <- gene.ids
    # 1st aggregation: take 'W'eighted average by summing up L and H over each gene subset
    Y.pre <- lapply(g.subsets, aggregator, Y, vars = c("L", "H"), fun = sum)
    # filter and derive the additional stats
    Y <- lapply(Y, function(y) derive.stats(filter.rc(y)))
    Y.pre <- lapply(Y.pre, function(y) derive.stats(filter.rc(y)))
    # 2nd aggregation: take 'U'nweigthed average of S and R over each gene subset
    Y.post <- lapply(g.subsets, aggregator, Y, vars = c("S", "R", "Q"), fun = mean)
    # name columns and return gene-wise data Y together with 'WA' and 'UA' aggregates
    names(Y.pre) <- paste0(rep("W", length(Y.pre)), names(g.subsets))
    names(Y.post) <- paste0(rep("U", length(Y.pre)), names(g.subsets))
    c(Y, Y.pre, Y.post)
}

             # 8 genes analyzed by Ifat
gene.ids <- c("PEG3", "INPP5F", "SNRPN", "PWAR6", "ZDBF2", "MEG3", "ZNF331", "GRB10",
             # 5 more genes analyzed by AGK 3/2/16
             "PEG10", "SNHG14", "NAP1L5", "KCNQ1OT1", "MEST",
             # 3 more genes present in data files
             "IGF2", "NLRP2", "UBE3A",
             # 'green' novel 1 MB imprinted genes; note that PWAR6 is already
             # included above
             "TMEM261P1", "AL132709.5", "RP11-909M7.3", "SNORD116-20", "RP13-487P22.1", "hsa-mir-335", "PWRN1")


# merge read counts with data on predictors into a data frame for fitting
#
# Arguments
# gene.ids: the ids for selected genes
# Y: read count data
# E: the wide format predictor data resulted from a call to the get.predictors function
#
# Value: a data frame, which contains also the 2 column matrix element H.N for
# fitting a logistic regression model
merge.data <- function(gene.ids = unlist(read.csv("../../data/genes.regression.new", as.is = TRUE)),
                       Y = get.readcounts(gene.ids)[gene.ids],
                       E = get.predictors()[scan(file = "../../data/e.vars", what = "character")]) {
    l <- lapply(gene.ids, function(g)
                cbind(y <- Y[[g]], data.frame(Individual = rownames(y), Gene = g), E))
    df <- do.call(rbind, l)
    cbind(data.frame(H.N = I(with(df, matrix(c(H, N), ncol = 2)))), df)
}


# reshape data on (transformed) read counts and predictors to long format
#
# Arguments
# gene.ids: the ids for selected genes
# Z: the wide format read count data resulted from a call to the get.readcounts function
# E: the wide format predictor data resulted from a call to the get.predictors function
# e.v: selected explanatory variables
#
# Value: a data frame in the long format, where the same variables for
# different genes (e.g. total read count N) are stacked on each other
reshape.Y2long <- function(gene.ids, Z = get.readcounts(gene.ids)[gene.ids], E = get.predictors(), e.v = e.vars) {
    stats <- names(Z[[1]])
    names(stats) <- stats
    l <- lapply(stats, function(s) stack(lapply(Z, getElement, s)))
    Y.long <- data.frame(sapply(l, getElement, 1))
    Y.long$Gene <- factor(l[[1]][ , 2, drop = TRUE])
    Y.long <- cbind(Y.long, E[ rep(seq_len(nrow(E)), length(Z)), e.v ])
    Y.long$Gene <- factor(Y.long$Gene, levels = names(Z), ordered = TRUE)
    return(Y.long)
}
