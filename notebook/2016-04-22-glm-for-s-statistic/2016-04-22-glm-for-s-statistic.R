################################# functions #################################

old.fun <- list(
# read data and set row names
get.data = function(filename, id.string='ID', ...) {
    df <- read.delim(filename, na.strings=c('NA', 'ND'), ...)
    row.names(df) <- df[, id.string]
    return(df)
},

# Copies and pastes the the columns 'select.col' from one data frame to
# another, then sets row names
paste.column = function(from.df, to.df, select.col=c("RNAseq_ID"),
                         row.n="RNAseq_ID") {
    from.df <- from.df[select.col]
    to.df <- merge(from.df, to.df, by="row.names")
    row.names(to.df) <- to.df[[row.n]]
    # remove annoying "Row.names" column inserted by the merge function
    to.df <- to.df[, ! grepl("Row.names", names(to.df))]
    return(to.df)
},

# clean.obs variables by removing contaminating other variables
clean.obs = function(df, remove.cols=c("age")) {
    # remove non-numeric columns
    df <- df[, sapply(df, is.numeric)]
    # remove columns specified by remove.cols
    df <- df[, setdiff(names(df), remove.cols)]
    # sort according to row names
    df <- df[sort.list(row.names(df)), ]
    return(df)
},

# replace elements in df.A where df.B <= thrs
filter.obs = function(df.A, df.B, thrs=50) {
    if(! all( dim(df.A) == dim(df.B))) return(NaN)
    A <- as.vector(as.matrix(df.A))
    B <- as.vector(as.matrix(df.B))
    df <- as.data.frame(matrix(ifelse(B > thrs, A, NA), nrow=nrow(df.A),
                         ncol=ncol(df.A)))
    names(df) <- names(df.A)
    row.names(df) <- row.names(df.A)
    return(df)
},

# transform 'S' statistics into LOI_R based on selected 'genes'
loir.tform = function(S, genes, do.rank=TRUE) {
    loir.rank <- function(x) {
        y <- rank(x, ties.method = "min", na.last = "keep")
        as.numeric(round(y / max(y, na.rm = TRUE) * 100, 0))
    }
    if(do.rank)
        foo <- loir.rank
    else
        foo <- identity
    S.mat <- as.matrix(as.data.frame(lapply(S[, genes, drop=FALSE], foo)))
    S.avg <- as.data.frame(apply(S.mat, 1, mean, na.rm=TRUE))
    row.names(S.avg) <- row.names(S)
    names(S.avg) <- "response"
    return(S.avg)
},

# high level function preparing data for fitting (filtering, column pasting)
prepare.fit.data = function(obs=S.stat, input=expl.var,
                             filter.df.B=tot.read.n, filter.thrs=50,
                             genes=genes[1:8], tform=TRUE,
                             weights=as.data.frame(rep(1, nrow(tot.read.n))))
{
    # filter
    obs <- filter.obs(obs, filter.df.B, thrs=filter.thrs)
    # merge weights into the response (i.e. the transformed obs)
    names(weights) <- "weights"
    row.names(weights) <- row.names(filter.df.B)
    tform.fun <- function(S) loir.tform(S, genes, do.rank=tform)
    response <- tform.fun(obs)
    response <- merge(response, weights, by="row.names")
    response <- response[-1]
    row.names(response) <- row.names(S.stat)
    # merge with expl.var
    fit.data <- paste.column(response, expl.var,
                             select.col=c("response", "weights"), row.n="RNAseq_ID")
    return(fit.data)
},

# fit glm for a selected 'gene'
do.fit = function(gene, family=gaussian, filter.thrs=50,
                   tform=TRUE, weights=tot.read.n[gene]) {
    dat <- prepare.fit.data(obs=S.stat, input=expl.var, filter.df.B=tot.read.n,
                            filter.thrs=filter.thrs, genes=gene, tform=tform,
                            weights=weights)
    glm(formula=lin.pred, family=family, data=dat, weights=weights)
}
)


################################# newer functions #################################

# make a formula for fit given a string 'response' character vector 'expl' of explanatory variables
mk.form <- function(response, expl = expl.var)
    as.formula(paste(response, "~", paste(expl, collapse = " + ")))

# reads data from files, cleans, filters and transforms them
data4fit <- function(S.f=files$S, N.f=files$N, X.f=files$X,
                     X2.f=files$X2, thrs=50, gns=genes) {
    # read data and set row names
    get.data <- function(filename, id.string='ID', ...) {
        df <- read.delim(filename, na.strings=c('NA', 'ND'), ...)
        row.names(df) <- df[, id.string]
        return(df)
    }
    # Copies and pastes the the columns 'select.col' from one data frame to
    # another, then sets row names according to points in column 'row.n'
    paste.column <- function(from.df, to.df, select.col=c("RNAseq_ID"),
                             row.n="RNAseq_ID") {
        from.df <- from.df[select.col]
        to.df <- merge(from.df, to.df, by="row.names")
        row.names(to.df) <- to.df[[row.n]]
        # remove annoying "Row.names" column inserted by the merge function
        to.df <- to.df[, ! grepl("Row.names", names(to.df))]
        return(to.df)
    }
    # clean.obs variables by removing contaminating other variables
    clean.obs <- function(df, remove.cols=c("age")) {
        # remove non-numeric columns
        df <- df[, sapply(df, is.numeric)]
        # remove columns specified by remove.cols
        df <- df[, setdiff(names(df), remove.cols)]
        return(df)
    }
    # replace elements in df.A where df.B <= thrs
    filter.S <- function(df.A, df.B, thrs=thrs) {
        if(! all( dim(df.A) == dim(df.B))) return(NaN)
        A <- as.vector(as.matrix(df.A))
        B <- as.vector(as.matrix(df.B))
        df <- as.data.frame(matrix(ifelse(B > thrs, A, NA), nrow=nrow(df.A),
                                   ncol=ncol(df.A)))
        names(df) <- names(df.A)
        row.names(df) <- row.names(df.A)
        return(df)
    }
    # transform 'S' statistics into LOI_R by ranking
    loir.rank <- function(S) {
        y <- rank(S, ties.method = "min", na.last = "keep")
        as.numeric(round(y / max(y, na.rm = TRUE) * 100, 0))
    }
    # averages any statistic Y based on selected 'gns'
    stat.Y <- function(Y, gns=gns, stat=mean, col.name="Y_stat", ...) {
        Y.stat <- as.data.frame(apply(as.matrix(Y[, gns]), 1, stat, ...))
        row.names(Y.stat) <- row.names(Y)
        names(Y.stat) <- col.name
        return(Y.stat)
    }
    # for scaled logistic function
    rescale <- function(mat) {
        n <- mat[, 1] + mat[, 2]
        p <- mat[, 1] / n
        x <- as.integer((p * 2 - 1) * n)
        c2 <- cbind(x[], as.integer(n - x[]))
        c2[ c2 < 0 & ! is.na(c2)] <- 0
        return(c2)
    }
    # get data & clean them
    N <- clean.obs(get.data(N.f, sep="\t"))
    S <- clean.obs(get.data(S.f, sep="\t"))
    X <- get.data(X.f, id.string="DLPFC_RNA_isolation..Sample.RNA.ID", sep="\t")
    X2 <- get.data(X2.f, id.string="DLPFC_RNA_ID", sep=",")
    # clean up explanatory variables
    X <- paste.column(X2, X, select.col=c("RNAseq_ID"), row.n="RNAseq_ID")
    # sort according to row names (i.e. "RNAseq_ID")
    S <- S[sort.list(row.names(S)), ]
    N <- N[sort.list(row.names(N)), ]
    X <- X[sort.list(row.names(X)), ]
    # filter S
    S <- filter.S(S, N, thrs=thrs)
    N <- filter.S(N, N, thrs=thrs)
    # ranks and their averages
    R <- as.data.frame(lapply(S, loir.rank))
    R.16 <- stat.Y(R, gns=gns, stat=mean, col.name="R_g16", na.rm=TRUE)
    R.13 <- stat.Y(R, gns=gns[1:13], stat=mean, col.name="R_g13", na.rm=TRUE)
    R.8 <- stat.Y(R, gns=gns[1:8], stat=mean, col.name="R_g8", na.rm=TRUE)
    # higher read counts and their sums
    H <- as.data.frame(lapply(names(S), function(g) as.integer(S[[g]] * N[[g]])))
    names(H) <- names(S)
    H.16 <- stat.Y(H, gns=gns, stat=sum, col.name="H_g16", na.rm=TRUE)
    H.13 <- stat.Y(H, gns=gns[1:13], stat=sum, col.name="H_g13", na.rm=TRUE)
    H.8 <- stat.Y(H, gns=gns[1:8], stat=sum, col.name="H_g8", na.rm=TRUE)
    # lower read counts and their sums
    L <- as.data.frame(lapply(names(S), function(g) N[[g]] - H[[g]]))
    names(L) <- names(S)
    L.16 <- stat.Y(L, gns=gns, stat=sum, col.name="L_g16", na.rm=TRUE)
    L.13 <- stat.Y(L, gns=gns[1:13], stat=sum, col.name="L_g13", na.rm=TRUE)
    L.8 <- stat.Y(L, gns=gns[1:8], stat=sum, col.name="L_g8", na.rm=TRUE)
    # total read count sums
    N.16 <- stat.Y(N, gns=gns, stat=sum, col.name="N_g16", na.rm=TRUE)
    N.13 <- stat.Y(N, gns=gns[1:13], stat=sum, col.name="N_g13", na.rm=TRUE)
    N.8 <- stat.Y(N, gns=gns[1:8], stat=sum, col.name="N_g8", na.rm=TRUE)
    # S averages
    S.16 <- H.16 / N.16; names(S.16) <- "S_g16"
    S.13 <- H.13 / N.13; names(S.13) <- "S_g13"
    S.8 <- H.8 / N.8; names(S.8) <- "S_g8"
    # give columns unique names
    names(S) <- paste("S", sep="_", names(S))
    names(N) <- paste("N", sep="_", names(N))
    names(R) <- paste("R", sep="_", names(R))
    names(H) <- paste("H", sep="_", names(H))
    names(L) <- paste("L", sep="_", names(L))
    df <- cbind(S, S.16, S.13, S.8,
                 R, R.16, R.13, R.8,
                 H, H.16, H.13, H.8,
                 L, L.16, L.13, L.8,
                 N, N.16, N.13, N.8,
                 X)
    # two-column matrices for fitting with binomial response distribution
    for(r in c(gns, "g8", "g13", "g16")) {
        df[[ paste0("C_", r) ]] <- cbind( df[[ paste0("H_", r) ]], df[[ paste0("L_", r) ]])
        df[[ paste0("C2_", r) ]] <- rescale(df[[ paste0("C_", r) ]])
    }
    return(df)
}

f <- list(
          nlm = function(y, d) { # normal linear model with rank R as response
              glm(formula = mk.form(paste0("R_", y)), family = gaussian, data = d)
          } ,
          nlm2 = function(y, d) { # as above but with the R's lm function instead of glm
              lm(formula = mk.form(paste0("R_", y)), data = d)
          } ,
          logi = function(y, d) { # logistic regression with the S statistic as response
              glm(formula = mk.form(paste0("C_", y)), family = binomial, data = d)
          } ,
          logi2 = function(y, d) { # as above but with rescaled and offset logistic link function
              glm(formula = mk.form(paste0("C2_", y)), family = binomial, data = d)
          } )

################################# analysis #################################

# data files
files <- list(S="pop_skew_3June15.txt",
              N="pop_cov_3June15.txt",
              X="DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv",
              X2="samples.csv")

# selected genes
             # 8 genes analyzed by Ifat
genes <- c("PEG3", "INPP5F", "SNRPN", "PWAR6", "ZDBF2", "MEG3", "ZNF331", "GRB10",
             # 5 more genes analyzed by AGK 3/2/16
             "PEG10", "SNHG14", "NAP1L5", "KCNQ1OT1", "MEST",
             # 3 more genes present in data files
             "IGF2", "NLRP2", "UBE3A")
names(genes) <- genes

# explanatory variables
expl.var <- c("`Age.of.Death`",
               "Institution",
               "Gender",
               "`PMI..in.hours.`",
               "Dx",
               "`DLPFC_RNA_isolation..RIN`", "`DLPFC_RNA_isolation..RIN.2`",
               "`DLPFC_RNA_report..Clustered.Library.Batch`",
               "`Ancestry.EV.1`", "`Ancestry.EV.2`", "`Ancestry.EV.3`", "`Ancestry.EV.4`", "`Ancestry.EV.5`" )

# response variables
responses <- c(genes, "g8", "g13", "g16")
names(responses) <- responses


# the data ready for fit
d <- data4fit(thrs=50, gns=genes)
# fitted models
m <- list()
for(y in responses) {
    m[[y]]$nlm <- f$nlm(y, d)
    m[[y]]$nlm2 <- f$nlm2(y, d)
    m[[y]]$logi <- f$logi(y, d)
    m[[y]]$logi2 <- f$logi2(y, d)
}
#m$nlm.ifat[["8"]] <- fit.lm(transform.data(genes[1:13], genes[1:8]), do.thrs=FALSE)

# using no filter
d.nof <- data4fit(thrs=0)
m.nof <- list()
for(y in responses) {
    m.nof[[y]]$nlm <- f$nlm(y, d=d.nof)
    m.nof[[y]]$nlm2 <- f$nlm2(y, d=d.nof)
    m.nof[[y]]$logi <- f$logi(y, d=d.nof)
    m.nof[[y]]$logi2 <- f$logi2(y, d=d.nof)
}
rm(list = c("y"))

# compare model 1 to model 2; ml1 and ml2 are lists of models fitted to
# different data; mf1 and mf2 are strings of model families ("nlm", "logi", "logi2")
# fun compares the models's AIC values, default is scaled difference
m2m <- function(ml1, ml2, mf1, mf2, fun = function(x, y) (x - y) / abs(x)) {
    sapply(responses, function(r) fun(ml1[[r]][[mf1]]$aic, ml2[[r]][[mf2]]$aic))
}
# test the effect of filtering
signif(sapply(c("nlm", "logi", "logi2"), function(x) m2m(m, m.nof, x, x)), 3)
# compare logi to nlm
signif(m2m(m, m, "logi", "nlm"), 3)
# compare logi2 to logi
signif(m2m(m, m, "logi2", "logi"), 3)
signif(m2m(m.nof, m.nof, "logi2", "logi"), 3)




################################# TODO #################################

## aggregating data from multiple genes by pooling together instead of averaging
if(FALSE) {
d.ag <- list()
m.ag <- list()
for(n in c(8, 13, 16)) {
    gn <- paste0("g", n)
    d.ag[[gn]] <- lapply(d, rep, n)
    d.ag[[gn]][[paste0("R_g", n)]] <- as.vector(unlist(d[ paste0( "R_", genes[1:n]) ]))
    d.ag[[gn]][[paste0("H_g", n)]] <- as.vector(unlist(d[ paste0( "H_", genes[1:n]) ]))
    d.ag[[gn]][[paste0("N_g", n)]] <- as.vector(unlist(d[ paste0( "N_", genes[1:n]) ]))
    d.ag[[gn]][[paste0("S_g", n)]] <- d.ag[[gn]][[paste0("H_g", n)]] / d.ag[[gn]][[paste0("N_g", n)]]
    m.ag[[gn]]$nlm <- f$nlm(gn, d.ag[[gn]])
    m.ag[[gn]]$nlm2 <- f$nlm2(gn, d.ag[[gn]])
    # IMPORTANT: CALLING f$logi THROWS AN ERROR (DEBUGGING FAILED):
    # E.G. THE FOLLOWING CALL
    # m.ag[[gn]]$logi3 <- f$logi(as.character(n), d.ag[[gn]])
    # PRODUCES THE FOLLOWING ERROR
    # Error in model.frame.default(formula = mk.form(paste0("S_g", y)), data = d,
    # : variable lengths differ (found for '(weights)')
    # USING DIRECTLY glm INSTEAD
    m.ag[[gn]]$logi <-
        glm(formula = mk.form(paste0("C_", gn)), family = binomial, data = d.ag[[gn]])
    m.ag[[gn]]$logi2 <-
        glm(formula = mk.form(paste0("C2_", gn)), family = binomial, data = d.ag[[gn]])
}
rm(list = c("n", "gn"))
}
