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
    R.16 <- stat.Y(R, gns=gns, stat=mean, col.name="R_avg16", na.rm=TRUE)
    R.13 <- stat.Y(R, gns=gns[1:13], stat=mean, col.name="R_avg13", na.rm=TRUE)
    R.8 <- stat.Y(R, gns=gns[1:8], stat=mean, col.name="R_avg8", na.rm=TRUE)
    # higher read counts and their sums
    H <- as.data.frame(lapply(names(S), function(g) as.integer(S[[g]] * N[[g]])))
    names(H) <- names(S)
    H.16 <- stat.Y(H, gns=gns, stat=sum, col.name="H_avg16", na.rm=TRUE)
    H.13 <- stat.Y(H, gns=gns[1:13], stat=sum, col.name="H_avg13", na.rm=TRUE)
    H.8 <- stat.Y(H, gns=gns[1:8], stat=sum, col.name="H_avg8", na.rm=TRUE)
    # lower read counts and their sums
    L <- as.data.frame(lapply(names(S), function(g) N[[g]] - H[[g]]))
    names(L) <- names(S)
    L.16 <- stat.Y(L, gns=gns, stat=sum, col.name="L_avg16", na.rm=TRUE)
    L.13 <- stat.Y(L, gns=gns[1:13], stat=sum, col.name="L_avg13", na.rm=TRUE)
    L.8 <- stat.Y(L, gns=gns[1:8], stat=sum, col.name="L_avg8", na.rm=TRUE)
    # total read count sums
    N.16 <- stat.Y(N, gns=gns, stat=sum, col.name="N_avg16", na.rm=TRUE)
    N.13 <- stat.Y(N, gns=gns[1:13], stat=sum, col.name="N_avg13", na.rm=TRUE)
    N.8 <- stat.Y(N, gns=gns[1:8], stat=sum, col.name="N_avg8", na.rm=TRUE)
    # S averages
    S.16 <- H.16 / N.16; names(S.16) <- "S_avg16"
    S.13 <- H.13 / N.13; names(S.13) <- "S_avg13"
    S.8 <- H.8 / N.8; names(S.8) <- "S_avg8"
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
    for(r in c(gns, "avg8", "avg13", "avg16")) {
        df[[ paste0("C_", r) ]] <- cbind( df[[ paste0("H_", r) ]], df[[ paste0("L_", r) ]])
        df[[ paste0("C2_", r) ]] <- scale4log2(df[[ paste0("C_", r) ]])
    }
    return(df)
}

# scale counts in two-column matrix for logistic regression under the logi2.S model
scale4log2 <- function(mat) {
    n <- mat[, 1] + mat[, 2]
    p <- mat[, 1] / n
    x <- as.integer((p * 2 - 1) * n)
    c2 <- cbind(x[], n - x[])
    c2[ c2 < 0 & ! is.na(c2)] <- 0
    return(c2)
}

# data aggregation by pooling points from data frame 'din' across all genes in gene set 'gns'
pool.gset <- function(gns, din) {
    n <- length(gns)
    resps <- "SNRHL" # string with initials of response variables
    # extract explanatory variables in din by omitting response variables
    expl <- grep(paste0("^[", resps, "C]2?_"), names(din), value=TRUE, invert=TRUE)
    # replicate points for each explanatory variable
    dout <- as.data.frame(lapply(din[ expl ], rep, n))
    # for each response pool points across genes
    for(s in unlist(strsplit(resps, ""))) {
        dout[[paste0(s, "_pool", n)]] <- as.vector(unlist(din[ paste0(s, "_", gns) ]))
    }
    # the higher and lower read counts go into 2 column matrices with optional # rescaling
    dout[[paste0("C_pool", n)]] <- C <- cbind(dout[[paste0("H_pool", n)]], dout[[paste0("L_pool", n)]]) 
    dout[[paste0("C2_pool", n)]] <- scale4log2(C) 
    dout
}

# Wrapper for various regression models; argument 'y' is the name of the gene
# or gene set, e.g 'PEG3' or 'avg13' and 'd' is the data frame of the fitted # data
f <- list(
          nlm.R = function(y, d) { # normal linear model with rank R as response
              glm(formula = mk.form(paste0("R_", y)), family = gaussian, data = d)
          } ,
          nlm2.R = function(y, d) { # as above but with the R's lm function instead of glm
              lm(formula = mk.form(paste0("R_", y)), data = d)
          } ,
          nlm.S = function(y, d) { # normal linear model with S statistic as response
              glm(formula = mk.form(paste0("S_", y)), family = gaussian, data = d)
          } ,
          nlm2.S = function(y, d) { # as above but with the R's lm function instead of glm
              lm(formula = mk.form(paste0("S_", y)), data = d)
          } ,
          logi.S = function(y, d) { # logistic regression with the S statistic as response
              glm(formula = mk.form(paste0("C_", y)), family = binomial, data = d)
          } ,
          logi2.S = function(y, d) { # as above but with rescaled and offset logistic link function
              glm(formula = mk.form(paste0("C2_", y)), family = binomial, data = d)
          } )

################################# analysis #################################
