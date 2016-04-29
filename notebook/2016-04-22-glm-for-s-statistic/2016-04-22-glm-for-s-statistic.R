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

################################# INTO RMARKDOWN #################################

if(FALSE){

# Fits simple models S ~ x2 (ignoring all other explanatory variables). Returns a list of fitted models
fit.S.x1 <- function(g = genes.or.gsets[1], dt = d, x1 = "Age.of.Death") {
    md.l <- list()
    md.l$nlm.R <- glm(mk.form(paste0("R_", g), x1), data=dt)
    md.l$nlm2.R <- lm(mk.form(paste0("R_", g), x1), data=dt)
    md.l$nlm.S <- glm(mk.form(paste0("S_", g), x1), data=dt)
    md.l$logi.S <- glm(mk.form(paste0("C_", g), x1), data=dt, family=binomial)
    md.l$logi2.S <- glm(mk.form(paste0("C2_", g), x1), data=dt, family=binomial)
    return(md.l)
}
# Plots simple models S ~ x1 in model list 'md.l', supposed to contain "nlm.S", "logi.S", "logi2.S"
# Data are also plotted
plot.S.x1 <- function(md.l, xlim = c(0, 600), ylim = c(0, 1), x1 = "Age.of.Death") {
    # ensure that nlm.S, logi.S, logi2.S are contained in md.l
    stopifnot(length(intersect(names(md.l), md.types <- c("nlm.S", "logi.S", "logi2.S"))) == 3)
    # create a list of link functions
    lf <- list()
    lf$nlm.S <-
        function(x) md.l[["nlm.S"]]$coefficients[1] + x * md.l[["nlm.S"]]$coefficients[x1]
    lf$logi.S <-
        function(x) 1 / (1 + exp( - md.l[["logi.S"]]$coefficients[1] - x * md.l[["logi.S"]]$coefficients[x1]) )
    lf$logi2.S <-
        function(x) 0.5 + 0.5 / (1 + exp( - md.l[["logi2.S"]]$coefficients[1] - x * md.l[["logi2.S"]]$coefficients[x1]) )
    # plot data
    fm <- mk.form(as.character(formula( md.l[[ md.types[1] ]] ))[2], x1)
    plot(fm, data=md.l[[1]]$data, pch="+", col="grey", xlim=xlim, ylim=ylim)
    # 'baseline' guide lines
    if (all(ylim == c(0, 1))) {
        abline(h=0.5, lty=3); abline(h=0, lty=3); abline(h=1, lty=3)
    }
    # plot model lines
    sapply(seq_along(lf), function(i)
           plot(lf[[i]], from=xlim[1], to=xlim[2], add=TRUE, col=c(3,4,2)[i], lty=1))
    return(lf)
}
# plot fitted lines for PEG3 using the simple Y ~ age model
g <- genes.or.gsets[1]
par(mfrow=c(1,2), mar=c(4,4,2,1))
invisible(plot.S.x1(md.l <- fit.S.x1(g)))
invisible(plot.S.x1(md.l, xlim=c(0,120), ylim=(c(quantile(md.l[[1]]$data[[ paste0("S_", g) ]],0.05,na.rm=TRUE),1))))

# calculate deviance and AIC
sapply(md.l, deviance)
sapply(md.l, function(x) x$aic)

# plot 4 types of responses pairwise against all explanatory variables
par(mfrow=c(3,4), mar=c(4,4,2,1))
plot(mk.form("S_pool16", expl = expl.var[1:12]), data=d.p$pool16, pch=".")

par(mfrow=c(3,4))
plot(mk.form("R_pool16", expl = expl.var[1:12]), data=d.p$pool16, pch=".")

par(mfrow=c(3,4))
plot(mk.form("S_avg16", expl = expl.var[1:12]), data=d, pch=".")

par(mfrow=c(3,4))
plot(mk.form("R_avg16", expl = expl.var[1:12]), data=d, pch=".")

# Plots response stat Y against explanatory variable x1.
# When 'ar' is TRUE an array of plots is produced presenting genes on separate
# plots; when 'ar' is FALSE data for all genes are shown in a single plot with
# color code
plot.smooth.Y.x1 <- function(Y = "S", ar = TRUE, x1 = "Age.of.Death") {
    bar <- function(g) {
        if(ar) {
            d.x <- d$Age.of.Death
            d.y <- d[[ paste0(Y, "_", g) ]]
        } else {
            d.x <- d.p$pool16$Age.of.Death
            d.y <- d.p$pool16[[ paste0(Y, "_pool16") ]]
        }
        plot(d.x, d.y,
             type="n", xlab = ifelse(ar, "", x1),
             ylab = ifelse(ar, "", Y),
             main = ifelse(ar, g, ""),
             mar = ifelse(ar, c(0,0,0,0), c(5,4,4,2)),
             ylim = c(quantile(d.y, 1e-2, na.rm=TRUE), max(d.y, na.rm=TRUE))
             )
    }
    foo <- function(i) {
        g <- genes[i]
        if(ar) bar(g)
        S.sm <- lowess(age <- d$Age.of.Death,
                       ifelse(is.na(s.g <- d[[ paste0(Y, "_", g) ]]),
                              mean(s.g, na.rm=TRUE), s.g), f = 1/2)
        lines(S.sm, col = col <- ifelse(ar, 1, i + 1))
        points(age, s.g, col = col, pch=".")
    }
    if(ar) {
        par(mfrow=c(4,4), mar=c(2,2,3,1))
        par(cex=1, cex.lab=1)
    }
    else {
        bar("")
    }
    sapply(seq_along(genes), foo)
    return(invisible())
}
plot.smooth.Y.x1(Y = "S", ar = TRUE)

plot.smooth.Y.x1(Y = "R", ar = TRUE)

# condition R_PEG3 vs age on a given institution, RIN value and ancestry
coplot(R_PEG3 ~ Age.of.Death | Institution, data=d)
coplot(R_PEG3 ~ Age.of.Death | DLPFC_RNA_isolation..RIN, data=d)
coplot(R_PEG3 ~ Age.of.Death | Ancestry.EV.2, data=d)

# compare model 1 to model 2; ml1 and ml2 are lists of models fitted to
# different data; mf1 and mf2 are strings of model families ("nlm.R", "logi.S", "logi2.S")
# fun compares the models's AIC values, default is scaled difference
m2m <- function(ml1, ml2, mf1, mf2, fun = function(x, y) (x - y) / abs(x)) {
    sapply(genes.or.gsets, function(r) fun(ml1[[r]][[mf1]]$aic, ml2[[r]][[mf2]]$aic))
}
# test the effect of filtering
rel.aic <- - t(sapply(c("nlm.R","nlm.S", "logi.S", "logi2.S"), function(x) m2m(m, m.nof, x, x)))
n.i <- t(cbind(retained <-
    sapply(genes.or.gsets, function(g) sum(! is.na(m[[g]]$nlm.R$data[[
                                                   paste0("N_", g) ]]))),
             filtered <-
                 sapply(genes.or.gsets,
                        function(g) sum(! is.na(m.nof[[g]]$nlm.R$data[[ paste0("N_", g) ]])))
             ))
par(mar = c(5, 6, 4, 2), mfrow = c(1, 2))
barplot(n.i, horiz=TRUE, beside=FALSE, names.arg=genes.or.gsets,
        las=1, xlab="# of points", main="size of dataset",
        args.legend=list(x = "right"), legend.text=c("retained", "filtered out"))
barplot(rel.aic,
        horiz=TRUE, beside=TRUE, names.arg=genes.or.gsets, las=1, col=c("red", "pink",
                                                               "green",
                                                               "blue"),
        xlab="minus relative AIC", main="improvement by filtering", legend.text=c("nlm.R",
                                                                     "nlm.S",
                                                                     "logi.S",
                                                                     "logi2.S"),
        args.legend=list(x = "topright"))

# barplot presentation of AIC for genes separately
par(mar = c(5, 6, 4, 2), mfrow = c(1, 2))
barplot(sapply(genes, function(r)
               sapply(c("nlm.R", "nlm.S", "logi.S", "logi2.S"), function(mod) m[[r]][[mod]]$aic)),
        horiz=TRUE, beside=TRUE, names.arg=genes, las=1, col=c("red", "pink",
                                                               "green",
                                                               "blue"),
        xlab="AIC", main="gene-wise fit", legend.text=c("nlm.R",
                                                                     "nlm.S",
                                                                     "logi.S",
                                                                     "logi2.S"),
        args.legend=list(x = "topright"))
# barplot presentation of AIC for data aggregated across genes
barplot(sapply(setdiff(genes.or.gsets, genes), function(r)
               sapply(c("nlm.R", "nlm.S", "logi.S", "logi2.S"), function(mod) m[[r]][[mod]]$aic)),
        horiz=TRUE, beside=TRUE, names.arg=setdiff(genes.or.gsets, genes), las=1,
        col=c("red", "pink", "green", "blue"), xlab="AIC", main="fit to
        aggregated data", legend.text=c("nlm.R", "nlm.S", "logi.S", "logi2.S"),
        args.legend=list(x = "bottomright"))

barplot( - t(sapply(c("nlm.R", "nlm.S", "logi.S", "logi2.S"), function(x) m2m(m, m.nof, x, x))),
        beside=TRUE, horiz=TRUE, names.arg=genes.or.gsets, las=1)

# compare logi.S to nlm.R
signif(m2m(m, m, "logi.S", "nlm.R"), 3)
# compare logi2.S to logi.S
signif(m2m(m, m, "logi2.S", "logi.S"), 3)
signif(m2m(m.nof, m.nof, "logi2.S", "logi.S"), 3)

plot.simple.cmp.multiple <- function(g = "PEG3") {
    par(mfrow=c(1,2))
    plot.S.x1(fit.S.x1(g))
    title(main=expression(paste(beta[0], ", ", beta[age], " from simple regr.")))
    plot.S.x1(m[[g]])
    title(main=expression(paste(beta[0], ", ", beta[age], " from multiple. regr.")))
}
plot.simple.cmp.multiple("PEG3") 

plot.simple.cmp.multiple("ZNF331")

# Gets property 'prop' of 'coef' from both md.l1 and md.l2 model lists
# The property 'prop' is an index for the vector c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
coefs.from.md.lists <- function(md.l1, md.l2,
                                prop = 1,
                                coef = "Age.of.Death") {
    md.types <- c("nlm.R", "nlm2.R", "nlm.S", "logi.S", "logi2.S")
    names(md.types) <- md.types
    lapply(md.types,
           function(ty) {
               m1 <- sapply(genes.or.gsets,
                            function(g) summary(md.l1[[g]][[ty]])$coefficients[coef, prop])
               m2 <- sapply(genes.or.gsets,
                            function(g) summary(md.l2[[g]][[ty]])$coefficients[coef, prop])
               cbind(m1, m2)
           })
}
# list of simple regression models
m.smpl <- lapply(genes.or.gsets[1:19], fit.S.x1, d)
for(g in genes.or.gsets[20:22])
    m.smpl[[g]] <- fit.S.x1(g, dt = d.p[[g]])
# get estimates and p-values for the age regression coefficient
betas <- coefs.from.md.lists(m.smpl, m, prop = 1)
pvals <- coefs.from.md.lists(m.smpl, m, prop = 4)
plot.beta.pval <- function(mod = "nlm.S", bts = betas, pvs = pvals,
                           xlim.b=c(-2e-3,2e-3), xlim.logp=c(0,20)) {
    par(mar = c(5, 6, 4, 2), mfrow = c(1, 2))
    barplot(t(bts[[ mod ]]), beside=TRUE, horiz=TRUE, las=1, xlim=xlim.b,
            xlab=expression(hat(beta)[age]), main="estim. effect of age")
    barplot(- log10(t(pvs[[ mod ]])), beside=TRUE, horiz=TRUE, las=1,
            xlim=xlim.logp, xlab=expression(-log[10](p)), main="significance",
            legend.text=c("simple", "multiple"), args.legend=list(x = "right"))
}
plot.beta.pval("nlm.R", xlim.b=c(-1e-0,1e-0), xlim.logp=c(0,15))

plot.beta.pval("nlm.S", xlim.logp=c(0,15))

plot.beta.pval("logi.S", xlim.b=c(-5e-2,5e-2), xlim.logp=c(0,40))

plot.beta.pval("logi2.S", xlim.b=c(-5e-2,5e-2), xlim.logp=c(0,60))

# plot fit lines from simple regression for all 16 genes
par(mfrow=c(4,4), mar=c(4,4,2,1))
lapply(genes, function(g) plot.S.x1(fit.S.x1(g)))
# plot fit lines from multipleiate regression for all 16 genes
par(mfrow=c(4,4), mar=c(4,4,2,1))
lapply(genes, function(g) plot.S.x1(fit.S.x1(g), xlim=c(0,120), ylim=(c(quantile(fit.S.x1(g)[[1]]$data[[ paste0("S_", g) ]],0.05,na.rm=TRUE),1))))

par(mfrow=c(4,4), mar=c(4,4,2,1))
lapply(genes, function(g) plot.S.x1(m[[g]]))

# further conditioning on RIN.2 is not necessary because of its high correlation with RIN
signif(with(d, cor(DLPFC_RNA_isolation..RIN, DLPFC_RNA_isolation..RIN.2)), 3)

}

################################# TODO #################################
