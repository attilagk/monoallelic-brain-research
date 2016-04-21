library(plyr)

################################# new functions #################################

# read data and set row names
get.data <- function(filename, id.string='ID', ...) {
    df <- read.delim(filename, na.strings=c('NA', 'ND'), ...)
    row.names(df) <- df[, id.string]
    return(df)
}

# Copies and pastes the the column with 'select.col' from one data frame to
# another, then sets row names
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
    df <- df[sort.list(row.names(df)), ]
    return(df)
}

# replace elements in df.A where df.B <= thrs
filter.obs <- function(df.A, df.B, thrs=50) {
    if(! all( dim(df.A) == dim(df.B))) return(NaN)
    A <- as.vector(as.matrix(df.A))
    B <- as.vector(as.matrix(df.B))
    df <- as.data.frame(matrix(ifelse(B > thrs, A, NA), nrow=nrow(df.A),
                         ncol=ncol(df.A)))
    names(df) <- names(df.A)
    row.names(df) <- row.names(df.A)
    return(df)
}

# transform 'S' statistics into LOI_R based on selected 'genes'
loir.tform <- function(S, genes) {
    loir.rank <- function(x) {
        y <- rank(x, ties.method = "min", na.last = "keep")
        as.numeric(round(y / max(y, na.rm = TRUE) * 100, 0))
    }
    S.mat <- as.matrix(as.data.frame(lapply(S[, genes, drop=FALSE], loir.rank)))
    S.avg <- as.data.frame(apply(S.mat, 1, mean, na.rm=TRUE))
    row.names(S.avg) <- row.names(S)
    names(S.avg) <- "response"
    return(S.avg)
}

# high level function preparing data for fitting (filtering, column pasting)
prepare.fit.data <- function(obs=S.stat, input=expl.var,
                             filter.df.B=tot.read.n, filter.thrs=50,
                             genes=genes[1:8], tform=TRUE,
                             weights=as.data.frame(rep(1, nrow(tot.read.n))))
{
    # filter
    obs <- filter.obs(obs, filter.df.B, thrs=filter.thrs)
    # merge weights into the response (i.e. the transformed obs)
    names(weights) <- "weights"
    row.names(weights) <- row.names(filter.df.B)
    if(tform)
        tform.fun <- function(x) loir.tform(x, genes)
    else
        tform.fun <- identity
    response <- tform.fun(obs)
    response <- merge(response, weights, by="row.names")
    response <- response[-1]
    row.names(response) <- row.names(S.stat)
    # merge with expl.var
    fit.data <- paste.column(response, expl.var,
                             select.col=c("response", "weights"), row.n="RNAseq_ID")
    return(fit.data)
}

# fit glm for a selected 'gene'
do.fit <- function(gene, family=gaussian,
                   tform=TRUE, weights=as.data.frame(rep(1, nrow(tot.read.n)))) {
    dat <- prepare.fit.data(obs=S.stat, input=expl.var, filter.df.B=tot.read.n,
                            filter.thrs=50, genes=gene, tform=tform,
                            weights=weights)
    glm(formula=lin.pred, family=family, data=dat, weights=weights)
}

################################# analysis #################################

# get data
S.stat <- get.data("pop_skew_3June15.txt", sep="\t")
tot.read.n <- get.data("pop_cov_3June15.txt", sep="\t")
tissue.smpl <- get.data("samples.csv", id.string="DLPFC_RNA_ID", sep=",")
expl.var <- get.data("DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv",
                       id.string="DLPFC_RNA_isolation..Sample.RNA.ID", sep="\t")

# clean up explanatory variables
expl.var <- paste.column(tissue.smpl, expl.var, select.col=c("RNAseq_ID"), row.n="RNAseq_ID")
expl.var <- expl.var[sort.list(row.names(expl.var)), ]
rm(tissue.smpl)

# clean observed (response) variables
S.stat <- clean.obs(S.stat)
tot.read.n <- clean.obs(tot.read.n)

# selected genes
genes <- c("PEG3", "INPP5F", "SNRPN", "PWAR6", "ZDBF2", "MEG3", "ZNF331",
           "GRB10", "PEG10", "SNHG14", "NAP1L5", "KCNQ1OT1", "MEST" )
names(genes) <- genes

# Linear predictor.  Note that not all available explanatory variables are
# selected (following Ifat's work)
lin.pred <- (response ~ `Age.of.Death` + Institution + Gender +
             `PMI..in.hours.` + Dx + `DLPFC_RNA_isolation..RIN` +
             `DLPFC_RNA_isolation..RIN.2` +
             `DLPFC_RNA_report..Clustered.Library.Batch` + `Ancestry.EV.1` +
             `Ancestry.EV.2` + `Ancestry.EV.3` + `Ancestry.EV.4` +
             `Ancestry.EV.5` )

# Reproducing Ifat's previous results.
# Should match closely but not precisely due to my omission of a rounding step
data.loir.8 <- prepare.fit.data(obs=S.stat, input=expl.var, filter.df.B=tot.read.n, filter.thrs=50, genes=genes[1:8])
nlm.loir.8 <- glm(formula=lin.pred, family=gaussian, data=data.loir.8)

nlm.loir <- lapply(genes, do.fit, family=gaussian)

################################# Ifat's code #################################

transform.data <- function(genes13, genes.subset=genes13[1:8]){
genes <- genes13

samples=read.csv("samples.csv", header = TRUE)
S = read.table("pop_skew_3June15.txt", header = TRUE, na.strings="ND", sep="\t")
covariates=read.csv("DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv", sep="\t",header = TRUE)
C = read.table("pop_cov_3June15.txt", header = TRUE, na.strings="ND", sep="\t")

row.names(C)<-C[,"ID"]

row.names(covariates)<-covariates[,"DLPFC_RNA_isolation..Sample.RNA.ID"] 
row.names(samples)<-samples$RNAseq_ID
row.names(S)<-S[,"ID"]

S0<-arrange(S,S$age) 
age<-S0$age
S1=S0[,2:18]
row.names(S1)<-S0[,"ID"]

C0<-arrange(C,C$age)
C1=C0[,2:18]
row.names(C1)<-C0[,"ID"]
C2<-matrix(data = NA, nrow = 579, ncol = 13)
colnames(C2)<-genes
row.names(C2)<-row.names(C1)

THR=50
S2<-matrix(data = NA, nrow = 579, ncol = 13)
colnames(S2)<-genes
row.names(S2)<-row.names(S1)

for (g in 1:13){
  gene=genes[g]
  gene_cov<-C1[,gene]
  THR_loc<-which(gene_cov>THR, arr.ind=TRUE)
  S2[THR_loc,gene]<-S1[THR_loc,gene]
  C2[THR_loc,gene]<-C1[THR_loc,gene]
}
S2<-data.frame(S2)
C2<-data.frame(C2)

genes <- genes.subset

R<-matrix(data = NA, nrow = 579, ncol = length(genes))
colnames(R)<-genes
for (g in seq_along(genes)){
  gene=genes[g]
  gene_data<-as.numeric(round(S2[,gene],3))
  gene_rank<-rank(gene_data,ties.method = "min",na.last = "keep")
  R[,gene]<-gene_rank
}

PR<-matrix(data = NA, nrow = 579, ncol = length(genes))
colnames(PR)<-genes
for (g in seq_along(genes)){
  gene=genes[g]
  rank_max<-max(R[,gene], na.rm = TRUE)
  for (s in 1:579){
    if (!is.na(R[s,gene])){
      PR[s,gene]<-as.numeric(round(R[s,gene]/rank_max*100,0))
    }
  }
} 


LOI_R<-matrix(data = NA, nrow = 579, ncol = 1)
for (s in 1:579){
  LOI_R[s,1]<-as.numeric(sprintf("%.0f",mean(PR[s,],na.rm = TRUE)))
}

genes <- genes13

QUA<-matrix(data = NA, nrow = 1, ncol = 13)
colnames(QUA)<-genes
for (g in 1:13){
  gene=genes[g]
  Tgene<-S2[,gene]
  Mgene<-length(which(!is.na(Tgene)))
  quar15<-round(Mgene*0.15)
  sortTgene<-sort(Tgene,na.last=TRUE)
  QUA[1,gene]<-sortTgene[quar15]
}
QUA<-as.matrix(QUA)

  ## transfer to 01matrix
X<-matrix(data = NA, nrow = 579, ncol = 13)
colnames(X)<-genes
for (g in 1:13){
  gene=genes[g]
  gene_med=QUA[1,gene]
  gene_skew<-S2[,gene]
  loc_NA<-which(is.na(gene_skew), arr.ind=TRUE)
  loc1<-which(gene_skew<gene_med, arr.ind=TRUE)
  loc0<-which(gene_skew>=gene_med, arr.ind=TRUE) 
  X[loc_NA,gene]<-NA
  X[loc1,gene]<-1
  X[loc0,gene]<-0
}
X<-data.frame(X)
row.names(X)<-row.names(S2)

X1<-cbind(LOI_R,X)
X2<-merge(samples,X1, by="row.names",all=FALSE)
row.names(X2)<-X2$DLPFC_RNA_ID
X3<-X2[,21:34]
row.names(X3)<-X2$DLPFC_RNA_ID
FULL<-merge(X3,covariates,by="row.names",all=FALSE)

return(FULL)
}

fit.lm <- function(FULL, do.thrs=TRUE, age.thrs=70) {
    if(do.thrs){
        glm(LOI_R ~ (`Age.of.Death` > age.thrs) + 
                      Institution+ Gender   + `PMI..in.hours.`+ Dx+ 
                      `DLPFC_RNA_isolation..RIN`+  `DLPFC_RNA_isolation..RIN.2`+ 
                      `DLPFC_RNA_report..Clustered.Library.Batch` + 
                      `Ancestry.EV.1`+ `Ancestry.EV.2`+ `Ancestry.EV.3`+ `Ancestry.EV.4`+ `Ancestry.EV.5`, 
                    data=FULL)
    } else {
        glm(LOI_R ~ `Age.of.Death` + 
                      Institution+ Gender   + `PMI..in.hours.`+ Dx+ 
                      `DLPFC_RNA_isolation..RIN`+  `DLPFC_RNA_isolation..RIN.2`+ 
                      `DLPFC_RNA_report..Clustered.Library.Batch` + 
                      `Ancestry.EV.1`+ `Ancestry.EV.2`+ `Ancestry.EV.3`+ `Ancestry.EV.4`+ `Ancestry.EV.5`, 
                    data=FULL)
    }
}

# m8 <- fit.lm(transform.data(genes13, genes13[1:8]))
# m13 <- fit.lm(transform.data(genes13, genes13))

