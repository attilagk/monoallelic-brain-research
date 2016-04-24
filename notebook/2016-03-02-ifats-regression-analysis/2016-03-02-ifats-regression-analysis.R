library(plyr)

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

