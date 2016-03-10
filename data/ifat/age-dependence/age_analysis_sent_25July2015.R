##############################
### Age analysis
##############################
#install.packages('gplots')
library(plyr)
library(RColorBrewer)
library(gplots)
library(ggplot2)
## read data
S = read.table("/Users/Ifat/Dropbox/writingRdata/pop_skew_3June15.txt", header = TRUE, na.strings="ND", sep="\t")
C = read.table("/Users/Ifat/Dropbox/writingRdata/pop_cov_3June15.txt", header = TRUE, na.strings="ND", sep="\t")
samples=read.csv("/Users/Ifat/Dropbox/writingRdata/samples.csv", header = TRUE)
H=read.table("/Users/Ifat/Dropbox/writingRdata/ID_hap_new1.txt", header = TRUE, na.strings="ND", sep="\t")
covariates=read.csv("/Users/Ifat/Dropbox/writingRdata/DLPFC.ensembl.KNOWN_AND_SVA.ADJUST.SAMPLE_COVARIATES.tsv", sep="\t",header = TRUE)
row.names(covariates)<-covariates[,"DLPFC_RNA_isolation..Sample.RNA.ID"] 
row.names(samples)<-samples$RNAseq_ID
  
row.names(S)<-S[,"ID"]
row.names(C)<-C[,"ID"]
S0<-arrange(S,S$age) 
C0<-arrange(C,C$age)
age<-S0$age
S1=S0[,2:18]
row.names(S1)<-S0[,"ID"]
C1=C0[,2:18]
row.names(C1)<-C0[,"ID"]
#genes=colnames(S1[,2:17])
#genes=sort(genes)
genes<-c("PEG3"     ,"INPP5F" ,  "SNRPN" , "PWAR6" , "ZDBF2" ,
        "MEG3"  ,  "ZNF331"   ,  "GRB10"    , "PEG10"  ,   "SNHG14"  ,       
        "NAP1L5"  ,  "KCNQ1OT1" ,"MEST" )
ID<-samples[,c("ID","RNAseq_ID","Dx")]

## filter by threshold
THR=50
S2<-matrix(data = NA, nrow = 579, ncol = 13)
C2<-matrix(data = NA, nrow = 579, ncol = 13)
colnames(S2)<-genes
colnames(C2)<-genes
row.names(C2)<-row.names(C1)
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


## scatter plot age skew
PDFname=paste("/Users/Ifat/Dropbox/writingRdata/age_skew_15July15_",THR,".pdf",sep="")
pdf(file=PDFname,width=25,height=15)
par(mfrow=c(3,5))
JPEGname=paste("/Users/Ifat/Dropbox/writingRdata/age_skew_15July15_",THR,".jpg",sep="")
jpeg(file=JPEGname,width=25,height=15)
N=13
for (i in 1:N ) {
  gene=genes[i]
  max_cov<-max(C2,na.rm = TRUE)
  x<-age
  y<-S2[,gene]
  col<-round((1-C2[,gene]/max_cov)*5)
  col[is.na(col)] <- 0
  col<-round(col,digits=0)
  z<-rgb(5-(col/2),col,col/2,4,max = 5)
  DF<-data.frame(x,y,z,stringsAsFactors = FALSE)
  #cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5
  with (DF,plot (x,y,ylim=c(0.5,1.05),main=gene, cex=1.5, cex.main=2.5, ylab="Skewness", xlab="Age",xlim=c(10,110),pch = 21,col="black"))
}
dev.off()

####################
## rank
####################
genes<-c("PEG3","INPP5F","SNRPN","PWAR6","ZDBF2","MEG3","ZNF331","GRB10")

R<-matrix(data = NA, nrow = 579, ncol = 8)
colnames(R)<-genes
for (g in 1:8){
  gene=genes[g]
  gene_data<-as.numeric(round(S2[,gene],3))
  gene_rank<-rank(gene_data,ties.method = "min",na.last = "keep")
  R[,gene]<-gene_rank
}
PR<-matrix(data = NA, nrow = 579, ncol = 8)
colnames(PR)<-genes
for (g in 1:8){
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
table(LOI_R)
## heatmap
label<-matrix(data = NA, nrow = 579, ncol = 1)
for (i in 0:10){
  d<-i*80+1
  label[d]<-age[d]
}
#my_palette <- colorRampPalette(c("black", "yellow"))(n = 100)

PDFname=paste("/Users/Ifat/Dropbox/writingRdata/age_image_9June15_",THR,".pdf",sep="")
title=paste("Age (coverage>",THR,")")
pdf(file=PDFname,width=10,height=15)
heatmap.2(R,
          scale = c("column"),
          Colv=FALSE,
          Rowv = FALSE,
          # cex.names=6.0
          dendrogram = c("none"),
          #col=my_palette,
          col=greenred(100),
          na.rm = TRUE,
          trace = "none",
          density.info = "none",
          cexRow=1.5,
          cexCol=1.1,
          labRow=label,
          offsetRow=0,
          offsetCol=-0.2,
          keysize =1,
          main = title)

dev.off()

## LOI
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
LOI<-rowSums(X,na.rm = TRUE)
## end LOI

LOI1<-matrix(data = NA, nrow = 579, ncol = 1)
LOI0<-matrix(data = NA, nrow = 579, ncol = 1)
LOIR<-matrix(data = NA, nrow = 579, ncol = 1)
for (i in 1:579){
  LOI1[i,1]<-sum(X[i,1:13]==1,na.rm = TRUE)
  LOI0[i,1]<-sum(X[i,1:13]==0,na.rm = TRUE)
  if (LOI0[i,1]>0||LOI1[i,1]){
    Ratio<-LOI1[i,1]/(LOI1[i,1]+LOI0[i,1])
    LOIR[i,1]<-Ratio 
  }
}


## GLM LOI
X1<-cbind(LOI_R,X)
X2<-merge(samples,X1, by="row.names",all=FALSE)
row.names(X2)<-X2$DLPFC_RNA_ID
X3<-X2[,21:34]
row.names(X3)<-X2$DLPFC_RNA_ID
FULL<-merge(X3,covariates,by="row.names",all=FALSE)
summary(glm(LOI_R ~ (`Age.of.Death`>70)+ 
              Institution+ Gender   + `PMI..in.hours.`+ Dx+ 
              `DLPFC_RNA_isolation..RIN`+  `DLPFC_RNA_isolation..RIN.2`+ 
              `DLPFC_RNA_report..Clustered.Library.Batch` + 
              `Ancestry.EV.1`+ `Ancestry.EV.2`+ `Ancestry.EV.3`+ `Ancestry.EV.4`+ `Ancestry.EV.5`, 
            data=MSSM))
summary(glm(LOI_R ~ `Age.of.Death`+ 
              Institution+ Gender   + `PMI..in.hours.`+ Dx+ 
              `DLPFC_RNA_isolation..RIN`+  `DLPFC_RNA_isolation..RIN.2`+ 
              `DLPFC_RNA_report..Clustered.Library.Batch` + 
              `Ancestry.EV.1`+ `Ancestry.EV.2`+ `Ancestry.EV.3`+ `Ancestry.EV.4`+ `Ancestry.EV.5`, 
            data=FULL))
MSSM<-subset(FULL,Institution=="MSSM")

summary(glm(LOI_R ~ `Age.of.Death`+ 
               Gender   + `PMI..in.hours.`+ Dx+ 
              `DLPFC_RNA_isolation..RIN`+  `DLPFC_RNA_isolation..RIN.2`+ 
              `DLPFC_RNA_report..Clustered.Library.Batch` + 
              `Ancestry.EV.1`+ `Ancestry.EV.2`+ `Ancestry.EV.3`+ `Ancestry.EV.4`+ `Ancestry.EV.5`, 
            data=MSSM))
summary(glm(LOI_R ~ (`Age.of.Death`>70)+ 
              Gender   + `PMI..in.hours.`+ Dx+ 
              `DLPFC_RNA_isolation..RIN`+  `DLPFC_RNA_isolation..RIN.2`+ 
              `DLPFC_RNA_report..Clustered.Library.Batch` + 
              `Ancestry.EV.1`+ `Ancestry.EV.2`+ `Ancestry.EV.3`+ `Ancestry.EV.4`+ `Ancestry.EV.5`, 
            data=MSSM))

## fisher's exact test
result<-matrix(1:39,ncol=3,nrow=13)
old_data<-subset(S2,age>=68,na.rm = TRUE)
young_data<-subset(S2,age<68,na.rm = TRUE)
F<-matrix(data = NA, nrow = 4, ncol = 13)
colnames(F)<-genes
for (g in 1:13 ) {
  gene=genes[g]
  gene_med=QUA[1,gene]  
  F[1,g]<-length(subset(old_data[,gene],old_data[,gene]>gene_med))
  F[2,g]<-length(subset(old_data[,gene],old_data[,gene]<=gene_med))
  F[3,g]<-length(subset(young_data[,gene],young_data[,gene]>gene_med))
  F[4,g]<-length(subset(young_data[,gene],young_data[,gene]<=gene_med))
  O<-c(F[1,g],F[2,g])
  Y<-c(F[3,g],F[4,g])
  x<-data.frame(O,Y)
  subjects<-F[1,g]+F[2,g]+F[3,g]+F[4,g]
  fisher_val<-fisher.test(x)  
  p<-format(fisher_val$p.value,scientific=TRUE,digits=3)
  
  result[g,1]<-gene
  result[g,2]<-p
  result[g,3]<-subjects
}

## hist and plot
PDFname=paste("/Users/Ifat/Dropbox/writingRdata/age_hist_LOI_9June15_",THR,".pdf",sep="")
pdf(file=PDFname,width=7,height=7)
OLD<-subset(X, age>68)
YOUNG<-subset(X,age<=68)
O<-rowSums(OLD,na.rm = TRUE)
Y<-rowSums(YOUNG,na.rm = TRUE)

bins<-c(0:14)
test<-c(0:13)
ybins=c(0,2,4,16)
ao<-hist(O,breaks=bins,plot=FALSE,right = FALSE,include.lowest=FALSE)
aao<-log2(ao$counts)
plot(c(1:14),aao, type="b", col="red", lwd=2, ylab="log2(subjects)",ylim=c(0,8),xlab="number of genes with LOI")
ay<-hist(Y,breaks=bins,plot=FALSE,right = FALSE,include.lowest=FALSE)
aay<-log2(ay$counts)
points(c(1:14),aay, type="b", col="black", lwd=2, ylab="log2(subjects)",xlab="number of genes with LOI")
legend(11.2, 8, c("Young","Old"), lty=c(1,1), # gives the legend appropriate symbols (lines)      
       lwd=c(2,2),col=c("black","red")) # gives the legend lines the correct color and width)

dev.off()


