#!/usr/bin/env Rscript

#Scoring CADD-SV
options(scipen = 999)

args = commandArgs(trailingOnly=TRUE)

name.caddsv=args[1]
cadd.sv.read=function(x,z){
  
  y=as.matrix(read.table(paste(x,"/matrix",z,".bed",sep=""),header=T))
  
  header=colnames(x)
  
  y[,1]=gsub('[chr\n]', '', y[,1])
  
  for(i in 1:131){
    y[,i]=as.numeric(y[,i])
  }
  
  y[is.na(y)] <- 0
  y=apply(y,2, as.numeric)
  y[is.na(y)] <- 0
  header=colnames(y)
  y=data.frame(y)
  
  ylen=y[,3]-y[,2]
  y$DNase.seq_max[which(y$DNase.seq_max>3000)]=3000 #DNAse outliers
  y$DDD_HaploInsuf[which(y$DDD_HaploInsuf>1)]=1 #DDD outliers
  
  tolog=c("EP_distance","A549_nested_dist","A549_tad_dist"
          ,"caki2_nested_dist","caki2_tad_dist","escTAD_distance"
          ,"microsyn_distance","exon_dist","gene_dist"
          ,"start_codon_dist","CADD_sum","CADD_count","gerp_count","PhastCons100_sum"
          ,"PhastCons30_sum","PhastCons20_sum","DI_min"
          ,"DI_max","DNase.seq_sum","H2AFZ_sum"
          ,"H3K27ac_sum","H3K27me3_sum","H3k36me3_sum"
          ,"H3K4me1_sum","H3K4me2_sum","H3K4me3_sum"
          ,"H3K79me2_sum","H3K9ac_sum","H3K9me3_sum"
          ,"H4k20me1_sum","totalRNA.seq_sum","LINSIGHT"
          ,"exon","transcript","gene","x3utr","x5utr","cds","nr_uc_bases")
  
  for(k in tolog){ #logs of distances and sums and counts
    i=which(header==k) 
    y[,i]=round(log10(abs(y[,i])+1),4)
  }
  
  return(y)
}

caddsv=function(x){
  k=list()
  y=list()
  k[[1]]=cadd.sv.read(x,z="")
  k[[2]]=cadd.sv.read(x,z="_100bpup")
  k[[3]]=cadd.sv.read(x,z="_100bpdown")
  z1=match(k[[3]][,3]-100,k[[1]][,3])
  z2=match(k[[3]][,3]-101,k[[1]][,3])
  z=cbind(z1,z2)
  z[is.na(z)] <- 0
  t=apply(z,1,max)
  
  k[[5]]=k[[3]][t,]##reordering sorted 100bpdown
  k[[4]]=k[[2]]+k[[5]]
  k[[4]][,114]=apply(cbind(k[[2]][,114],k[[3]][,121]),1,min)
  k[[4]][,115]=apply(cbind(k[[2]][,115],k[[3]][,122]),1,min)
  k[[4]][,116]=apply(cbind(k[[2]][,116],k[[3]][,123]),1,min)
   y[[1]]=k[[1]]
  y[[2]]=k[[4]]
  return(y)
}


gnomad.scale=readRDS("models/2021scale.RDS")
gnomad.rank=readRDS("models/2021gnomadrank20k.RDS")
gnomad.rank2=readRDS("models/2021gnomadrank20k_2.RDS")

cdel.model=readRDS("models/cdelmodelRF.RDS")
hdel.model=readRDS("models/hdelmodelRF.RDS")
cins.model=readRDS("models/cinsmodelRF.RDS")
hins.model=readRDS("models/hinsmodelRF.RDS")
#contains models and scaling parameters


tbs=caddsv(name.caddsv) #to be scored
id=as.matrix(read.table(paste("input/id_",name.caddsv,".bed",sep=""))) #DEL INS DUP identifyer



DEL=list()
DEL[[1]]=as.matrix(tbs[[1]][which(id[,4]=="DEL"),])
DEL[[2]]=as.matrix(tbs[[2]][which(id[,4]=="DEL"),])
INS=list()
INS[[1]]=as.matrix(tbs[[1]][which(id[,4]=="INS"),])
INS[[2]]=as.matrix(tbs[[2]][which(id[,4]=="INS"),])
DUP=list()
DUP[[1]]=as.matrix(tbs[[1]][which(id[,4]=="DUP"),])
DUP[[2]]=as.matrix(tbs[[2]][which(id[,4]=="DUP"),])

##fixing single DUP DEL INS transformation bug
#if(dim(DEL[[1]])[2]==1){
#DEL[[1]]=t(DEL[[1]])
#DEL[[2]]=t(DEL[[2]])
#}
#if(dim(INS[[1]])[2]==1){
# INS[[1]]=t(INS[[1]])
#  INS[[2]]=t(INS[[2]])
#}
#if(dim(DUP[[1]])[2]==1){
#  DUP[[1]]=t(DUP[[1]])
#  DUP[[2]]=t(DUP[[2]])
#}


dels=DEL
inss=INS
dups=DUP
#scaling according to gnomad healthy population distribution

for(i in 4:131){
  dels[[1]][,i]=scale(c(DEL[[1]][,i]),center=gnomad.scale[[1]][[i]][[2]],scale=gnomad.scale[[1]][[i]][[3]])
  dels[[2]][,i]=scale(c(DEL[[2]][,i]),center=gnomad.scale[[2]][[i]][[2]],scale=gnomad.scale[[2]][[i]][[3]])
  inss[[1]][,i]=scale(c(INS[[1]][,i]),center=gnomad.scale[[3]][[i]][[2]],scale=gnomad.scale[[3]][[i]][[3]])
  inss[[2]][,i]=scale(c(INS[[2]][,i]),center=gnomad.scale[[4]][[i]][[2]],scale=gnomad.scale[[4]][[i]][[3]])
  dups[[1]][,i]=scale(c(DUP[[1]][,i]),center=gnomad.scale[[5]][[i]][[2]],scale=gnomad.scale[[5]][[i]][[3]])
  dups[[2]][,i]=scale(c(DUP[[2]][,i]),center=gnomad.scale[[6]][[i]][[2]],scale=gnomad.scale[[6]][[i]][[3]])
  
}

dels[[1]][is.na(dels[[1]])] <- 0
dels[[2]][is.na(dels[[2]])] <- 0
inss[[1]][is.na(inss[[1]])] <- 0
inss[[2]][is.na(inss[[2]])] <- 0
dups[[1]][is.na(dups[[1]])] <- 0
dups[[2]][is.na(dups[[2]])] <- 0

#scoring DEL INS and DUPs according to the appropriate model
ranker=function(x,gnomad) {
  k=x
  for(i in 1:length(k))
  {
    k[i]=length(which(gnomad<x[i]))/length(gnomad)
  }
  return(k)
}

library(randomForest)
if(dim(dels[[1]])[1]>0){
  del1=predict(cdel.model,dels[[1]][,4:131])*(-1)
  del2=predict(hdel.model,dels[[2]][,4:131])*(-1)
  rank.del1=ranker((-1)*del1,gnomad.rank[[1]])
  rank.del2=ranker((-1)*del2,gnomad.rank[[2]])
  del=apply(cbind(rank.del1,rank.del2),1,min)
  rank.del=1-ranker(del,gnomad.rank2[[1]])
  
  
  
}
if(dim(inss[[1]])[1]>0){
  ins1=predict(hins.model,inss[[1]][,4:131])*(-1)
  ins2=predict(cins.model,inss[[2]][,4:131])*(-1)
  rank.ins1=ranker((-1)*ins1,gnomad.rank[[3]])*(-1)
  rank.ins2=ranker((-1)*ins2,gnomad.rank[[4]])*(-1)
  ins=apply(cbind(rank.ins1,rank.ins2),1,min)*(-1)
  rank.ins=1-ranker(ins,gnomad.rank2[[2]])
  
  
  
}
if(dim(dups[[1]])[1]>0){
  dup1=predict(cdel.model,dups[[1]][,4:131])*(-1)
  dup2=predict(hdel.model,dups[[2]][,4:131])*(-1)
  rank.dup1=ranker((-1)*dup1,gnomad.rank[[5]])*(-1)
  rank.dup2=ranker((-1)*dup2,gnomad.rank[[6]])*(-1)
  dup=apply(cbind(rank.dup1,rank.dup2),1,min)*(-1)
  rank.dup=1-ranker(dup,gnomad.rank2[[3]])

}





##fixing single DUP DEL INS transformation bug

dels2=data.frame(dels[[1]])
inss2=data.frame(inss[[1]])
dups2=data.frame(dups[[1]])

if(dim(id)[2]==5){
  z.del=id[which(id[,4]=="DEL"),5]
  z.ins=id[which(id[,4]=="INS"),5]
  z.dup=id[which(id[,4]=="DUP"),5]
} else {
  z.del=paste(name.caddsv,paste(id[,1],paste(as.numeric(id[,2]),as.numeric(id[,3]),sep="-"),sep=":"),sep="_")[which(id[,4]=="DEL")]
  z.ins=paste(name.caddsv,paste(id[,1],paste(as.numeric(id[,2]),as.numeric(id[,3]),sep="-"),sep=":"),sep="_")[which(id[,4]=="INS")]
  z.dup=paste(name.caddsv,paste(id[,1],paste(as.numeric(id[,2]),as.numeric(id[,3]),sep="-"),sep=":"),sep="_")[which(id[,4]=="DUP")]

}

header=c(colnames(tbs[[1]])[1:3],"type","name","CADDSV-score","raw-score-span","raw-score-flank",colnames(dels[[1]])[4:131])


if(dim(dels[[1]])[1]>0){
d=cbind(dels2[,1:3],"DEL",z.del,rank.del,del1,del2,dels2[,4:131])
colnames(d)=header
d[which(d[,1]==0),1]="X"

}else(d=data.frame())

if(dim(dups[[1]])[1]>0){
dupli=cbind(dups2[,1:3],"DUP",z.dup,rank.dup,dup1,dup2,dups2[,4:131])
colnames(dupli)=header
dupli[which(dupli[,1]==0),1]="X"
}else(dupli=data.frame())

if(dim(inss[[1]])[1]>0){
inserts=cbind(inss2[,1:3],"INS",z.ins,rank.ins,ins1,ins2,inss2[,4:131])
colnames(inserts)=header
inserts[which(inserts[,1]==0),1]="X"
}else(inserts=data.frame())



write.table(rbind(d,inserts,dupli),paste("output/",name.caddsv,".score",sep=""),sep="\t",
            row.names = F, 
            col.names = F, 
            quote = F)

