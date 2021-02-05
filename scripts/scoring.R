#!/usr/bin/env Rscript

#Scoring CADD-SV


args = commandArgs(trailingOnly=TRUE)

name.caddsv=args[1]

#######################################################################################################################
#functions and loading models/scaling parameters of gnomad/ rank distribution of gnomad
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
  k[[2]]=cadd.sv.read(x,z="_100kbup")
  k[[3]]=cadd.sv.read(x,z="_100kbdown")
  k[[4]]=k[[2]]+k[[3]]
  k[[4]][,114]=apply(cbind(k[[2]][,114],k[[3]][,114]),1,min)
  k[[4]][,115]=apply(cbind(k[[2]][,115],k[[3]][,115]),1,min)
  k[[4]][,116]=apply(cbind(k[[2]][,116],k[[3]][,116]),1,min)
  y[[1]]=k[[1]]
  y[[2]]=k[[4]]
  return(y)
}


gnomad.scale=readRDS("/dependencies/2021scale.rds")
gnomad.rank=readRDS("/dependencies/2021gnomadrank.RDS")
cdel.model=readRDS("/dependencies/2021cdelmodel.rds")
hdel.model=readRDS("/dependencies/2021hdelmodel.rds")
cins.model=readRDS("/dependencies/2021cinsmodel.rds")
hins.model=readRDS("/dependencies/2021hinsmodel.rds")
 #contains models and scaling parameters


tbs=caddsv(name.caddsv) #to be scored
id=as.matrix(read.table(paste("../beds/",name.caddsv,"_id.bed",sep=""))) #DEL INS DUP identifyer

#######################################################################################################################
#splitting according to SV-type

DEL=list()
DEL[[1]]=as.matrix(tbs[[1]][which(id[,4]=="DEL"),])
DEL[[2]]=as.matrix(tbs[[2]][which(id[,4]=="DEL"),])
INS=list()
INS[[1]]=as.matrix(tbs[[1]][which(id[,4]=="INS"),])
INS[[2]]=as.matrix(tbs[[2]][which(id[,4]=="INS"),])
DUP=list()
DUP[[1]]=as.matrix(tbs[[1]][which(id[,4]=="DUP"),])
DUP[[2]]=as.matrix(tbs[[2]][which(id[,4]=="DUP"),])

#######################################################################################################################
##fixing single DUP DEL INS transformation bug
if(dim(DEL[[1]])[2]==1){
DEL[[1]]=t(DEL[[1]])
DEL[[2]]=t(DEL[[2]])
}
if(dim(INS[[1]])[2]==1){
  INS[[1]]=t(INS[[1]])
  INS[[2]]=t(INS[[2]])
}
if(dim(DUP[[1]])[2]==1){
  DUP[[1]]=t(DUP[[1]])
  DUP[[2]]=t(DUP[[2]])
}


dels=DEL
inss=INS
dups=DUP
#######################################################################################################################
#scaling according to gnomad healthy population distribution

for(i in 4:124){
     dels[[1]][,i]=scale(c(DEL[[1]][,i]),center=gnomad.scale[[1]][[i]][[2]],scale=gnomad.scale[[1]][[i]][[3]])
     dels[[2]][,i]=scale(c(DEL[[2]][,i]),center=gnomad.scale[[2]][[i]][[2]],scale=gnomad.scale[[2]][[i]][[3]])
     inss[[1]][,i]=scale(c(INS[[1]][,i]),center=gnomad.scale[[3]][[i]][[2]],scale=gnomad.scale[[3]][[i]][[3]])
     inss[[2]][,i]=scale(c(INS[[2]][,i]),center=gnomad.scale[[4]][[i]][[2]],scale=gnomad.scale[[4]][[i]][[3]])
     dups[[1]][,i]=scale(c(DUP[[1]][,i]),center=gnomad.scale[[5]][[i]][[2]],scale=gnomad.scale[[5]][[i]][[3]])
     dups[[2]][,i]=scale(c(DUP[[2]][,i]),center=gnomad.scale[[6]][[i]][[2]],scale=gnomad.scale[[6]][[i]][[3]])
     
}
 
#######################################################################################################################
#scoring DEL INS and DUPs according to the appropriate model
library(randomForest)
del1=predict(cdel.model,dels[[1]][,4:124])
del2=predict(hdel.model,dels[[2]][,4:124])
ins1=predict(hins.model,inss[[1]][,4:124])
ins2=predict(cins.model,inss[[2]][,4:124])
dup1=predict(cdel.model,dups[[1]][,4:124])
dup2=predict(cins.model,dups[[2]][,4:124])


#######################################################################################################################
#ranking according to gnomad healthy cohort distribution
ranker=function(x,gnomad) {
  k=x
  for(i in 1:length(k))
  {
    k[i]=length(which(gnomad<x[i]))/length(gnomad)
  }
  return(k)
}

rank.del1=ranker(del1,gnomad.rank[[1]])
rank.del2=ranker(del2,gnomad.rank[[2]])
rank.ins1=ranker(ins1,gnomad.rank[[3]])
rank.ins2=ranker(ins2,gnomad.rank[[4]])
rank.dup1=ranker(dup1,gnomad.rank[[5]])
rank.dup2=ranker(dup2,gnomad.rank[[6]])
#######################################################################################################################
# weighing model scores

del=apply(cbind(rank.del1,rank.del2),1,min)
ins=apply(cbind(rank.ins1,rank.ins2),1,min)
dup=apply(cbind(rank.dup1,rank.dup2),1,min)

#######################################################################################################################
#putting it back together
cadd=cbind(id,2)
cadd[which(id[,4]=="DEL"),5]=del
cadd[which(id[,4]=="DUP"),5]=dup
cadd[which(id[,4]=="INS"),5]=ins

#######################################################################################################################
#output
write.table(cadd,paste(name.caddsv,".score",sep=""),sep="\t",
            row.names = F, 
            col.names = T, 
            quote = F)

