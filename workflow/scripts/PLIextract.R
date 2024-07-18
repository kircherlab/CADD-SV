#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)
##Annotation Script for SV-human-and-chimp dataset for PLI file



gn2=read.table(args[1],sep="\t",header=F)
pli=read.table(args[2],sep=",",header=T)

gn=gn2[,4]
#gn=droplevels(gn)
gn=as.character(gn)
gene=as.vector(pli[,1])
pli=pli[,2]

score=list() ## gene to PLI score asignment
for(i in which(gn>0)){
  a=strsplit(gn[i],split=",")
  score[[i]]=list()
  for(j in 1:length(a[[1]])){
    score[[i]][[j]]=pli[which(gene==a[[1]][j])]
  }
}

res=vector()
if (!is.null(unlist(score))){
  for(i in which(as.numeric(summary(score)[,1])>0)){
    if(length(unlist(score[[i]]))>0){
      res[i]=max(unlist(score[[i]]))
    }else(res[i]=0)
  }
}

a1=which(res>0)
a2=1:length(gn)
a3=setdiff(a2,a1) ## vector of empty entries -> non genic / no Pli

res[a3]=0
res2=cbind(gn2[,(-4)],res)

write.table(res2,args[3],sep="\t",
            row.names = F, 
            col.names = F, 
            quote = F)


