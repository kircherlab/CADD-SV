#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

sink(file=NULL)


cadd.sv.read=function(x,z){
  
y=as.matrix(read.table(paste(x,"/matrix",z,".bed",sep=""),header=T))

header=colnames(x)

y[,1]=gsub('[chr\n]', '', y[,1])

for(i in 1:130){
  y[,i]=as.numeric(y[,i])
}

y[is.na(y)] <- 0
y=apply(y,2, as.numeric)
y[is.na(y)] <- 0

ylen=y[,3]-y[,2]
y[which(y[,43]>3000),44]=3000 #DNAse outliers
y[which(y[,126]>1),126]=1 #DDD outliers

for(i in c(73,98,102,106,110,114,118,121,122,123,5,7,9,11,42,43,45,47,49,51,53,55,57,59,61,63,65,67,69,131,85,86,87)){ #logs of distances and sums and counts
  y[,i]=round(log10(abs(y[,i])+1),4)
}

y=y[,-c(73,98,102,106,110,114,118)]

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

##############################################################################################################

load("model.RData")
set1=caddsv(args[1]) #read annotation file

f=dim(set1)[2] #features -3

predict1.cdel=predict(cdel.model,set1[[1]][,4:f])
predict2.hdel=predict(hdel.model,set1[[2]][,4:f])
pred=apply(cbind(predict1.cdel,predict2.hdel),1,min)

##############################################################################################################

ranker=function(x,gnomad) {
  k=x
  for(i in 1:length(k))
  {
    k[i]=length(which(gnomad<x[i]))/length(gnomad)
  }
  return(k)
}

rank.set1=ranker(pred,gnomad) #within gnomad scores
perc.rank <- function(x) trunc(rank(x))/length(x)

score=perc.rank(rank.set1)

##############################################################################################################

write.table(cbind(set1[,1:3],score),args[2],sep="\t",
            row.names = F, 
            col.names = F, 
            quote = F)
