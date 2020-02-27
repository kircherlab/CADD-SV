#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

sink(file=NULL)
load("model.RData")
x=as.matrix(read.table(args[1],header=T))

for(i in 1:108){
  x[,i]=as.numeric(x[,i])
}

x[is.na(x)] <- 0
x=apply(x,2, as.numeric)
x[is.na(x)] <- 0


data=as.data.frame(x[,4:108])

pred=predict(model, data,type="response")

write.table(cbind(x[,1:3],pred),args[2],sep="\t",
            row.names = F, 
            col.names = F, 
            quote = F)
