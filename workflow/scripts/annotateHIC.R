#!/usr/bin/env Rscript


args = commandArgs(trailingOnly=TRUE)
##Annotation Script for SV-human-and-chimp dataset for HIC domains in bed format

sub=read.table(args[1],sep="\t")


l=length(sub[,1])

x1=vector()
x2=vector()
x3=vector()
x4=vector()

for(i in 1:l){
  x1[i]=sub[i,5]-sub[i,2]
  x2[i]=sub[i,6]-sub[i,3]
  x3[i]=sub[i,5]-sub[i,3]
  x4[i]=sub[i,6]-sub[i,2]
  
}

xx=cbind(x1,x2,x3,x4)
p1=which(xx[,1]>0)
p2=which(xx[,2]>0)
p3=which(xx[,3]>0)
p4=which(xx[,4]>0)

n1=which(xx[,1]<0)
n2=which(xx[,2]<0)
n3=which(xx[,3]<0)
n4=which(xx[,4]<0)

btad=intersect(intersect(intersect(p1,p2),p3),p4) #before a TAD
atad=intersect(intersect(intersect(n1,n2),n3),n4) #after a TAD
stad=intersect(intersect(intersect(n1,p2),n3),p4) #spanning a TAD
itad=intersect(intersect(intersect(p1,n2),n3),p4) #intra a TAD
lbtad=intersect(intersect(intersect(p1,p2),n3),p4) #left boundary a TAD
rbtad=intersect(intersect(intersect(n1,n2),n3),p4) #right boundary a TAD

intad=rep(0,l)
intad[stad]=1
intad[itad]=1

notad=rep(0,l)
notad[btad]=1
notad[atad]=1

boundarytad=rep(0,l)
boundarytad[lbtad]=1
boundarytad[rbtad]=1

m=vector()
for(i in 1:l){
  m[i]=min(abs(xx[i,]))
}


bed=cbind(as.character(sub[,1]),sub[,2],sub[,3],intad,notad,boundarytad,m) #chr start end binary: intad/ outtad / boundary / distance to boundary

colnames(bed)=c("chr","start","end","in_tad","out_tad","boundary","distance")

write.table(bed,args[2],sep="\t",
            row.names = F, 
            col.names = F, 
            quote = F)

