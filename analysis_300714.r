what="nothing"
if (what=="format")
 {
 for(i in 1: length(dir(pattern="RPKM")))
 {
  load(dir(pattern="RPKM")[i])
 }
 load("~/smarguer/SAM_DATA/MATOS_revision.rda")
 test=R1807_1_50[,2]
 for(i in 1:12)
 {
  toDo=get(paste("R1807_",i,"_50",sep=''))
  test=cbind(test,toDo[,"RPKM_exonic"])
 }
 test=test[,2:13]
 row.names(test)=R1807_1_50[,1]
 rm(toDo)
 
 test=cbind(test,BIG[,c(8,9,33,34)])
 
 CDC2=test
 
 CDC2N=CDC2
 
 for(i in 1:12)
 {
  CDC2N[,i]=CDC2N[,i]/CDC2[,1]
  CDC2N[,i]=round(CDC2N[,i],2)
  CDC2N[which(CDC2N[,i]=="Inf"),i]=NA 
  CDC2N[which(CDC2N[,i]=="-Inf"),i]=NA 
  CDC2N[which(CDC2N[,i]==0),i]=NA 
 }
 round(CDC2[,1:12],2)
 nCDC2=CDC2
 RNA=read.delim("CDC2_RNA.txt",stringsAsFactors=F,header=F)
 AREA=read.delim("CDC2_AREA.txt",stringsAsFactors=F,header=F)
 
 for(i in 1:12)
 {
  nCDC2[,i]=nCDC2[,i]*RNA[i,1]
 }
 
 nCDC2N=nCDC2
 
 for(i in 1:12)
 {
  nCDC2N[,i]=nCDC2N[,i]/nCDC2[,1]
  nCDC2N[,i]=round(nCDC2N[,i],2)
  nCDC2N[which(nCDC2N[,i]=="Inf"),i]=NA
  nCDC2N[which(nCDC2N[,i]=="-Inf"),i]=NA
  nCDC2N[which(nCDC2N[,i]==0),i]=NA
 }
 
 CDC2[,"common.BIG"]=as.character(CDC2[,"common.BIG"])
 CDC2N[,"common.BIG"]=as.character(CDC2N[,"common.BIG"])
 nCDC2[,"common.BIG"]=as.character(nCDC2[,"common.BIG"])
 nCDC2N[,"common.BIG"]=as.character(nCDC2N[,"common.BIG"])
 colnames(CDC2)[1:12]=c("R0","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11")
 colnames(CDC2N)[1:12]=c("R0","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11")
 colnames(nCDC2)[1:12]=c("R0","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11")
 colnames(nCDC2N)[1:12]=c("R0","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11")
 save(list=c("CDC2","CDC2N","nCDC2","nCDC2N","RNA","AREA"),file="RPKM_CDC2_1807.rda")
}

if(what=="PROT")
{
 conv=read.delim("~/smarguer/SAM_DATA/PomBase2UniProt_061014.TXT2",stringsAsFactors=F,header=F)
 colnames(conv)=c("Systematic","Uniprot")
 P=read.delim("LF_031014.txt", stringsAsFactors=F,row.names=1)
 P=P[,13:(ncol(P)-2)]
 Pn=P
 for(i in 1:12)
 {
  Pn[,i]=P[,i]/P[,1]
 }
 P.raw=P
 P=Pn

 for(i in 1:nrow(P))
 {
  if (row.names(P)[i] %in% conv[,2])
  {
   if(nrow(conv[which(conv[,2]==row.names(P)[i]),])>1)
   {
    print(conv[which(conv[,2]==row.names(P)[i]),])
    print("###")
   }
   else
   {
    row.names(P)[i]=conv[which(conv[,2]==row.names(P)[i]),1]
   }
  }
  else
  {
   print(row.names(P)[i])
  }
 }

 Prot=matrix(NA,nrow(CDC2N),13)
 Prot=data.frame(Prot)
 for(i in 1:nrow(CDC2N))
 {
  #print(i)
  if(row.names(CDC2N)[i] %in% row.names(P))
  {
   if(length(which(row.names(P)==row.names(CDC2N)[i])) > 1)
   {
    print(row.names(P)[which(row.names(P)==row.names(CDC2N)[i])])
    Prot[i,1]=row.names(CDC2N)[i]
    Prot[i,2:13]=rep(NA,12)
   }
   else
   {
    Prot[i,1]=row.names(P)[which(row.names(P)==row.names(CDC2N)[i])]
    Prot[i,2:13]=P[which(row.names(P)==row.names(CDC2N)[i]),1:12]
   }
  }
  else
  {
   Prot[i,1]=row.names(CDC2N)[i]
   Prot[i,2:13]=rep(NA,12)
  }
 }
 CDC2N=cbind(CDC2N,Prot)
 CDC2N=CDC2N[,-17]
 CDC2N=cbind(CDC2N,BIG[,"MM_prot"])
 CDC2N[,17:28]=round(CDC2N[,17:28],2)
 colnames(CDC2N)[1:12]=c("R0","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11")
 colnames(CDC2N)[17:29]=c("P0","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11","MM_prot")
 save(list=c("CDC2","CDC2N","nCDC2","nCDC2N","RNA","AREA"),file="RPKM_CDC2_1807.rda")
}
####
if(what=="stress")
{
 Chen=read.delim("~/smarguer/SAM_DATA/DATA_Chen_2003.txt",stringsAsFactors=F,skip=8)
 Chen1=read.delim("~/smarguer/SAM_DATA/DATA_Chen_2003.txt",stringsAsFactors=F,skip=5)
 colnames(Chen)=colnames(Chen1)
 rm(Chen1)
 heat=data.frame()
 for(i in 1:nrow(CDC2N))
 {
  if(toupper(row.names(CDC2N)[i]) %in% toupper(Chen[,3]))
  {
   print(Chen[which(toupper(Chen[,3])==toupper(row.names(CDC2N)[i])),3])
   heat[i,1]=Chen[which(toupper(Chen[,3])==toupper(row.names(CDC2N)[i])),"Heat"][1]
   heat[i,2]=Chen[which(toupper(Chen[,3])==toupper(row.names(CDC2N)[i])),"Heat.1"][1]
   heat[i,3]=Chen[which(toupper(Chen[,3])==toupper(row.names(CDC2N)[i])),"Heat.2"][1]
   row.names(heat)[i]=Chen[which(toupper(Chen[,3])==toupper(row.names(CDC2N)[i])),3][1]
  }
  else
  {
   heat[i,1]=NA
   heat[i,2]=NA
   heat[i,3]=NA
   row.names(heat)[i]=row.names(CDC2N)[i]
  }
 }
 colnames(heat)=c("H0","H15","H60")
}
####
#load("RPKM_CDC2_1807.rda")
if(what=="slope")
{
 slope=vector()
 rsquared=vector()
 for(i in 1:nrow(nCDC2))
 {
  print(i)
  dat=nCDC2[i,1:12]
  dat=as.numeric(dat)
  dat[which(dat==0)]=NA
 
  if(all(is.na(dat)==T))
  {
   slope[i]=NA
   rsquared[i]=NA
  }
  else
  {
   mod=lm(log2(dat) ~ log2(AREA[,1]))
   rsquared[i]=summary(mod)$adj.r.squared
   slope[i]=mod$coefficient[2]
  }
 }
}

plot.cdc=function(gene,dat.raw=CDC2,dat.norm=CDC2G,dat.scaled=nCDC2,dat.scaled.norm=nCDC2N,type,complete=F,Coef=slope)
{
 heat=c("H0","H15","H60")
 rna=c("R0","R1","R2","R3","R4","R5","R6","R7","R8","R9","R10","R11")
 prot=c("P0","P1","P2","P3","P4","P5","P6","P7","P8","P9","P10","P11")
 tot=c(heat,rna,prot)
 if(length(which(dat.norm[,"common.BIG"]==gene))==0)
 {
  if((length(grep(substr(gene,1,3),dat.norm[,"common.BIG"])) == 1) & (nchar(gene) < 7))
  {
   return(paste("Did you mean: ",dat.norm[grep(substr(gene,1,3),dat.norm[,"common.BIG"]),"common.BIG"],"?",sep=""))
  }
  else if(length((dat.norm[grep(substr(gene,1,3),dat.norm[,"common.BIG"]),"common.BIG"]) > 1) & (nchar(gene) < 7))
  {
   toDo.col=c(heat,rna,prot)
   toDo=dat.norm[grep(substr(gene,1,3),dat.norm[,"common.BIG"]),tot]
   toDo=Zero2Na(toDo)
   if(complete==T)
   {
    toDo=toDo[which(toDo[,"P0"]!=0),]
   }
   heatmap.2(as.matrix(log2(toDo)),Colv=F,breaks=10,symbreaks=T,col=blue2yellow,symkey=T,trace="n",na.rm=T,dendrogram="n",key=F,cexRow=0.8,labRow=dat.norm[grep(substr(gene,1,3),dat.norm[,"common.BIG"]),"common.BIG"])

   print("can't find it, what about one of these:")
   return(dat.norm[grep(substr(gene,1,3),dat.norm[,"common.BIG"]),"common.BIG"])
  }
  else if (grepl("SP",gene)==T)
  {
   gene=dat.norm[which(row.names(dat.norm)==gene),"common.BIG"]
   print(paste("Using: ",gene))
  }
  else if (nrow(dat.norm[grep(gene,dat.norm[,"Annot"]),tot])> 1)
  {
   toDo.col=c(heat,rna,prot)
   toDo=dat.norm[grep(gene,dat.norm[,"Annot"]),tot]
   toDo=Zero2Na(toDo)
   if(complete==T)
   {
    toDo=toDo[which(toDo[,"P0"]!=0),]
   }
   heatmap.2(as.matrix(log2(toDo)),Colv=F,breaks=10,symbreaks=T,col=blue2yellow,symkey=T,trace="n",na.rm=T,dendrogram="n",key=F,cexRow=0.8,labRow=dat.norm[grep(gene,dat.norm[,"Annot"]),"common.BIG"])
   print("can't find anything...")
   return(cbind(dat.scaled.norm[grep(gene,dat.scaled.norm[,"Annot"]),"common.BIG"],dat.scaled.norm[grep(gene,dat.scaled.norm[,"Annot"]),"Annot"]))
  }
  else
  {
   return("Nothing, really nothing...")
  }
 }
###
 par(mfrow=c(3,2))
###plots scaled RNA-seq scores vs Cell Area###
 plot(log2(AREA[,1]),log2(dat.scaled[which(dat.scaled[,"common.BIG"]==gene),rna]),xlim=c(4,8),xlab="area",ylab="scaled RPKM",cex=0.7,col="red")
 legend(x="topleft",legend=paste("slope=",signif(slope[which(dat.scaled[,"common.BIG"]==gene)],2),sep=""),bty="n",cex=1.5)
###plots unscaled RNA-seq scores###
 plot(0:11,log2(dat.raw[which(dat.raw[,"common.BIG"]==gene),rna]),xlim=c(0,11),ylim=c(-2,15),xlab="time",ylab="unscaled RPKM",cex=0.7,col="blue",type="b")
 abline(h=c(0))
 abline(h=log2(c(quantile(dat.raw[,rna[1]],prob=seq(0,1,0.1),na.rm=T))),col="green",lwd=0.5)
###plots unscaled RNA-seq scores normalised to time point 0###
 plot(0:11,log2(dat.norm[which(dat.norm[,"common.BIG"]==gene),rna]),xlim=c(0,11),ylim=c(-4,10),xlab="time",ylab="log unscaled mRNA fold increase",cex=0.7,col="blue",type="b")
 abline(h=0)
 abline(h=c(-1,1),lty=2)
 abline(h=c(-2.32,2.32),lty=2,col="red")
###plots unscaled proteomics scores normalised to time point 0###
 plot(0:11,log2(dat.norm[which(dat.norm[,"common.BIG"]==gene),prot]),xlim=c(0,11),ylim=c(-4,10),xlab="time",ylab="log unscaled protein fold increase",cex=0.7,col="purple",type="b")
 abline(h=0)
 abline(h=c(-1,1),lty=2)
 abline(h=c(-2.32,2.32),lty=2,col="red")
 if(all(is.na(dat.norm[which(dat.norm[,"common.BIG"]==gene),17:28]))==T)
 legend(x="topleft",legend="NO DATA",bty="n",cex=2)
###plots scaled RNA-seq scores normalised to time point 0 alongside increase in cell area###
 plot(0:11,log2(dat.scaled.norm[which(dat.scaled.norm[,"common.BIG"]==gene),rna]),xlim=c(0,11),ylim=c(-4,10),xlab="time",ylab="log scaled fold increase",cex=0.7,col="blue",type="b")
 lines(0:11,log2(AREA[,1]/AREA[1,1]),xlim=c(0,11),ylim=c(-4,10),cex=0.7,col="red",type="b",lty=2)
 legend(x="topleft",legend=c("area",gene),text.col=c("red","blue"),bty="n",cex=2)
 abline(h=0)
 abline(h=c(-1,1),lty=2)
#abline(h=c(-2.32,2.32),lty=2,col="red")
###plots log ratio between increase in area and increase in RNA-seq scores###
 plot(0:11,(-log2(AREA[,1]/AREA[1,1])+log2(dat.scaled.norm[which(dat.scaled.norm[,"common.BIG"]==gene),rna])),ylim=c(-8,8),xlab="Time",ylab="log2 deviation from area")
 D=cbind(0:11,t(-log2(AREA[,1]/AREA[1,1])+log2(dat.scaled.norm[which(dat.scaled.norm[,"common.BIG"]==gene),rna])))
 lines(D[which(abs(D[,2])>=1),1],D[which(abs(D[,2])>=1),2],col="red",pch=20,cex=1.5,type="p")
 abline(h=0)
 abline(h=c(-1,1),lty=2,col="red")
 print(dat.scaled.norm[which(dat.scaled.norm[,"common.BIG"]==gene),"Annot"])
 print(row.names(dat.scaled.norm)[which(dat.scaled.norm[,"common.BIG"]==gene)])
}

filter.cdc=function(dat=CDC2N,cu=lhei, lwid=lwidt,n)
{
 out=0
 for(i in 1:nrow(dat))
 {
  if(cut > 1)
  {
   if(length(which(CDC2G[i,c(1:12,17:28)] > cut))>n)
   {
    out=c(out,i)
   }
  }
  else
  {
   if(length(which(CDC2G[i,c(1:12,17:28)] < cut))>n)
   {
   out=c(out,i)
   }
  }
 }
 out=out[2:length(out)]
 return(out)
}

#dev=log2(AREA[,1]/AREA[1,1])-log2(nCDC2N[,1:12])
#dev.pos=vector()
#for (i in 1:nrow(CDC2))
#{
# dev.pos[i]=length(which(abs(test[i,]) > 2))
#}



