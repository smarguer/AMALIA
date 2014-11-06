conv=read.delim("P:\\SAM_DATA/PomBase2UniProt_061014.TXT2",stringsAsFactors=F,header=F)
colnames(conv)=c("Systematic","Uniprot")

P=read.delim("../AMALIA/LF_031014.txt", stringsAsFactors=F,row.names=1)
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
  row.names(P)[i]=conv[which(conv[,2]==row.names(P)[i]),1]
 }
 else
 {
  print(row.names(P)[i])
 }
}





