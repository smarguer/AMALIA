Zero2Na=function(dat)
{
 for(i in 1:ncol(dat))
 {
  dat[which(dat[,i]==0),i]=NA
 } 
 return(dat)
}







