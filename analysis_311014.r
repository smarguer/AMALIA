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
    #row.names(P)[i]=conv[which(conv[,2]==row.names(P)[i]),1]
   }
  }
  else
  {
   #print(row.names(P)[i])
  }
 }

