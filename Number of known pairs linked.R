pathname1<-c()
myfiles0 <- Sys.glob("*shortest path gene*.txt")
gene=as.matrix(read.table("DCM_module1.txt",fill=T,header=F))
for(j in 1:nrow(gene)){
  time=0
  print(paste("Run to the ",i,"th gene",sep=""))
  gene1<-gene[j,]
  for (i in 1:length(myfiles0)){
    path<-read.table(myfiles0[i],sep="\t",header = F,stringsAsFactors = F)
    index=which(gene1 %in% path[,1])
    if(length(index)>0){
      time=time+1
    }
    
  }  
  a<-cbind(gene1,time)
  pathname1<-rbind(pathname1,a)
}
write.table(pathname1,file="module_gene_times.txt",quote=F,sep="\t",col.names=T,row.names=F)