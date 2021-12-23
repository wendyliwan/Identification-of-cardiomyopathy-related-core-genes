library(plyr)
library(igraph)
net<-read.table("recon_net.txt",header = T,sep = "\t")
gene<-read.table("know gene.txt",sep = "\t")
g1<- graph_from_data_frame(net, directed = F)
#plot(g1)
for (i in 1:(nrow(gene)-1)) {
  print(paste("i=",i,sep=""))
  for (j in (i+1):nrow(gene)) {
    print(paste("j=",j,sep=""))
    a=get.all.shortest.paths(g1, gene[i,1],gene[j,1])
    df <- as.data.frame(t(sapply(a$res, as_ids)))
    path_gene<-c()
    for (k in 1:ncol(df)) {
      path_gene<-c(path_gene,df[,k])
      
    }
    path_gene<-unique(path_gene)
    write.table(df,file = paste(gene[i,1],gene[j,1],"shortest path.txt"),sep="\t",row.names = F,col.names = F,quote=F)
    write.table(path_gene,file = paste(gene[i,1],gene[j,1],"shortest path gene.txt"),sep="\t",row.names = F,col.names = F,quote=F)
    
  }
  
}



