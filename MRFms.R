library(parmigene)
expre=read.table("GSE116250_rpkm.txt",header=T)##expression profile data
index<-grep("ICM",colnames(expre))
expre=expre[,-index]
net_degree<-as.matrix(read.table("recon_net_symbol_degree.txt",sep = "\t",header = T))##A file of all genes and their degrees in the metabolic network
net_degree<-cbind(net_degree[,2],net_degree[,1])
net=as.matrix(read.table("DCM-module1.txt"))##candidate module
hub=as.matrix(read.table("DCM-module1_hub.txt",sep=""))##DEGs in the candidate module.
mod=unique(as.character(net))
#####Fold Change#####
FC=function(x){
  case=x[15:51]  #disease samples
  normal=x[1:14] #normal samples 
  mean_1=mean(case)
  mean_2=mean(normal)
  result=mean_1/mean_2
  return(result)
}
####Calculate the expression differences of gene.####
function1=function(x,expre){
  name=intersect(rownames(expre),as.character(x))
  location=NULL
  for(i in 1:length(name)){
    loca=which(rownames(expre) == name[i])
    location=c(location,loca)
  }
  gse_matrix<-expre[location,]
  f_hub=apply(gse_matrix,1,FC)
  return(f_hub)
}
g=length(mod)
#the expression differences of DEGs in the module
fi_hub=function1(as.character(hub),expre)
result1=(sum(fi_hub))/sqrt(g)
no_hub=setdiff(mod,as.character(hub))#non-DEGs in the module 
result2=0 
if(length(no_hub) > 0){
      no_hub=as.matrix(no_hub)
      nohub_expre=merge(no_hub,expre,by.x="V1",by.y="row.names")
      rownames(nohub_expre)=nohub_expre[,1]
      nohub_expre=nohub_expre[,2:ncol(nohub_expre)]
      #the expression differences of non-DEGs in the module
      fi_nohub=as.matrix(function1(no_hub,expre))
      ##degree of non-DEGs 
      degree_nohub=merge(no_hub,net_degree,by="V1")

      b=nrow(fi_nohub)
      d=nrow(degree_nohub)
      fi=as.numeric(fi_nohub[,1])
      data1=NULL
      if(b == d){
        #If the number of non-DeGs in the module is greater than 1, Mutual Information is calculated.
        if(nrow(nohub_expre) > 1){
          mi_nohub=knnmi.all(nohub_expre,3)#Mutual Information of non-DEGs.
          for(j in 1:b-1){
            loca=which(degree_nohub[,1] == rownames(fi_nohub)[j])
            c1=fi[j]/sqrt(as.numeric(degree_nohub[loca,2]))
            for (k in (j+1):b){
              loca1=which(degree_nohub[,1] == rownames(fi_nohub)[k])
              c2=(c1-(fi[k]/sqrt(as.numeric(degree_nohub[loca1,2]))))^2
              data1=c(data1,(c2*mi_nohub[j,k]))
            }
          }
        }
        #If the number of non-DeGs in the module is less than 2, the mutual information is 0.
        if(nrow(nohub_expre) < 2){
          mi_nohub=0
          for(j in 1:b-1){
            loca=which(degree_nohub[,1] == rownames(fi_nohub)[j])
            c1=fi[j]/sqrt(as.numeric(degree_nohub[loca,2]))
            for (k in (j+1):b){
              loca1=which(degree_nohub[,1] == rownames(fi_nohub)[k])
              c2=(c1-(fi[k]/sqrt(as.numeric(degree_nohub[loca1,2]))))^2
              data1=c(data1,(c2*mi_nohub))
            }
          }
        }
      }
      result=sum(data1)
      result2=result/nrow(net)
    }


Result=result1-result2 #MRFms of the candidate module.
write.table(Result,"MRFms_DCM-module1.txt",sep="",row.names=F,col.names=F,quote=F)
