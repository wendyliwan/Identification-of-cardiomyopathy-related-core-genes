library(parmigene)
MRFms<-read.table("MRFms_ICM-module1.txt",sep = "\t",header = T)#MRFms of the candidate module  
net_degree<-as.matrix(read.table("recon_net_symbol_degree.txt",sep = "\t",header = T))##A file of all genes and their degrees in the metabolic network
net_degree<-cbind(net_degree[,2],net_degree[,1])
rownames(net_degree)<-net_degree[,1]
all_hub=as.matrix(read.table("DEGs_ICM.txt"),header = F)#ICM_NF DEGs
net=as.matrix(read.table("ICM-module1.txt"))##candidate module
mod=unique(as.character(net))
hub_mod=intersect(all_hub,mod)#DEGs in the candidate module.
nohub_mod=setdiff(mod,hub_mod)#non-DEGs in the candidate module.
all_gene=setdiff(net_degree[,1],mod)
all_hub1=intersect(all_hub,all_gene)
all_hub1_degree<-net_degree[all_hub1,]
all_nohub=setdiff(all_gene,all_hub)
all_nohub_degree<-net_degree[all_nohub,]

expre1=read.table("GSE116250_rpkm.txt",header=T)##expression profile data
name=rownames(expre1)
index<-grep("DCM",colnames(expre1))
expre1=expre1[,-index]

#####Fold Change#####
FC=function(x){
  case=x[15:27]  #disease samples
  normal=x[1:14]  #normal samples
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
####Randomly generate 1000 size-conserved random modules.####
turn<-1000
rand_module1000<-c()
rand_modhub1000<-c()
rand_modnohub1000<-c()
for(i in 1:turn){
  new_modhub<-c()
  new_modnohub<-c()
  new_modhub<-sample(all_hub1,length(hub_mod));
  new_modnohub<-sample(all_nohub,length(mod)-length(hub_mod));
  rand_modhub1000<-rbind(rand_modhub1000,new_modhub)
  rand_modnohub1000<-rbind(rand_modnohub1000,new_modnohub)
  rand_module1000<-cbind(rand_modhub1000,rand_modnohub1000)
}
####Compute the p-value of the permutation test of the candidate module####
times=dim(rand_module1000)[1]
time<-0
all_Result=c()
for(i in 1:times){
  expre=expre1
  print(paste("Run to the",i,"random module",sep=""))
  random_module<-as.matrix(rand_module1000[i,])
  hub=as.matrix(rand_modhub1000[i,])
  mod=unique(as.character(random_module))
  module=t(random_module)
  loca=which(rownames(expre1) %in% module)
  GSE_matrix1=expre1[loca,]
  mean=t(as.matrix(apply(GSE_matrix1,2,mean)))
  index=which(!(module %in% rownames(expre1)))
  NAgene=module[,index]
  GSE_matrix2=c()
  if(length(NAgene) > 0){
    for(n in 1:ncol(module)){
      if(!(module[,n] %in% rownames(expre1))==TRUE){
        GSE_matrix2=rbind(GSE_matrix2,mean)
      }
    }
    rownames(GSE_matrix2)=NAgene
  }
  expre=rbind(expre,GSE_matrix2)
  
  g=length(mod)#The number of genes in the random module
  #the expression differences of DEGs in the random module
  fi_hub=function1(as.character(hub),expre)
  
  result1=(sum(fi_hub))/sqrt(g)
  
  no_hub=setdiff(mod,as.character(hub))#non-DEGs in the random module
  if(length(no_hub) > 0){
    no_hub=as.matrix(no_hub)
    nohub_expre=merge(no_hub,expre,by.x="V1",by.y="row.names")
    rownames(nohub_expre)=nohub_expre[,1]
    nohub_expre=nohub_expre[,2:ncol(nohub_expre)]
    #the expression differences of non-DEGs in the random module
    fi_nohub=as.matrix(function1(no_hub,expre))
    ##degree of non-DEGs 
    degree_nohub=merge(no_hub,net_degree,by="V1")
    
    b=nrow(fi_nohub)
    d=nrow(degree_nohub)
    fi=as.numeric(fi_nohub[,1])
    data1=NULL
    if(b == d){
      #If the number of non-DeGs in the random module is greater than 1, Mutual Information is calculated.
      if(nrow(nohub_expre) > 1){
        mi_nohub=knnmi.all(nohub_expre,3) #Mutual Information of non-DEGs.
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
      #If the number of non-DeGs in the random module is less than 2, the mutual information is 0.
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
  else{
    result2=0
  }
  
  Result=result1-result2 #MRFms of the random module.
  all_Result=c(all_Result,Result)
  
  if(Result < MRFms[1,1]){
    time<-time+1
  }
  
  
}

P_value<-1-(time/1000)
write.table(P_value,file="ICM-module1 MRFms Ptest pvalue(size conserved).txt",quote=F,sep="\t",col.names=F,row.names=F)