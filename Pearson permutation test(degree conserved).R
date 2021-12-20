GSE<-as.matrix(read.table("GSE116250_rpkm.txt",header=T))##expression profile data
index<-grep("ICM",colnames(GSE))
GSE=GSE[,-index]
all_gene_id<-as.character(rownames(GSE))
trans=function(x){
  num=as.numeric(as.character(x))
  return(num)
}
GSE=apply(GSE,2,trans)
rownames(GSE)=all_gene_id
GSE_tumor<-GSE[,15:51]#disease samples
GSE_normal<-GSE[,1:14]#normal samples 
###Read the degree data of all genes in the metabolic network, and set the first column to the row name.  #####
network_degree<-as.matrix(read.table("recon_net_symbol_degree.txt",sep = ",",header = T))
network_degree<-cbind(network_degree[,1],network_degree[,2])
rownames(network_degree)<-network_degree[,1]

####Design function: Calculate the Pearson average value of the module according to the expression values of all gene pairs.####
my_function<-function(module,GSE_matrix)
{   
  
  loca=which(rownames(GSE_matrix) %in% module)
  GSE_matrix2=t(GSE_matrix[loca,])
  corps<-c();
  for(i in 1:(dim(GSE_matrix2)[1]-1))
  {
    for(j in (i+1):(dim(GSE_matrix2)[1]))
    {
      x<-as.numeric(GSE_matrix2[i,])
      y<-as.numeric(GSE_matrix2[j,])
      corp<-cor.test(x,y,alternative = "two.sided",method = "pearson")$estimate;
      corps<-c(corps,corp);
    }
    module_corp<-mean(corps);
    return(module_corp);
  }
}

####Calculate the p-value of the permutation tests for each initial module.####

times<-c()
turn<-1000 #Enter a random number of times.#
c<-c()
P<-c()
modules=as.matrix(read.table("DCM-all-module.txt",fill=T,header=F)) ##Enter data for all initial modules of the DCM/ICM/D_I.
#Calculate the Pearson difference score between samples of different states for each module.
for(j in 1:nrow(modules)){
  print(paste("Run to the ",j,"module"))
  module<-unique(as.character(unlist(strsplit(as.vector(modules[j,]),"\t"))))
  #Randomly generate 1000 degree-conserved random modules for the module.
  rand_module1000<-c()
  temp_module2<-intersect(module,rownames(network_degree));
  module_degree<-as.matrix(t(subset(t(network_degree),select = temp_module2))) 
  degree_temp<-unique(module_degree[,2])
  for(i in 1:turn){
    module_change<-c()
    for(j in 1:length(degree_temp)){
      m_rand<-as.matrix(subset(as.data.frame(module_degree), V2==degree_temp[j], select = c(V1,V2)))
      rand<-as.matrix(subset(as.data.frame(network_degree), V2==degree_temp[j],  select = c(V1,V2)))
      if(length(rand[,1])>0){
        change<-sample(rand[,1],length(m_rand[,1]),replace=T)
      }
      module_change<-c(module_change,change)
    }
    rand_module1000<-rbind(rand_module1000,module_change)
  }
  
  times<-dim(rand_module1000)[1]
  module_variablity<-abs(my_function(module,GSE_tumor)-my_function(module,GSE_normal))
  time<-0
  for(k in 1:times){
    print(paste("Run to the",k,"random module",sep=""))
    random_module<-as.matrix(rand_module1000[k,])
    random_module_variablity<-my_function(random_module,GSE_tumor)-my_function(random_module,GSE_normal);
    if(abs(module_variablity) >abs(random_module_variablity))
    {time<-time+1
    }
    else{
      time<-time
    }
  }
  P_value<-1-time/1000
  
  
  c[j]<-j
  P[j]<-P_value
}
M<-data.frame(c,P)
write.table(M,file="DCM_all_modules_pearson_Ptest_pvalue(degree conserved).txt",quote=F,sep="\t",col.names=F,row.names=F)
