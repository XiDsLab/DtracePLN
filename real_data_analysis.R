library(Seurat)
library(PLNet)
library(ggplot2)
library(PLNmodels)
library(GMPR)
library(SeuratData)
library(enrichR)
library(reshape2)

###read data "ifnb"
data("ifnb")
count_mat<-(ifnb@assays$RNA@counts)
anno_vec<-ifnb$seurat_annotations
anno_vec_2<-ifnb$orig.ident

###choose CD14
CD14_stim<-(which(anno_vec == names(table(anno_vec))[1] & anno_vec_2 == names(table(anno_vec_2))[2]))
count_CD14_stim<-count_mat[,CD14_stim]

###choosing 200 high variable genes

P3se_ifnb_stim = CreateSeuratObject(counts = count_CD14_stim,min.cells = 3)
P3se_ifnb_stim <- NormalizeData(P3se_ifnb_stim, normalization.method = "LogNormalize", scale.factor = 10000)
P3se_ifnb_stim <- FindVariableFeatures(P3se_ifnb_stim,nfeatures = 200)
variable_gene_stim<-Seurat::VariableFeatures(P3se_ifnb_stim)

count_CD14_stim_ifnb<-count_CD14_stim[which(row.names(count_CD14_stim)%in%variable_gene_stim),]
datafinal.stim<-count_CD14_stim_ifnb[which(rowSums(ifelse(as.matrix(count_CD14_stim_ifnb)>1,1,0))>=1),]
datafinal.stim<-as.matrix(datafinal.stim)

###Do PLNet
A.stim<-PLNet(t(datafinal.stim),Sd_est =GMPR(t(datafinal.stim),1,1))
###Do VPLN
Y.stim=t(datafinal.stim)
row.names(Y.stim)<-1:nrow(Y.stim)
original_list<-list(data_1=as.data.frame(Y.stim),
                    Covariate=as.data.frame(matrix(rep(0,dim(Y.stim)[1]),ncol = 1)))
pre_data<-prepare_data(counts = original_list$data_1,covariates = original_list$Covariate,offset = "TSS")

fits <- PLNnetwork(Abundance ~ log(Offset), data = pre_data,control_init = list(nPenalties=80,min.ratio=0.000001))
bm.stim<-fits$models

######Network choosing
###Network density setting
den<-0.07
###Choosing network by density
###PLNet
density<-c()
for(i in 1:length(A.stim$lambda_vec)){
  density[i]<-(sum(ifelse(A.stim$Omega_est[[i]]!=0,1,0))-200)/200/199
}
index<-which.min(abs(density-den))
omega.stim<-as.matrix(A.stim$Omega_est[[index]])
###VPLN
density<-c()
for(i in 1:80){
  density[i]<-bm.stim[[i]]$density
}
index<-which.min(abs(density-den))
pre_V_choose.stim<-bm.stim[[index]]$model_par$Omega
omega.stim<-pre_V_choose.stim

###Choosing graph by BIC
omega.stim<-as.matrix(A.stim$Omega_chooseB)
omega.stim<-bm.stim[[which.max(fits$criteria$BIC)]]$model_par$Omega

###calculate partial correlation
pomega.stim<-as.matrix(omega.stim)
n<-dim(pomega.stim)[1]
for (i in 1:n){
  for (j in 1:n){
    pomega.stim[i,j]<-omega.stim[i,j]/sqrt(omega.stim[i,i]*omega.stim[j,j])*(-1)
  }
}

#########All the Follow-up analysis for networks are based on the "pomega.stim" or "omega.stim" 
#########chosen by density or BIC aforehand.  
######heatmap plot
###GO cluster
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018" ,"KEGG_2019_Human")
enriched.D <- enrichr(row.names(omega.stim), dbs)
bp.stim <- enriched.D[["GO_Biological_Process_2018"]]
geneclu1<-strsplit(bp.stim$Genes[3],";")[[1]]
clu1index<-which(row.names(omega.stim)%in%geneclu1)

enriched.D1 <- enrichr(row.names(omega.stim)[-clu1index], dbs)
bp.stim1 <- enriched.D1[["GO_Biological_Process_2018"]]
geneclu2<-strsplit(bp.stim1$Genes[1],";")[[1]]
clu2index<-which(row.names(omega.stim)%in%geneclu2)

enriched.D2 <- enrichr(row.names(omega.stim)[-c(clu1index,clu2index)], dbs)
bp.stim2 <- enriched.D2[["GO_Biological_Process_2018"]]
geneclu3<-strsplit(bp.stim2$Genes[7],";")[[1]]
clu3index<-which(row.names(omega.stim)%in%geneclu3)

enriched.D3 <- enrichr(row.names(omega.stim)[-c(clu1index,clu2index,clu3index)], dbs)
bp.stim3 <- enriched.D3[["GO_Biological_Process_2018"]]
geneclu4<-union(strsplit(bp.stim3$Genes[c(26)],";")[[1]],strsplit(bp.stim3$Genes[c(755)],";")[[1]])
clu4index<-which(row.names(omega.stim)%in%geneclu4)

cluindex<-c(clu1index,clu2index,clu3index,clu4index)

###calculate within-between connection ratios of 4 modules
###pomega.stim<-ifelse(pomega.stim!=0,1,0) ##add it for calculating unweighted ratios
sum1<-sum(abs(pomega.stim[clu1index,clu1index]))-length(clu1index)
sum2<-sum(abs(pomega.stim[clu2index,clu2index]))-length(clu2index)
sum3<-sum(abs(pomega.stim[clu3index,clu3index]))-length(clu3index)
sum4<-sum(abs(pomega.stim[clu4index,clu4index]))-length(clu4index)
sum1all<-sum(abs(pomega.stim[clu1index,cluindex]))-length(clu1index)
sum2all<-sum(abs(pomega.stim[clu2index,cluindex]))-length(clu2index)
sum3all<-sum(abs(pomega.stim[clu3index,cluindex]))-length(clu3index)
sum4all<-sum(abs(pomega.stim[clu4index,cluindex]))-length(clu4index)
ratio<-c(sum1,sum2,sum3,sum4)/c(sum1all,sum2all,sum3all,sum4all)

###plot
pomega.stim<-pomega.stim[cluindex,cluindex]
omega.melt<-melt(as.matrix(pomega.stim))
omega.melt$value[which(omega.melt$value==-1)]<-0
ggplot(omega.melt, aes(x=Var2,y=Var1))+geom_tile(aes(fill=value))+
  scale_fill_gradient2(low =  "blue" ,mid="white", high =  "red" )+
  #theme(axis.text.x = element_blank())+
  #theme(axis.text.y = element_text(size = 5,vjust = 0.3, hjust = 0.5))+
  geom_rect(aes(xmin = 1 - 0.5, xmax = 35 + 0.5, ymin = 0, ymax = 0),
            fill = "transparent", color = '#E64E00', size = 3)+
  geom_rect(aes(xmin = 36 - 0.5, xmax = 50 + 0.5, ymin = 0, ymax = 0),
            fill = "transparent", color = '#FF7F00', size = 3)+
  geom_rect(aes(xmin = 51 - 0.5, xmax = 61 + 0.5, ymin = 0, ymax = 0),
            fill = "transparent", color = '#65B48E', size = 3)+
  geom_rect(aes(xmin = 62 - 0.5, xmax = 73 + 0.5, ymin = 0, ymax = 0),
            fill = "transparent", color = '#3E5CC5', size = 3)+
  geom_rect(aes(ymin = 1 - 0.5, ymax = 35 + 0.5, xmin = 0, xmax = 0),
            fill = "transparent", color = '#E64E00', size = 3)+
  geom_rect(aes(ymin = 36 - 0.5, ymax = 50 + 0.5, xmin = 0, xmax = 0),
            fill = "transparent", color = '#FF7F00', size = 3)+
  geom_rect(aes(ymin = 51 - 0.5, ymax = 61 + 0.5, xmin = 0, xmax = 0),
            fill = "transparent", color = '#65B48E', size = 3)+
  geom_rect(aes(ymin = 62 - 0.5, ymax = 73 + 0.5, xmin = 0, xmax = 0),
            fill = "transparent", color = '#3E5CC5', size = 3)+
  geom_rect(aes(xmin = 1 - 0.5, xmax = 35 + 0.5, ymin = 1 - 0.5, ymax = 35 + 0.5),
            fill = "transparent", color = '#E64E00', size = 0.5)+
  geom_rect(aes(xmin = 36 - 0.5, xmax = 50 + 0.5, ymin = 36 - 0.5, ymax = 50 + 0.5),
            fill = "transparent", color = '#FF7F00', size = 0.5)+
  geom_rect(aes(xmin = 51 - 0.5, xmax = 61 + 0.5, ymin = 51 - 0.5, ymax = 61 + 0.5),
            fill = "transparent", color = '#65B48E', size = 0.5)+
  geom_rect(aes(xmin = 62 - 0.5, xmax = 73 + 0.5, ymin = 62 - 0.5, ymax = 73 + 0.5),
            fill = "transparent", color = '#3E5CC5', size = 0.5)+
  labs(y=element_blank(),x=element_blank())+
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(),axis.text.y = element_blank())

######genes connecting to IFNB1
rownames(datafinal.stim)[which(omega.stim[193,]!=0)]

######degree plot for genes which only expressed in the stim cells
###Choosing genes which only expressed in the stim cells
CD14_ctrl<-(which(anno_vec == names(table(anno_vec))[1] & anno_vec_2 == names(table(anno_vec_2))[1]))
count_CD14_ctrl<-count_mat[,CD14_ctrl]
count_CD14_ctrl_ifnb<-count_CD14_ctrl[which(row.names(count_CD14_ctrl)%in%variable_gene_stim),]
zeroindex<-which(rowSums(as.matrix(count_CD14_ctrl_ifnb))==0)

###plot
featureD<-c()
degreeD<-c()
featureV<-c()
degreeV<-c()
for (i in 1:100){
  omega<-as.matrix(A.stim$Omega_est[[i]])
  featureD[i]<-(sum(ifelse(omega!=0,1,0))-200)/2
  degreeD[i]<-sum(rowSums(ifelse(omega!=0,1,0))[zeroindex])-14
}
for (i in 1:80){
  omega<-as.matrix(bm.stim[[i]]$model_par$Omega)
  featureV[i]<-(sum(ifelse(omega!=0,1,0))-200)/2
  degreeV[i]<-sum(rowSums(ifelse(omega!=0,1,0))[zeroindex])-14
}

dt2<-data.frame(feature1=c(featureD,featureV,featureD,featureV),degree=c(degreeD,degreeV,0.14*featureD,0.14*featureV),Method=rep(c('PLNet','VPLN','Random'),times = c(100,80,180)))
ggplot(data = dt2, mapping = aes(x = feature1, y = degree, colour=Method)) + geom_line(size = 1.5) +
  xlab('Edge counts')+ylab('Total degree')+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "black"))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size = 20))+
  theme(axis.text.x = element_text(size = 20),axis.text.y = element_text(size = 20))+
  theme(legend.text = element_text(size = 20),legend.title = element_text(size = 20))+
  scale_x_continuous(limits=c(0, 2000))+
  scale_y_continuous(limits=c(0, 400))

######ChIP-seq
###Read ChIP-seq data
TF<-read.table("./real_data_analysis_ref/attribute_list_entries.txt",header=T,na.strings = c("NA"))
TF2<-read.table("./real_data_analysis_ref/attribute_list_entries2.txt",header=T,na.strings = c("NA"))
TF_target<-read.table("./real_data_analysis_ref/gene_attribute_edges.txt",header=T,na.strings = c("NA"))
TF_target2<-read.table("./real_data_analysis_ref/gene_attribute_edges2.txt",header=T,na.strings = c("NA"))

###Establish the silver standard from ChIP-seq data
TFset<-row.names(omega.stim)[which(row.names(omega.stim)%in%union(TF$GeneSym,TF2$GeneSym))]
TF_target_sub<-TF_target[which(TF_target$target%in%TFset),]
TF_target_sub<-TF_target_sub[which(TF_target_sub$source%in%row.names(omega.stim)),]
TF_target_sub2<-TF_target2[which(TF_target2$target%in%TFset),]
TF_target_sub2<-TF_target_sub2[which(TF_target_sub2$source%in%row.names(omega.stim)),]
TF_target_sub<-rbind(TF_target_sub,TF_target_sub2)
TF_target_sub<-TF_target_sub[,c(1,4)]
TF_target_sub<-unique(TF_target_sub)

###Calculting findcount and truecount
omega<-as.matrix(omega.stim)
omega[omega!=0]<-1
TF_target_sub1<-TF_target_sub[which(TF_target_sub$source%in%row.names(omega)&TF_target_sub$target%in%row.names(omega)),]
TF_target_sub1<-subset(TF_target_sub1,as.vector(source)!=as.vector(target))
truereg<-TF_target_sub1
totalcount<-dim(TF_target_sub1)[1]
index<-which(row.names(omega)%in%TF_target_sub1$target)
findcount<-sum(rowSums(omega[index,index])-1)/2+sum(rowSums(omega[index,-index]))
truecount<-0
for(j in 1:totalcount){
  if(omega[which(row.names(omega)%in%TF_target_sub1$target[j]),which(row.names(omega)%in%TF_target_sub1$source[j])]!=0){
    truecount<-truecount+1
  }
}
