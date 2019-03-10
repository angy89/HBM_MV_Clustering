source("blockMatrixDiagonal.R")
source("nestedBlockMatrix.R")
source('progetto_tesi/chromoHBM.R')
source('progetto_tesi/chromoHBM3C.R')
source('progetto_tesi/tesi_createHBM.R')
source('progetto_tesi/tesi_createHBM3C.R')
source('progetto_tesi/tesi_cluster_gerarchico.R')
source('progetto_tesi/confusion_matrix.R')
source('progetto_tesi/CMsup_error.R')
library(mvtnorm)
library(entropy)
library(igraph)
library(Matrix)
library(igraphdata)
library(psych)
library(gplots)

source("random_graph_functions.R")
source("vonNewmanEntropy.R")
source("local_reaching_centrality.R")
source("utilities.R")

# I perform 7 different types of simulations 
# HSim is a hierarchical symmetric matrix
# NHSim is a matrix with no hierarchical structure
# TCGABreast is the analsis performed on the TCGA breast cancer database
# nestedHierarchical 
# BTree is a binary tree structure
# Karate is the karate networks
# yeast 
Types = c("HSim","NHSim","nestedHierarchical","BTree","karate","TCGABreast","syntren")#""
#Type = Types[5]
nOBJ = 100

EntropyMat = matrix(0,nrow=6,ncol=2)
rownames(EntropyMat) = Types[1:6]
EntropyMat[,1]=Types[1:6]
nEigenvalues = c()
HBM_rank = c()

for(Type in Types[1:6]){
  
  
  
  # 
  # if(type == "yeast"){
  #   library(igraphdata)
  #   data(yeast)
  #   
  #   ADJ = get.adjacency(yeast)
  #   ADJ = as.matrix(ADJ)
  #   ADJ[ADJ>0] = 1
  #   
  #   plot_gray_image(ADJ,main = "GE Correlation Matrix")
  #   
  # }
  
  if(Type=="syntren"){
    gen_data = gen_syntren_data(toPlot = TRUE)
    res = complete_analysis(gen_data,nFeat1=1000,nFeat2=500,toSave=FALSE,filename="",
                            toClust = FALSE,plotLocation="figures/syntren.png")
    dev.off()
    
    nEigenvalues = c(nEigenvalues, evaluate_svd(res$HBM3C_mat))
    HBM_rank = c(HBM_rank, rankMatrix(res$HBM3C_mat)[1])
    
    
    EntropyMat[Type,2]=res$net_res$VNE
  }
  if(Type=="karate"){
    gen_data = gen_karate_data()
    res = complete_analysis(gen_data,nFeat1=1000,nFeat2=500,toSave=FALSE,filename="",
                                 toClust = TRUE,groupCodes=V(karate)$color,colorCodes=rainbow(2),
                            plotLocation= "figures/karate.png")
    plot(karate)
    dev.off()
    
    nEigenvalues = c(nEigenvalues, evaluate_svd(res$HBM3C_mat))
    HBM_rank = c(HBM_rank, rankMatrix(res$HBM3C_mat)[1])
    
    
    EntropyMat[Type,2]=res$net_res$VNE
  }
  if(Type=="BTree"){
    gen_data = gen_btree_data(nOBJ = 151,nChild1 = 2,nChild2 = 4)
    
    res = complete_analysis(gen_data,nFeat1=1000,nFeat2=500,toSave=FALSE,filename="",
                      toClust = TRUE,groupCodes=NULL,colorCodes=NULL,
                      plotLocation= "figures/btree.png")
    plot(gen_data$G,vertex.size=5,vertex.label.cex = 0.01,edge.arrow.size = 0.01)
    
    dev.off()
    
    nEigenvalues = c(nEigenvalues, evaluate_svd(res$HBM3C_mat))
    HBM_rank = c(HBM_rank, rankMatrix(res$HBM3C_mat)[1])
    
    
    EntropyMat[Type,2]=res$net_res$VNE
    
  }
  if(Type == "nestedHierarchical"){
    gen_data = gen_nestedHierarchical_data(nIter = 3,n_nodes = 5,nChild =2,toPlot=FALSE)
    res = complete_analysis(gen_data,nFeat1=1000,nFeat2=500,toSave=FALSE,filename="",
                            toClust = TRUE,groupCodes=gen_data$vertex_leg,
                            colorCodes=rainbow(length(gen_data$vertex_leg)),
                            plotLocation= "figures/nestedH.png")
    plot(gen_data$g, vertex.label = NA,vertex.size = 5)
    
    dev.off()
    
    nEigenvalues = c(nEigenvalues, evaluate_svd(res$HBM3C_mat))
    HBM_rank = c(HBM_rank, rankMatrix(res$HBM3C_mat)[1])
    
    EntropyMat[Type,2]=res$net_res$VNE
    
  }
  if(Type=="HSim"){
    blockParams1 = list(list(nBlocks = 4,nOBJ=100,minCorr = 0.8,maxCorr = 0.99,minNoise=0.4,maxNoise = 0.6),
                       list(nBlocks = 1,nOBJ=51,minCorr = 0.7,maxCorr = 0.8,minNoise=0.3,maxNoise = 0.6))
    blockParams2 = list(list(nBlocks = 3,nOBJ=75,minCorr = 0.8,maxCorr = 0.9,minNoise=0.5,maxNoise = 0.6),
                       list(nBlocks = 2,nOBJ=76,minCorr = 0.65,maxCorr = 0.85,minNoise=0.2,maxNoise = 0.3))
    gen_data =  gen_Hsim_data(blockParams1 =blockParams1 ,blockParams2 = blockParams2)
    res = complete_analysis(gen_data,nFeat1=1000,nFeat2=500,toSave=FALSE,filename="",
                            toClust = TRUE,groupCodes=gen_data$res1$groupCodes,
                            colorCodes=gen_data$res1$colorCodes,
                            plotLocation= "figures/HSIM.png")
    dev.off()
    
    nEigenvalues = c(nEigenvalues, evaluate_svd(res$HBM3C_mat))
    HBM_rank = c(HBM_rank, rankMatrix(res$HBM3C_mat)[1])
    
    EntropyMat[Type,2]=res$net_res$VNE
    
  }
  if(Type=="NHSim"){
    gen_data = gen_NHsim_data(nOBJ=151,toPlot=TRUE)  
    res = complete_analysis(gen_data,nFeat1=1000,nFeat2=500,toSave=FALSE,filename="",
                            toClust = TRUE,groupCodes=NULL,
                            colorCodes=NULL,plotLocation= "figures/NHSim.png") 
    dev.off()
    
    nEigenvalues = c(nEigenvalues, evaluate_svd(res$HBM3C_mat))
    HBM_rank = c(HBM_rank, rankMatrix(res$HBM3C_mat)[1])
    
    
    EntropyMat[Type,2]=res$net_res$VNE
    
  }
  if(Type=="TCGABreast"){
    load("TCGABRCA.RData")
     d1 = abs(cor(GE,method = "spearman"))
     d2 = abs(cor(miE,method = "spearman"))
    
    HBM3C_mat = HBM3C_function(d1,d2)
    net_res = network_analysis(HBM3C_mat)
    
    heat_map =     heatmap.2(HBM3C_mat)
    gg = graph.adjacency(adjmatrix =heat_map$carpet,mode = "undirected",weighted = TRUE)
    gg = simplify(gg) 
    
    vonNewmanEntropy(gg)$VNEntropy
    
    
    png("figures/TGCA.png",width = 800,height = 600)
    
    par(mfrow=c(3,3))
    plot(1,type="n")
    plot(1,type="n")
    
    plot_gray_image(d1,main="Distance View 1")
    plot_gray_image(d2,main="Distance View 2")
    
    plot_gray_image(HBM3C_mat,main="HBM matrix")
    plot_gray_image(heat_map$carpet)
    
    #distMat = local_reaching_centrality(HBM3C_mat)
    hist(HBM3C_mat,main="HBM values distribution")
    #boxplot(distMat,main = "Local Reaching Centrality")
    plot(net_res$EIG[1:10],ylab ="eigenvalue",xlab = "",main="HBM matrix eigenvalues")
    clust_analysis(HBM3C_mat,sample_names = colnames(GE),
                   groupCodes=NA,
                   colorCodes =NA) 
    dev.off()
    nEigenvalues = c(nEigenvalues, evaluate_svd(res$HBM3C_mat))
    HBM_rank = c(HBM_rank, rankMatrix(res$HBM3C_mat)[1])
    
    EntropyMat[Type,2]=net_res$VNE
    
   }
}

colnames(EntropyMat) = c("SimulationType","VonNewmanEntropy")
EntropyMat = as.data.frame(EntropyMat)
EntropyMat = EntropyMat[order(as.numeric(as.vector(EntropyMat$VonNewmanEntropy)),decreasing = TRUE),]
