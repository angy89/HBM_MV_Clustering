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

source("random_graph_functions.R")
source("vonNewmanEntropy.R")
source("local_reaching_centrality.R")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
par(mar=c(0,0,0,0))
# I perform 7 different types of simulations 
# HSim is a hierarchical symmetric matrix
# NHSim is a matrix with no hierarchical structure
# TCGABreast is the analsis performed on the TCGA breast cancer database
# nestedHierarchical 
# BTree is a binary tree structure
# Karate is the karate networks
# yeast 
Types = c("HSim","NHSim","TCGABreast","nestedHierarchical","BTree","karate","yeast")
#Type = Types[5]
nOBJ = 100

EntropyMat = matrix(0,nrow=5,ncol=2)
rownames(EntropyMat) = Types[1:5]
EntropyMat[,1]=Types[1:5]

for(Type in Types[1:5]){
  


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
if(Type=="karate"){
  data(karate)
  
  pdf(file=paste("immagini/karate/","original_graph.pdf",sep=""))
  plot(karate)
  dev.off()
  
  ADJ = get.adjacency(karate)
  ADJ = as.matrix(ADJ)
  plot_gray_image(ADJ)
  
  ADJ[ADJ>0] = 1
  pdf(file=paste("immagini/karate/","net_adj.pdf",sep=""))
  plot_gray_image(ADJ)
  dev.off()
  
  nOBJ = nrow(ADJ)
  nGenes = 1000
  nMirna = 500
  
  # The karate ADJ matrix will be used to as covariance matrix of two different multivariate normal distributions
  Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 2000)
  # This wil be a data matrix wit nGenes features and nOBJ samples where the covariance of the sample will follow
  # the structure of the karate network adjacency matrix
  GE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma1$mat))
  # This wil be a data matrix wit nmiRNA features and nOBJ samples where the covariance of the sample will follow
  # the structure of the karate network adjacency matrix
  miE <- rmvnorm(nMirna,rep(0,nOBJ),as.matrix(Sigma1$mat))
  
  colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
  colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
  d1 = abs(cor(GE)) #Correlation matrix between the samples in the two dataset
  d2 = abs(cor(miE))
  
  pdf(file=paste("immagini/karate/","GE_similarity.pdf",sep=""))
  plot_gray_image(d1) 
  dev.off()
  
  pdf(file=paste("immagini/karate/","mir_similarity.pdf",sep=""))
  plot_gray_image(d2)
  dev.off()
  
  plot_gray_image(d1)
  plot_gray_image(d2)
  
  listMatrixSimilarityTag = list(d1,d2)
  listDissimilarityMatrix = list(1-d1,1-d2)
  
  cluster.args <-c()
  cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
  
  print("HBM construction...\n")
  result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
  plot_gray_image(result_chromo_hbm[[2]]$hbm)
  plot_gray_image(result_chromo_hbm[[1]]$hbm)
  
  
  print("HBM3C algorithm...\n")
  par(mfrow=c(2,2))
  result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)

  plot_gray_image(result_chromo_hbm_3c)
  distMat = local_reaching_centrality(result_chromo_hbm_3c)
  hist(result_chromo_hbm_3c)
  boxplot(distMat)
  
  result_chromo_hbm_3c = range01(result_chromo_hbm_3c)
  result_chromo_hbm_3c = 1-result_chromo_hbm_3c
 
  plot_gray_image(result_chromo_hbm_3c)
  
  pdf(file=paste("immagini/karate/","HBM3C.pdf",sep=""))
  plot_gray_image(result_chromo_hbm_3c)
  dev.off()
  plot_gray_image(result_chromo_hbm_3c)
  
  
  
  
  
  gg = graph.adjacency(adjmatrix =result_chromo_hbm_3c,mode = "undirected",weighted = TRUE)
  gg = simplify(gg) 
  
  pdf(file=paste("immagini/karate/","HBM3C_graph.pdf",sep=""))
  
  plot(gg, vertex.label = NA,vertex.size = 5,vertex.color = V(karate)$color)
  
  dev.off()
  
  plot(gg, vertex.label = NA,
       vertex.size = 5,vertex.color = V(karate)$color,
       layout=layout_with_fr,
       edge.width = E(gg)$weigth*2)
  
  VNE = vonNewmanEntropy(gg)$VNEntropy
  write.table(x = VNE,file = "immagini/karate/VNE.txt")
  EntropyMat[Type,2]=VNE
  cat("Von Newman Entropy: ",VNE,"\n")
  EIG = vonNewmanEntropy(gg)$eigenvalues
  
  pdf(file=paste("immagini/karate/","eigenvalues.pdf",sep=""))
  plot(EIG,ylab ="eigenvalues",xlab = "")
  dev.off()
  
  plot(EIG,ylab ="eigenvalue",xlab = "")
  
  resultClusterGer <- cluster_gerarchico(list(result_chromo_hbm_3c), "complete", 4, 
                                         colnames(GE), 
                                         colnames(GE) )
  
  pdf(file=paste("immagini/karate/","hclust_HBM3C.pdf",sep=""))
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=V(karate)$color,
                             colorCodes = rainbow(2))
  dev.off()
  
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=V(karate)$color,
                             colorCodes = rainbow(2))
}
if(Type=="BTree"){
  library(igraph)
  G <- graph.tree(n=nOBJ,children=2)
  G = as.undirected(G)
  plot(G,vertex.size=5,vertex.label.cex = 0.01,edge.arrow.size = 0.01)
  x = walktrap.community(graph = G)
  ADJ = get.adjacency(G)
  hls = hclust(as.dist(distances(G)),method = "ward")
  plot(hls)
  
  pdf(file=paste("immagini/binaryTree/","initial_graph.pdf",sep=""))
  plot(G,vertex.size=5,vertex.label.cex = 0.01,edge.arrow.size = 0.01)
  dev.off()
  
  
  ADJ = get.adjacency(G)
  pdf(file=paste("immagini/binaryTree/","initial_ADJ.pdf",sep=""))
  plot_gray_image(as.matrix(ADJ))
  dev.off()
  
  plot_gray_image(as.matrix(ADJ))
  
  
  G1 <- graph.tree(n=nOBJ,children=4)
  G1 = as.undirected(G1)
  plot(G1,vertex.size=5,vertex.label.cex = 0.01,edge.arrow.size = 0.01)
  
  pdf(file=paste("immagini/binaryTree/","initial_graph2.pdf",sep=""))
  plot(G1,vertex.size=5,vertex.label.cex = 0.01,edge.arrow.size = 0.01)
  dev.off()
  
  ADJ1 = get.adjacency(G1)
  pdf(file=paste("immagini/binaryTree/","initial_ADJ2.pdf",sep=""))
  plot_gray_image(as.matrix(ADJ1))
  dev.off()
  
  plot_gray_image(as.matrix(ADJ1))
  
  library(Matrix)
  print("Finding near positive definite matrix...\n")
  Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 1000)
  Sigma2 = nearPD(x = ADJ1,corr = TRUE,maxit = 1000)
  
  nGenes = 1000
  nMirna = 500
  print("Generating GE data...\n")
  GE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma1$mat))
  print("Generating miE data...\n")
  miE <- rmvnorm(nMirna,rep(0,nOBJ),as.matrix(Sigma2$mat))
  colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
  colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
  d1 = abs(cor(GE))
  d2 = abs(cor(miE))
  diag(d1) = 0
  diag(d2) = 0
  
  pdf(file=paste("immagini/binaryTree/","GEsim.pdf",sep=""))
  plot_gray_image(d1)
  dev.off()
  
  pdf(file=paste("immagini/binaryTree/","MIRsim.pdf",sep=""))
  plot_gray_image(d2)
  dev.off()
  
  plot_gray_image(d1)
  plot_gray_image(d2)
  
  listMatrixSimilarityTag = list(d1,d2)
  listDissimilarityMatrix = list(1-d1,1-d2)
  
  cluster.args <-c()
  cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
  
  print("HBM construction...\n")
  result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
  print("HBM3C algorithm...\n")
  result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
  
  pdf(file=paste("immagini/binaryTree/","HBM3C.pdf",sep=""))
  plot_gray_image(result_chromo_hbm_3c)
  dev.off()
  
  plot_gray_image(result_chromo_hbm_3c,main ="miE Correlation Matrix")
  
  diag(result_chromo_hbm_3c) = 0
  gg = graph.adjacency(adjmatrix = max(result_chromo_hbm_3c)-result_chromo_hbm_3c,mode = "undirected",weighted = TRUE)
  gg = simplify(gg)
  y = walktrap.community(gg)
  co <- layout.reingold.tilford(gg, params=list(root=1)) 
  plot(gg, layout=co)
  
  VNE = vonNewmanEntropy(gg)$VNEntropy
  write.table(x = VNE,file = "immagini/binaryTree/VNE.txt")
  EntropyMat[Type,2]=VNE
  
  cat("Von Newman Entropy: ",VNE,"\n")
  EIG = vonNewmanEntropy(gg)$eigenvalues
  
  pdf(file=paste("immagini/binaryTree/","eigenvalues.pdf",sep=""))
  plot(EIG,ylab ="eigenvalues",xlab = "")
  dev.off()
  
  plot(EIG,ylab ="eigenvalue",xlab = "")
  
  print("Hierarchical clustering...\n")
  resultClusterGer <- cluster_gerarchico(list(result_chromo_hbm_3c), "complete", 4, 
                                         colnames(GE), 
                                         colnames(GE) )
  
  pdf(file=paste("immagini/binaryTree/","hclust_HBM3C.pdf",sep=""))
  plot(resultClusterGer[[1]]$cluster_gerarchico)
                            
  dev.off()
  
  plot(resultClusterGer[[1]]$cluster_gerarchico)
  
  
  HBM_dist = cophenetic(resultClusterGer[[1]]$cluster_gerarchico)
  
  save.image(paste(getwd(),"/","binaryTree","_back.RData",sep=""))
  
    

}
if(Type == "nestedHierarchical"){
  library(igraph)
  nGenes = 1000
  nMirna = 500
  
  nIter <- 2
  nChild <- 4
  
  # g = graph.adjacency(1,mode = "undirected")
  # root = simplify(g)
  # g=make_initial_graph()
  # g <- attach_to_center(1,root, g)
  # plot(g)
  
  g=make_initial_graph()
  
  vertex_leg = rep(1,vcount(g))
  index = 2
  for (j in 1:nIter) {
    g0 <- g
    for (i in 1:nChild){
      g <- attach_to_center(1,g, g0)
      #g = attach_to_center(1,g,make_initial_graph())
      plot(g, vertex.label = NA,vertex.size = 5,vertex.color = vertex_leg)
      
      vertex_leg = c(vertex_leg,rep(index,vcount(g0)))
      index = index +1
    }
    
  }
  plot(g, vertex.label = NA,vertex.size = 5,vertex.color = vertex_leg)
  
  
  
  pdf(file=paste("immagini/stars_net/","original_graph.pdf",sep=""))
  plot(g, vertex.label = NA,vertex.size = 5,vertex.color=vertex_leg)
  
  dev.off()
  
  plot(g, vertex.label = NA,vertex.size = 5)
  ADJ = get.adjacency(g,attr = "weight")
  ADJ = as.matrix(ADJ)
  
  ADJ[ADJ>0] = 1

  pdf(file=paste("immagini/stars_net/","original_ADJ.pdf",sep=""))
  plot_gray_image(ADJ)
  dev.off()  
  
  plot_gray_image(ADJ)
  
  nOBJ = nrow(ADJ)
  
  Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 2000)
  plot_gray_image(as.matrix(Sigma1$mat))
  
  
  GE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma1$mat))
  miE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma1$mat))
  
  
  colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
  colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
  d1 = abs(cor(GE))
  d2 = abs(cor(miE))
  
  
  pdf(file=paste("immagini/stars_net/","GE_ADJ.pdf",sep=""))
  plot_gray_image(d1)
  dev.off()  
  
  
  pdf(file=paste("immagini/stars_net/","miRNA_ADJ.pdf",sep=""))
  plot_gray_image(d2)
  dev.off()  
  
  plot_gray_image(d1)
  plot_gray_image(d2)
  
  listMatrixSimilarityTag = list(d1,d2)
  listDissimilarityMatrix = list(1-d1,1-d2)
  
  cluster.args <-c()
  cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
  
  print("HBM construction...\n")
  result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
  print("HBM3C algorithm...\n")
  result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
  
  
  pdf(file=paste("immagini/stars_net/","HBM3C_ADJ.pdf",sep=""))
  plot_gray_image(result_chromo_hbm_3c)
  dev.off()  
  
  plot_gray_image(result_chromo_hbm_3c,main ="miE Correlation Matrix")
  
  diag(result_chromo_hbm_3c) = 0
  gg = graph.adjacency(adjmatrix = max(result_chromo_hbm_3c)-result_chromo_hbm_3c,mode = "undirected",weighted = TRUE)
  plot(gg, vertex.label = NA,vertex.size = 5)
  
  pdf(file=paste("immagini/stars_net/","HBM3C_graph.pdf",sep=""))
  plot(gg, vertex.label = NA,vertex.size = 5)
  dev.off()  
  
  VNE = vonNewmanEntropy(gg)$VNEntropy
  write.table(x = VNE,file = "immagini/stars_net/VNE.txt")
  EntropyMat[Type,2]=VNE
  
  cat("Von Newman Entropy: ",VNE,"\n")
  EIG = vonNewmanEntropy(gg)$eigenvalues
  
  pdf(file=paste("immagini/stars_net/","eigenvalues.pdf",sep=""))
  plot(EIG,ylab ="eigenvalues",xlab = "")
  dev.off()
  
  plot(EIG,ylab ="eigenvalue",xlab = "")
  
  resultClusterGer <- cluster_gerarchico(list(result_chromo_hbm_3c), "complete", 4, 
                                         colnames(GE), 
                                         colnames(GE) )
  xx = cutree(resultClusterGer[[1]]$cluster_gerarchico,9)
  pdf(file=paste("immagini/stars_net/","hclust_HBM3C.pdf",sep=""))
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=vertex_leg,
                             colorCodes = rainbow(length(unique(vertex_leg))))
  dev.off()
  
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=vertex_leg,
                             colorCodes = rainbow(length(vertex_leg)))
  
  HBM_dist = cophenetic(resultClusterGer[[1]]$cluster_gerarchico)
  
  # cat("correlation between GE block matrix and HBM hierarchical clustering :", cor(as.dist(1-res1$blockMatrix),HBM_dist,method = "pearson"),"\n")
  # cat("correlation between miE block matrix and HBM hierarchical clustering :",cor(as.dist(1-res2$blockMatrix),HBM_dist,method = "pearson"),"\n")
  # 
  #save.image(paste(getwd(),"/",Type,"_back.RData",sep=""))
  
}
if(Type=="HSim"){
  blockParams = list(list(nBlocks = 2,nOBJ=30,minCorr = 0.8,maxCorr = 0.95,minNoise=0.6,maxNoise = 0.7),
                     list(nBlocks = 3,nOBJ=70,minCorr = 0.7,maxCorr = 0.85,minNoise=0.5,maxNoise = 0.6))
  res1 = nestedBlockMatrix(nOutBlocks=2,blockParams,minOutNoise = 0.1,maxOutNoise = 0.2)
  
  blockParams = list(list(nBlocks = 2,nOBJ=70,minCorr = 0.8,maxCorr = 0.95,minNoise=0.4,maxNoise = 0.7),
                     list(nBlocks = 2,nOBJ=30,minCorr = 0.9,maxCorr = 0.95,minNoise=0.4,maxNoise = 0.8))
  res2 = nestedBlockMatrix(nOutBlocks=2,blockParams,minOutNoise = 0.1,maxOutNoise = 0.2)
  
  pdf(file=paste("immagini/block_mat/","blockMat1.pdf",sep=""))
  plot_gray_image(res1$blockMatrix)
  dev.off()
  
  ADJ = res1$blockMatrix
  diag(ADJ) = 0
  ADJ[ADJ<0.18] = 0
  g = graph.adjacency(ADJ,mode="undirected",weighted = T)
  plot(g, vertex.label = NA,vertex.size = 5)
  
  pdf(file=paste("immagini/block_mat/","blockMat2.pdf",sep=""))
  plot_gray_image(res2$blockMatrix)
  dev.off()
  
  
  library(Matrix) 
  print("Finding near positive definite matrix...\n")
  Sigma1 = nearPD(x = res1$blockMatrix,corr = TRUE,maxit = 1000)
  Sigma2 = nearPD(x = res2$blockMatrix,corr = TRUE,maxit = 1000)
  
  nGenes = 1000
  nMirna = 500
  print("Generating GE data...\n")
  GE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma1$mat))
  print("Generating miE data...\n")
  miE <- rmvnorm(nMirna,rep(0,nOBJ),as.matrix(Sigma2$mat))
  colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
  colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
  d1 = abs(cor(GE))
  d2 = abs(cor(miE))
  
  pdf(file=paste("immagini/block_mat/","GEsim.pdf",sep=""))
  plot_gray_image(d1)
  dev.off()
  
  pdf(file=paste("immagini/block_mat/","MIRsim.pdf",sep=""))
  plot_gray_image(d2)
  dev.off()
  
  plot_gray_image(d1)
  plot_gray_image(d2)
  
  listMatrixSimilarityTag = list(d1,d2)
  listDissimilarityMatrix = list(1-d1,1-d2)
  
  cluster.args <-c()
  cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
  
  print("HBM construction...\n")
  result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
  print("HBM3C algorithm...\n")
  result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
  
  pdf(file=paste("immagini/block_mat/","HBM3C.pdf",sep=""))
  plot_gray_image(result_chromo_hbm_3c)
  dev.off()
  
  plot_gray_image(result_chromo_hbm_3c,main ="miE Correlation Matrix")
  
  hbm3c_adj = (max(result_chromo_hbm_3c)-result_chromo_hbm_3c)
  hbm3c_adj = hbm3c_adj/max(hbm3c_adj)
  diag(hbm3c_adj) = 0

  plot_gray_image(hbm3c_adj)
  
  gg = graph.adjacency(adjmatrix = hbm3c_adj,mode = "undirected",weighted = TRUE)
  plot(gg,vertex.label.cex = 0.000001, layout = layout_with_fr,vertex.size=5)
  
  VNE = vonNewmanEntropy(gg)$VNEntropy
  write.table(x = VNE,file = "immagini/block_mat/VNE.txt")
  EntropyMat[Type,2]=VNE
  
  cat("Von Newman Entropy: ",VNE,"\n")
  EIG = vonNewmanEntropy(gg)$eigenvalues
  
  pdf(file=paste("immagini/block_mat/","eigenvalues.pdf",sep=""))
  plot(EIG,ylab ="eigenvalues",xlab = "")
  dev.off()
  
  plot(EIG,ylab ="eigenvalue",xlab = "")
  
  print("Hierarchical clustering...\n")
  resultClusterGer <- cluster_gerarchico(list(result_chromo_hbm_3c), "complete", 4, 
                                         colnames(GE), 
                                         colnames(GE) )
  
  pdf(file=paste("immagini/block_mat/","hclust_HBM3C.pdf",sep=""))
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=res1$groupCodes,
                             colorCodes = res1$colorCodes)  
  dev.off()
  
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=res1$groupCodes,
                             colorCodes = res1$colorCodes)
  
  HBM_dist = cophenetic(resultClusterGer[[1]]$cluster_gerarchico)
  HBM_dist2 = cophenetic(hclust(as.dist(1-res1$blockMatrix),method = "complete"))
  
  cat("correlation between GE block matrix and HBM hierarchical clustering :", cor(as.dist(1-res1$blockMatrix),HBM_dist,method = "pearson"),"\n")
  cat("correlation between miE block matrix and HBM hierarchical clustering :",cor(as.dist(1-res2$blockMatrix),HBM_dist,method = "pearson"),"\n")
 # save.image(paste(getwd(),"/","HSIM","_back.RData",sep=""))
  
}
if(Type=="NHSim"){
  M1 = myCorrMat(nOBJ,0.6,0.9)
  M2 = myCorrMat(nOBJ,0.7,0.8)
  
  ADJ = M1
  diag(ADJ) = 0
  ADJ[ADJ<0.3] = 0
  g = graph.adjacency(ADJ,mode="undirected",weighted = T)
  plot(g, vertex.label = NA,vertex.size = 5)

  pdf(file=paste("immagini/random_mat/","RandMat1.pdf",sep=""))
  plot_gray_image(M1)
  dev.off()
  pdf(file=paste("immagini/random_mat/","RandMat2.pdf",sep=""))
  plot_gray_image(M2)
  dev.off()
  
  plot_gray_image(M1,main = "GE Correlation Matrix")
  plot_gray_image(M2,main ="miE Correlation Matrix")
  nGenes = 1000
  nMirna = 500
  print("Generating GE data...\n")
  GE <- rmvnorm(nGenes,rep(0,nOBJ),M1)
  print("Generating miE data...\n")
  miE <- rmvnorm(nMirna,rep(0,nOBJ),M2)
  colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
  colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
  d1 = cor(GE)
  d2 = cor(miE)
  
  pdf(file=paste("immagini/random_mat/","GEsim.pdf",sep=""))
  plot_gray_image(d1)
  dev.off()
  
  pdf(file=paste("immagini/random_mat/","MIRsim.pdf",sep=""))
  plot_gray_image(d2)
  dev.off()
  
  plot_gray_image(d1,main = "GE Correlation Matrix")
  plot_gray_image(d2,main ="miE Correlation Matrix")
  
  listMatrixSimilarityTag = list(d1,d2)
  listDissimilarityMatrix = list(1-d1,1-d2)
  
  cluster.args <-c()
  cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
  
  print("HBM construction...\n")
  result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
  print("HBM3C algorithm...\n")
  result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
  plot_gray_image(result_chromo_hbm_3c)
  
  pdf(file=paste("immagini/random_mat/","BHM3C.pdf",sep=""))
  plot_gray_image(result_chromo_hbm_3c)
  dev.off()
  
  diag(result_chromo_hbm_3c) = 0
  gg = graph.adjacency(adjmatrix = max(result_chromo_hbm_3c)-result_chromo_hbm_3c,mode = "undirected",weighted = TRUE)
  
  VNE = vonNewmanEntropy(gg)$VNEntropy
  write.table(x = VNE,file = "immagini/random_mat/VNE.txt")
  EntropyMat[Type,2]=VNE
  
  cat("Von Newman Entropy: ",VNE,"\n")
  EIG = vonNewmanEntropy(gg)$eigenvalues
  
  pdf(file=paste("immagini/random_mat/","eigenvalues.pdf",sep=""))
  plot(EIG,ylab ="eigenvalues",xlab = "")
  dev.off()
  
  plot(EIG,ylab ="eigenvalue",xlab = "")
  
  
  print("Hierarchical clustering...\n")
  resultClusterGer <- cluster_gerarchico(list(result_chromo_hbm_3c), "complete", 3, 
                                         colnames(GE), 
                                         colnames(GE) )
  
  pdf(file=paste("immagini/random_mat/","hclust_HBM3C.pdf",sep=""))
  plot(resultClusterGer[[1]]$cluster_gerarchico)
  dev.off()
  
  
  plot(resultClusterGer[[1]]$cluster_gerarchico)
  HBM_dist = cophenetic(resultClusterGer[[1]]$cluster_gerarchico)
  
  cat("correlation between GE block matrix and HBM hierarchical clustering :", cor(as.dist(1-d1),HBM_dist,method = "pearson"),"\n")
  cat("correlation between miE block matrix and HBM hierarchical clustering :",cor(as.dist(1-d2),HBM_dist,method = "pearson"),"\n")
  
}
if(Type=="TCGABreast"){
  load("TCGABRCA.RData")
  nOBJ = ncol(GE)
  print("Evaluating GE correlations\n")
  d1 = cor(GE)
  print("Evaluating mER correlations \n")
  d2 = cor(miE)
  
  
  pdf(file=paste("immagini/TCGA/","GE_sim.pdf",sep=""))
  plot_gray_image(d1,main = "GE Correlation Matrix")
  dev.off()
  
  pdf(file=paste("immagini/TCGA/","mir_sim.pdf",sep=""))
  plot_gray_image(d2,main = "GE Correlation Matrix")
  dev.off()
  
  plot_gray_image(d1,main = "GE Correlation Matrix")
  
  plot_gray_image(d2,main ="miE Correlation Matrix")
  
  listMatrixSimilarityTag = list(d1,d2)
  listDissimilarityMatrix = list(1-d1,1-d2)
  library(igraph)
  cluster.args <-c()
  cluster.method <- function(graph) cluster_louvain(graph, weights = abs(E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
  
  print("HBM construction...\n")
  result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
  print("HBM3C algorithm...\n")
  result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
  
  pdf(file=paste("immagini/TCGA/","HBM3C.pdf",sep=""))
  plot_gray_image(result_chromo_hbm_3c)
  dev.off()
  
  
  plot_gray_image(result_chromo_hbm_3c)
  
  diag(result_chromo_hbm_3c) = 0
  gg = graph.adjacency(adjmatrix = max(result_chromo_hbm_3c)-result_chromo_hbm_3c,mode = "undirected",weighted = TRUE)
  gg = simplify(gg)

  VNE = vonNewmanEntropy(gg)$VNEntropy
  write.table(x = VNE,file = "immagini/TCGA/VNE.txt")
  EntropyMat[Type,2]=VNE
  
  cat("Von Newman Entropy: ",VNE,"\n")
  EIG = vonNewmanEntropy(gg)$eigenvalues
  
  pdf(file=paste("immagini/TCGA/","eigenvalues.pdf",sep=""))
  plot(EIG,ylab ="eigenvalues",xlab = "")
  dev.off()
  
  plot(EIG,ylab ="eigenvalue",xlab = "")
  
  pdf(file=paste("immagini/TCGA/","HBM3C_graph.pdf",sep=""))
  plot(gg, vertex.label = NA,vertex.size = 5)
  dev.off()

  
  print("Hierarchical clustering...\n")
  resultClusterGer <- cluster_gerarchico(list(result_chromo_hbm_3c), "average", 5, 
                                         colnames(GE), 
                                         colnames(GE) )
  
  groupCodes = PC$x
  colorCodes <- c(Basal="red", Her2="green",LumA = "purple",LumB = "yellow")
  
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=groupCodes,
                             colorCodes = colorCodes)
  
  pdf(file=paste("immagini/TCGA/","hclust_HBM3C.pdf",sep=""))
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=groupCodes,
                             colorCodes = colorCodes)  
  dev.off()
  
  HBM_dist = cophenetic(resultClusterGer[[1]]$cluster_gerarchico)
  
  cat("correlation between GE block matrix and HBM hierarchical clustering :", cor(as.dist(1-d1),HBM_dist,method = "pearson"),"\n")
  cat("correlation between miE block matrix and HBM hierarchical clustering :",cor(as.dist(1-d2),HBM_dist,method = "pearson"),"\n")
}
}

colnames(EntropyMat) = c("SimulationType","VonNewmanEntropy")
EntropyMat = as.data.frame(EntropyMat)
EntropyMat = EntropyMat[order(EntropyMat$VonNewmanEntropy),]
