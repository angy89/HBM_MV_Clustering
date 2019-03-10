#Utilities functions

evaluate_svd = function(HBM3C_mat){
  svd_res = svd(HBM3C_mat)
  par(mfrow=c(2,3))
  plot_gray_image(HBM3C_mat)
  for(j in 2:9){
    plot_gray_image(svd_res$u[,1:j] %*% diag(svd_res$d[1:j]) %*% t(svd_res$v[,1:j]))
  }
  return(sum(svd_res$d>1))
}


gen_karate_data = function(toPlot=FALSE){
  data(karate)
  #plot(karate)
  ADJ = get.adjacency(karate)
  ADJ = as.matrix(ADJ)
  ADJ[ADJ>0] = 1
  if(toPlot)plot_gray_image(ADJ)
  
  # The karate ADJ matrix will be used to as covariance matrix of two different multivariate normal distributions
  Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 2000)
  
  toRet = list(ADJ1 = ADJ,ADJ2=ADJ,Sigma1=Sigma1,Sigma2=Sigma1)
  return(toRet)
}

gen_syntren_data = function(toPlot=FALSE){
  require(grndata)
  
  data = syntren300.data
  data = as.matrix(data)
  gold = syntren300.net
  gold = as.matrix(gold)
  
  #plot(karate)
  ADJ = gold
 
  if(toPlot)plot_gray_image(ADJ)
  
  # The karate ADJ matrix will be used to as covariance matrix of two different multivariate normal distributions
  Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 2000)
  
  toRet = list(ADJ1 = ADJ,ADJ2=ADJ,Sigma1=Sigma1,Sigma2=Sigma1)
  return(toRet)
}


gen_Hsim_data = function(blockParams1,blockParams2,toPlot=FALSE){
  res1 = nestedBlockMatrix(nOutBlocks=2,blockParams1,minOutNoise = 0.2,maxOutNoise = 0.4)
  res2 = nestedBlockMatrix(nOutBlocks=2,blockParams2,minOutNoise = 0.2,maxOutNoise = 0.4)
  print("Finding near positive definite matrix...\n")
  Sigma1 = nearPD(x = res1$blockMatrix,corr = TRUE,maxit = 1000)
  Sigma2 = nearPD(x = res2$blockMatrix,corr = TRUE,maxit = 1000)
  if(toPlot)plot_gray_image(Sigma1$mat)
  
  toRet = list(ADJ1 = res1$blockMatrix,ADJ2=res2$blockMatrix,Sigma1=Sigma1,Sigma2=Sigma2,res1=res1,res2=res2)
  return(toRet)
}


gen_NHsim_data = function(nOBJ,toPlot=FALSE){
  M1 = myCorrMat(nOBJ,0.6,0.9)
  M2 = myCorrMat(nOBJ,0.7,0.8)
  if(toPlot)plot_gray_image(M1)
  
  toRet = list(ADJ1 = M1,ADJ2=M2,Sigma1=list(mat=M1),Sigma2=list(mat=M2))
  return(toRet)
}


gen_btree_data = function(nOBJ = 100,nChild1 = 2,nChild2 = 4,toPlot=FALSE){
  G <- graph.tree(n=nOBJ,children=nChild1)
  G = as.undirected(G)
  ADJ = get.adjacency(G)
  G1 <- graph.tree(n=nOBJ,children=nChild2)
  G1 = as.undirected(G1)
  ADJ1 = get.adjacency(G1)

  if(toPlot)plot_gray_image(ADJ1)
  
  print("Finding near positive definite matrix...\n")
  Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 1000)
  Sigma2 = nearPD(x = ADJ1,corr = TRUE,maxit = 1000)
  
  toRet = list(ADJ1 = ADJ,ADJ2=ADJ1,Sigma1=Sigma1,Sigma2=Sigma2,G=G,G1=G1)
  return(toRet)
}

gen_nestedHierarchical_data = function(nIter = 2,nChild = 2,n_nodes = 3,toPlot=TRUE){
  

  g=make_initial_graph(n_nodes = n_nodes)
  
  vertex_leg = rep(1,vcount(g))
  index = 2
  pb = txtProgressBar(min=1,max = nIter,style=3)
  for (j in 1:nIter) {
    g0 <- g
    for (i in 1:nChild){
      g <- attach_to_center(g, g0)
      vertex_leg = c(vertex_leg,rep(index,n_nodes))
      index = index +1
    }
    setTxtProgressBar(pb,j)
    
  }
  close(pb)
  if(toPlot)plot(g, vertex.label = NA,vertex.size = 5)
  
  ADJ = get.adjacency(g,attr = "weight")
  ADJ = as.matrix(ADJ)
  
  ADJ[ADJ>0] = 1
  nOBJ = nrow(ADJ)
  
  Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 2000)
  
  toRet = list(ADJ1 = ADJ,ADJ2=ADJ,g=g,Sigma1=Sigma1,Sigma2=Sigma1,vertex_leg=vertex_leg)
  return(toRet)
}

gen_multi_view_dataset = function(nOBJ,nFeat1, nFeat2,S1,S2){
  # This wil be a data matrix wit nFeat1 features and nOBJ samples where the covariance of the sample will follow
  # the structure of the karate network adjacency matrix
  Exp1 <- rmvnorm(nFeat1,rep(0,nOBJ),as.matrix(S1))
  # This wil be a data matrix wit nFeat2 features and nOBJ samples where the covariance of the sample will follow
  # the structure of the karate network adjacency matrix
  Exp2 <- rmvnorm(nFeat2,rep(0,nOBJ),as.matrix(S2))
  
  colnames(Exp1)  = paste("Sample",1:ncol(Exp1),sep="_")
  colnames(Exp2) = paste("Sample",1:ncol(Exp2),sep="_")
  d1 = abs(cor(Exp1)) #Correlation matrix between the samples in the two dataset
  d2 = abs(cor(Exp2))
  
  toRet = list(Exp1=Exp1,Exp2=Exp2,d1=d1,d2=d2)
  return(toRet)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}


clust_analysis = function(HBM3C_mat,sample_names = colnames(Exp1),
                          groupCodes=V(karate)$color,
                          colorCodes = rainbow(2)){
  resultClusterGer <- cluster_gerarchico(list(HBM3C_mat), "complete", 4, 
                                         sample_names, 
                                         sample_names )
  
  plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                             groupCodes=groupCodes,
                             colorCodes=colorCodes)
}

# cluster.method Crea la funzione del cluster da utilizzare per la generazione dell'HBM
HBM3C_function = function(d1,d2,cluster.args = c(),
                           cluster.method  = function(graph) cluster_louvain(graph, weights = (E(graph)$weight)) ){
  listMatrixSimilarityTag = list(d1,d2)
  listDissimilarityMatrix = list(1-d1,1-d2)
  
  print("HBM construction...\n")
  result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
  
  print("HBM3C algorithm...\n")
  HBM3C_mat <- createHBM3C(result_chromo_hbm)
  
  return(HBM3C_mat)
}

plot_HBM3C_mat = function(HBM3C_mat,EIG,toRet,gen_data){
  par(mfrow=c(3,3))
  plot_gray_image(as.matrix(gen_data$ADJ1),main="Original Matrix View1")
  plot_gray_image(as.matrix(gen_data$ADJ2),main="Original Matrix View1")
  
  plot_gray_image(toRet$d1,main="Distance View 1")
  plot_gray_image(toRet$d2,main="Distance View 2")
  
  plot_gray_image(HBM3C_mat,main="HBM matrix")
  
  #distMat = local_reaching_centrality(HBM3C_mat)
  hist(HBM3C_mat,main="HBM values distribution")
  #boxplot(distMat,main = "Local Reaching Centrality")
  plot(EIG[1:10],ylab ="eigenvalue",xlab = "",main="HBM matrix eigenvalues")
  
}

network_analysis = function(HBM3C_mat){
  HBM3C_mat = range01(HBM3C_mat)
  HBM3C_mat = 1-HBM3C_mat
  
  gg = graph.adjacency(adjmatrix =HBM3C_mat,mode = "undirected",weighted = TRUE)
  gg = simplify(gg) 
  
  VNE = vonNewmanEntropy(gg)$VNEntropy
  #EntropyMat[Type,2]=VNE
  cat("Von Newman Entropy: ",VNE,"\n")
  EIG = vonNewmanEntropy(gg)$eigenvalues
  
  toRet = list(gg=gg,VNE=VNE,EIG=EIG)
  return(toRet)
}

complete_analysis = function(gen_data,nFeat1=1000,nFeat2=500,toSave=FALSE,filename,
                             toClust = TRUE,groupCodes,colorCodes,plotLocation){
  
  toRet = gen_multi_view_dataset(nOBJ = nrow(gen_data$ADJ1),nFeat1 = nFeat1,nFeat2 = nFeat2,
                                 S1 = gen_data$Sigma1$mat,S2 = gen_data$Sigma2$mat)
  
  HBM3C_mat = HBM3C_function(toRet$d1,toRet$d2)
  net_res = network_analysis(HBM3C_mat)
  
  png(plotLocation,width=800,height = 600)
  plot_HBM3C_mat(HBM3C_mat,net_res$EIG,toRet,gen_data)
 # par(mfrow=c(1,1))
  if(toClust){
    clust_analysis(HBM3C_mat,sample_names = colnames(toRet$Exp1),
                   groupCodes=groupCodes,
                   colorCodes =colorCodes) 
  }
  
  if(toSave){
    save(toRet,HBM3C_mat,net_res,gen_data,file=filename)
  }
  
  return(list(MVData = toRet,HBM3C_mat=HBM3C_mat,net_res=net_res))

}
