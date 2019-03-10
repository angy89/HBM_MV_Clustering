MVDA_clustering<-function(miRNA, genes, info){
  
  miRNA_km <- kmeans_preprocessing(DB = miRNA, nCenters = 10)
  
  miRNA_km_summary <- clustering_summary(DB=miRNA, cluster=miRNA_km$clustering)
  
  genes_pamk <- pamk_preprocessing(DB = genes,nCenters = 50)
  
  genes_pamk_summary <- clustering_summary(DB=genes, cluster=genes_pamk$clustering)
  
  
  # 2 step prototype ranking
  
  sda_km_miRNA <- sda_ranking(prototype = miRNA_km$centers,info = info$x,ranking.score = "avg")
  
  plot(sda_km_miRNA)
  
  sda_genes <- sda_ranking(prototype = genes_pamk$centers,info = info$x,ranking.score = "avg")
  
  plot(sda_genes)
  
  
  #3 step single view clustering
  
  prototypes_miRNA <- miRNA_km$centers[rownames(sda_km_miRNA),]
  
  nCenter <- 4
  
  miRNA_pat_km <- kmeans_sv(nCenters = nCenter,prototype = prototypes_miRNA)
  
  cm_km <- confusion_matrix(etichette = info$x,
                            
                            clustering = miRNA_pat_km$clustering,
                            
                            pazienti_geni = t(prototypes_miRNA),
                            
                            nCluster = nCenter)
  
  miRNAerror <- CMsup_error(CM_sup = cm_km,nPat = dim(prototypes_miRNA)[2])
  
  prototypes_genes <- genes_pamk$centers[rownames(sda_genes),]
  
  nCenter <- 4
  
  genes_pat_km <- kmeans_sv(nCenters = nCenter,prototype = prototypes_genes)
  
  cm_km <- confusion_matrix(etichette = info$x,clustering = genes_pat_km$clustering,
                            
                            pazienti_geni = t(prototypes_genes),nCluster = nCenter)
  
  geneserror <- CMsup_error(CM_sup = cm_km,nPat = dim(prototypes_genes)[2])
  
  nrows <- dim(info)[1]
  
  ncols <- length(table(miRNA_pat_km$clustering))
  
  rows_id <- rownames(info)
  
  cols_id <- paste("Clustering",names(table(miRNA_pat_km$clustering)),sep=" ")
  
  X1 <- supervised_matrix(etichette = miRNA_pat_km$clustering,nRow = nrows,
                          
                          nCol = ncols,row_id = rows_id,col_id  = cols_id)
  
  X2 <- supervised_matrix(etichette = genes_pat_km$clustering,nRow = nrows,
                          
                          nCol = ncols,row_id = rows_id,col_id  = cols_id)
  
  SUP <- supervised_matrix(etichette = as.numeric(info$x),nRow = nrows,
                           
                           nCol = length(table(info$x)),
                           
                           row_id = rows_id,col_id  = cols_id)
  
  X <- cbind(X1,X2,SUP)
  
  K <- dim(X)[2]
  
  patients <- t(cbind(t(prototypes_genes),t(prototypes_miRNA)))
  
  # mf_res <- MF(A = X,k = K,eps = 0.01,iter_max = 100,nView = 2,V1.lc1 = 4,
  #              
  #              V1.lc2 = 4,V1.lc3 = 4,V1.lc4 = 0,V1.lc5 = 0,
  #              
  #              info = info,nclust_s = 4,pazienti = patients)
  
  GLI_res <- general_late_integration(A = X,k = K,alfa = 1, eps = 0.001,pazienti = patients,info = info)
  
  return(GLI_res)
  
}
