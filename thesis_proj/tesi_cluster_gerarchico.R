#effettuo la tecnica del cluster gerarchico per ogni matrice appartenente alla lista
cluster_gerarchico<-function(listMatrix,methodHClust,nCluster, feature, etichette){
  errore <-c()
  confusion_matrix <- list()
  cluster_MI <- c()
  clsGer<-c()
  for(i in 1:length(listMatrix)){
    dist_clust<- as.dist(listMatrix[[i]])
    clust_ger<-hclust(dist_clust, method = methodHClust)
    pred<-cutree(clust_ger, k = nCluster)       #il numero dei cluster generati uguale al numero di cluster generati dal metodo MVDA                
    pred <- setNames(pred, feature)
    CM<-confusion_matrix(etichette = etichette, clustering = pred, pazienti_geni = NULL ,nCluster = dim(table(pred)))
    errore <- CMsup_error(CM_sup = CM,nPat = dim(listMatrix[[i]])[2])
    clsGer[[i]] <- list('cluster_gerarchico'=clust_ger, 'pred'=pred, 'confusion_matrix' = confusion_matrix, 'errore'= errore)
  }
  
  return(clsGer)
  
}