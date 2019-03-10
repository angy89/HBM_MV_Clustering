createMatrixSimilarity<-function(listMatrix, method,indexSim){
  listSimMatrix <-c()
  listMatrixT <- c()
  for(i in 1:length(listMatrix)){
    listSimMatrix[[i]]<-abs(cor(listMatrix[[i]],method = method))
    listMatrixT[[i]] = listSimMatrix[[i]]
    listMatrixT[[i]][listMatrixT[[i]]<indexSim] = 0
  }
  result_matrix_similarity <- list('listSimMatrix'=listSimMatrix,'listMatrixT'=listMatrixT)
  return(result_matrix_similarity)
}