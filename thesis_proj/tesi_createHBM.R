createListHBM<-function(listMatrix,clusterMethod,clusterArgs){
  listHBM <- c()
  for(i in 1:length(listMatrix)){
    listHBM[[i]]<-chromoHBM(listMatrix[[i]], clusterMethod, clusterArgs)
  }
  return(listHBM)
}