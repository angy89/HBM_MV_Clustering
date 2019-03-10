local_reaching_centrality = function(HBM3C_mat = result_chromo_hbm_3c ) {
  levels = sort(unique(as.vector(HBM3C_mat)))
  distMat = matrix(0,nrow(HBM3C_mat),ncol=max(HBM3C_mat))
  diag(HBM3C_mat) = 0
  
  for(i in 1:nrow(HBM3C_mat)){
    for(j in levels){
      distMat[i,j] = sum(HBM3C_mat[i,]==j)/nrow(HBM3C_mat)
    }
  }
  
  return(distMat)
}