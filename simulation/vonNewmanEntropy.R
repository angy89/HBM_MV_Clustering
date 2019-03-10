vonNewmanEntropy = function(gg){
  #compute the normalized laplacian of the graph

  # c = 1/(2*ecount(gg))
  # D = diag(x = degree(gg),nrow = vcount(gg),ncol = vcount(gg))
  # A = get.adjacency(graph = gg,sparse = FALSE,attr = "weight")
  # A[A>0] = 1
  # 
  # L = c * (D - A)
  # tr(L)
  # eigenvalues = eigen(L)$values
  # eigenvalues = eigenvalues[eigenvalues>0]
  # VNEntropy = - sum(eigenvalues * log2(eigenvalues))
  # 
  # return(list(laplacian = L, eigenvalues = eigenvalues, VNEntropy = VNEntropy))
  
  lap = graph.laplacian(gg,normalized = TRUE)
  glap = graph.adjacency(lap,weighted = T,mode="undirected") 
  #compute the spectrum of the laplacian graph
  eigenvalues = spectrum(glap,which = list(pos = "LA",howmany = vcount(gg)-1))$values
  eigenvalues = eigenvalues[eigenvalues>0]

  vne = 0
  for(i in eigenvalues){
    vne = vne + (i * log2(i))
  }
  return(list(laplacian = lap, eigenvalues = eigenvalues, VNEntropy = vne))
  
}