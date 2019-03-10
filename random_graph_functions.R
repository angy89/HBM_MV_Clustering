range01 <- function(x){(x-min(x))/(max(x)-min(x))}


gene_data_given_cov_mat = function(g){
  # Generate data with a given covaraince matrix:
  # http://stats.stackexchange.com/questions/120179/generating-data-with-a-given-sample-covariance-matrix
  # This is easy to do by generating samples from a standard Gaussian and multiplying them 
  # by a square root of the covariance matrix, e.g. by chol(Σ)
  
  ADJ = get.adjacency(g,attr = "weight")
  library(Matrix)
  n = 100;
  d = nrow(ADJ);
  
  # Make adjacency matrix semi definite positive:
  # http://cstheory.stackexchange.com/questions/6486/condition-to-make-an-adjacency-matrix-positive-semidefinite
  
  Sigma = nearPD(x = ADJ,corr = TRUE)
  SigmaMat = Sigma$mat
  
  # NMat contains d samples (long n) coming from a standard Gaussian distribution
  NMat = matrix(0,n,d)
  for(i in 1:d){
    NMat[,i] = rnorm(n,0,0.1)
  }
  
  CL = chol(SigmaMat)
  
  # X is the matrix with n (genes) and d samples
  X =  NMat %*% CL ;
  
  RCor = cor(as.matrix(X),method = "pearson")
  
  image(RCor)
  image(X)
  toRet = list(X=X,RCor = RCor,CL = CL,ADJ = ADJ,NMat = NMat,Sigma=Sigma)
  return(toRet)
}


hierarchical_star_net_rec = function(nLevels = 1,module_size = 5,toplot=TRUE){
  if(nLevels==1){
    g = make_initial_graph(module_size)
    if(toplot) plot(g,layout = layout.reingold.tilford(g, root=3))
    return(g)
  } 
  
  g = make_initial_graph(module_size)
  end = 8
  starting_child = c()
  to_rem = c()
  
  for(i in 2:nLevels){
    nsn = 4^(i-1)
    nodes = 1:vcount(g)
    to_rem = c(to_rem,starting_child)
    sum = 3 + length(to_rem)
    for(j in 1:(nsn/4)){
      to_rem = c(to_rem,sum)
      sum = sum + module_size
    }
    
    starting_child = nodes[-to_rem]
    cat(starting_child,"\n")
    for(sc in starting_child){
      start = sc
      g = g + make_initial_graph(module_size) + edges(c(start,end))
      end = end + module_size
    }
  }
  if(toplot) plot(g,layout = layout.reingold.tilford(g, root=3))
  
  return(g)
}

hierarchical_triangle_net_rec2 = function(nLevels = 1,module_size = 3,toplot=TRUE){
  corrs = seq(0.5, 1, length.out = nLevels+1)
  
  if(nLevels==1){
    g = make_initial_graph(module_size,parent_cor = list(min = corrs[1],max=corrs[2]),
                           between_cor = list(min=corrs[1],max = corrs[2]))
    if(toplot) plot(g,layout = layout.reingold.tilford(g, root=2))
    return(g)
  } 
  
  g = make_initial_graph(module_size,parent_cor = list(min = corrs[1],max=corrs[2]),
                         between_cor = list(min=corrs[1],max = corrs[2]))
  end = 5
  starting_child = c()
  to_rem = c()
  
  for(i in 2:nLevels){
    nsn = 2^(i-1)
    nodes = 1:vcount(g)
    to_rem = c(to_rem,starting_child)
    sum = 2 + length(to_rem)
    for(j in 1:(nsn/2)){
      to_rem = c(to_rem,sum)
      sum = sum + module_size
    }
    
    starting_child = nodes[-to_rem]
    cat(starting_child,"\n")
    for(sc in starting_child){
      start = sc
      gadd = make_initial_graph(module_size,parent_cor = list(min = corrs[i],max=corrs[i+1]),
                                between_cor = list(min=corrs[i],max = corrs[i+1]))
      g = g + gadd + edges(c(start,end),weight=runif(1,corrs[i],corrs[i+1]))
      
      #edge_mat = rbind(start, gorder(g) + 1:gorder(gadd))
      #n_edges = ncol(edge_mat)
      
      # edge_mat2 = rbind(to_rem, gorder(g) + 1:gorder(gadd))
      #n_edges2 = ncol(edge_mat2)
      #g <- g + gadd + edges(as.vector(edge_mat),weight=runif(n_edges2,corrs[i],corrs[i+1]))# + 
      # edges(as.vector(edge_mat2),weight=runif(n_edges2,corrs[i-1],corrs[i]))
      
      end = end + module_size
    }
  }
  if(toplot) plot(g,layout = layout.reingold.tilford(g, root=2))
  
  return(g)
}

hierarchical_triangle_net_rec = function(nLevels = 1,module_size = 3,toplot=TRUE){
  corrs = seq(0.7, 1, length.out = nLevels+1)
  
  if(nLevels==1){
    g = make_initial_graph(module_size,parent_cor = list(min = corrs[1],max=corrs[2]),
                           between_cor = list(min=corrs[1],max = corrs[2]))
    if(toplot) plot(g,layout = layout.reingold.tilford(g, root=2))
    return(g)
  } 
  
  g = make_initial_graph(module_size,parent_cor = list(min = corrs[1],max=corrs[2]),
                              between_cor = list(min=corrs[1],max = corrs[2]))
  end = 5
  starting_child = c()
  to_rem = c()
  
  for(i in 2:nLevels){
    nsn = 2^(i-1)
    nodes = 1:vcount(g)
    to_rem = c(to_rem,starting_child)
    sum = 2 + length(to_rem)
    for(j in 1:(nsn/2)){
      to_rem = c(to_rem,sum)
      sum = sum + module_size
    }
    
    starting_child = nodes[-to_rem]
    cat(starting_child,"\n")
    for(sc in starting_child){
      start = sc
      gadd = make_initial_graph(module_size,parent_cor = list(min = corrs[i],max=corrs[i+1]),
                                between_cor = list(min=corrs[i],max = corrs[i+1]))
      g = g + gadd + edges(c(start,end),weight=runif(1,corrs[i],corrs[i+1]))
      
      #edge_mat = rbind(start, gorder(g) + 1:gorder(gadd))
      #n_edges = ncol(edge_mat)
      
     # edge_mat2 = rbind(to_rem, gorder(g) + 1:gorder(gadd))
      #n_edges2 = ncol(edge_mat2)
      #g <- g + gadd + edges(as.vector(edge_mat),weight=runif(n_edges2,corrs[i],corrs[i+1]))# + 
                     # edges(as.vector(edge_mat2),weight=runif(n_edges2,corrs[i-1],corrs[i]))
      
      end = end + module_size
    }
  }
  if(toplot) plot(g,layout = layout.reingold.tilford(g, root=2))
  
  return(g)
}








#' This adds the gadd graph to the main graph, g, and wires all of its vertices
#' to the central vertex of g
attach_to_center <- function(idx = 1,g, gadd,min_cor=0.4,max_cor=0.5) {
  #runif(25,0.4,0.5)
  edge_mat = rbind(idx, gorder(g) + 1:gorder(gadd))
  n_edges = ncol(edge_mat)
  g <- g + gadd + edges(as.vector(edge_mat),weight=runif(n_edges,min_cor,max_cor))
}


random_graph = function(min_cor = 0.7,n_islands = 3,min_island_size = 3,max_island_size = 7){
  island_size = ceiling(sample(size = n_islands,x = min_island_size:max_island_size))
  
  #island_size = ceiling(runif(n = n_islands,min = min_island_size,max = max_island_size))
  n_nodes = sum(island_size)
  
  edges_between = 0.1
  
  #For preferential attachment layout
  #n_edges_between_islands = ceiling(((n_nodes * n_nodes) - (n_nodes - n_islands))*edges_between)
  
  
  #Each islands of the network is generated as a preferential attachment model
  ADJ = matrix(0,nrow = n_nodes,ncol = n_nodes)
  start = 1
  end = island_size[1]
  islands = list()
  
  for(i in 1:n_islands){
    #G_pa= sample_pa(island_size[i],directed = TRUE)
    #E(G_pa)$weight = runif(n = length(E(G_pa)), min = min_cor,max = 1)
    #sub_adj = get.adjacency(G_pa,attr = "weight",type = "both",sparse = FALSE)
    
    sub_adj = matrix(0,island_size[i],island_size[i])
    for(row_ in 1:island_size[i]){
      for(col_ in row_:island_size[i]){
        sub_adj[row_,col_] = sub_adj[col_,row_] = runif(1,min = min_cor,max = 1)
      }
    }
    
    diag(sub_adj) = 0
    
    islands[[i]]= start:end
    ADJ[start:end,start:end] = sub_adj
    cat(start, " ",end ,"\n")
    start = end+1
    end = start + island_size[i+1] -1
  }
  
  GG = graph.adjacency(ADJ,mode="undirected",weighted = TRUE)
  return(GG)
  
  
}

#Create initial graph of 5 nodes
make_initial_graph = function(n_nodes=5, parent_cor  = list(min=0.6,max=0.7), 
                                       between_cor = list(min = 0.9,max=0.9)){
  ADJ = matrix(0,n_nodes,n_nodes)
  ADJ[1,2:n_nodes] = ADJ[2:n_nodes,1] = runif(n = n_nodes -1, min = parent_cor$min,max = parent_cor$max)
  
  for(i in 2:(n_nodes-1)){
    for(j in (i+1):n_nodes){
      ADJ[i,j] = ADJ[j,i] = runif(n = 1,min = between_cor$min,max = between_cor$max)
    }
  }
  
  GG = graph.adjacency(ADJ,mode="undirected",weighted = TRUE)
  return(GG)
  
}
