#hierarchical simulator

library(igraph)
min_cor = 0.7
between_min_cor = 0.2
between_max_cor = 0.3
min_island_size = 3
max_island_size = 7
n_islands = 3

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


n_edges_between_islands = ceiling(sum(ADJ[lower.tri(ADJ,diag = FALSE)] == 0) * edges_between)
  

GG = graph.adjacency(ADJ,mode="undirected",weighted = TRUE)
plot(GG)

#Some random edges are attacched between the islands;
# The edges inside the islands and the edges between the islands are added as strong edges

for(i in 1:n_edges_between_islands){
  idx_island_from = sample(x = 1:length(islands),size = 1)
  from = islands[[idx_island_from]][sample(1:length(islands[[idx_island_from]]),1)]
  sample_to=1:length(islands)
  sample_to = sample_to[-idx_island_from]
  idx_island_to = sample(x = sample_to,size = 1)
  to = islands[[idx_island_to]][sample(1:length(islands[[idx_island_to]]),1)]
  ADJ[from, to]= ADJ[to,from] = rnorm(1, between_min_cor, between_max_cor)
}


for(i in 1:(n_nodes-1)){
  for(j in (i+1):n_nodes){
    if(ADJ[i,j] == 0){
      prob = runif(1,0,1)
      if(prob>0.5){
        ADJ[i,j] = ADJ[j,i] = runif(1,min = 0,max = 0.1)
      }
    }
  }
}

GG = graph.adjacency(ADJ,mode="undirected",weighted = TRUE)
par(mar = c(0,0,0,0))
GG$layout=layout.fruchterman.reingold(GG,weights=E(GG)$weight)

plot(GG,edge.arrow.size = 0.001,edge.width = E(GG)$weight * 5,mark.groups = islands,vertex.size = 5,vertex.label.cex = 0.1)

# Generate data with a given covaraince matrix:
# http://stats.stackexchange.com/questions/120179/generating-data-with-a-given-sample-covariance-matrix
# This is easy to do by generating samples from a standard Gaussian and multiplying them 
# by a square root of the covariance matrix, e.g. by chol(Σ)

library(Matrix)
n = 100;
d = nrow(ADJ);

# Make adjacency matrix semi definite positive:
# http://cstheory.stackexchange.com/questions/6486/condition-to-make-an-adjacency-matrix-positive-semidefinite

Sigma = nearPD(x = ADJ,corr = TRUE)

# NMat contains d samples (long n) coming from a standard Gaussian distribution
NMat = matrix(0,n,d)
for(i in 1:d){
  NMat[,i] = rnorm(n)
}

CL = chol(Sigma$mat)

# X is the matrix with n (genes) and d samples
X =  NMat %*% CL ;

RCor = cor(as.matrix(X),method = "pearson")

