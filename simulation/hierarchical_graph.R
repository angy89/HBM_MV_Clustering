library(igraph)
source("random_graph_functions.R")

nIter <- 2
nChild <- 4

# The initial graph
#g <- make_empty_graph(5, directed = FALSE) + edges(1,2,1,3,1,4,1,5,2,3,3,4,4,5,5,2)
g=make_initial_graph()

for (j in 1:nIter) {
  g0 <- g
  for (i in 1:nChild){
    g <- attach_to_center(g, g0)
  }
}

plot(g, vertex.label = NA,vertex.size = 5)

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

# NMat contains d samples (long n) coming from a standard Gaussian distribution
NMat = matrix(0,n,d)
for(i in 1:d){
  NMat[,i] = rnorm(n)
}

CL = chol(Sigma$mat)

# X is the matrix with n (genes) and d samples
X =  NMat %*% CL ;

RCor = cor(as.matrix(X),method = "pearson")

image(X)



