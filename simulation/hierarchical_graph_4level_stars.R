library(igraph)
source("random_graph_functions.R")
source('progetto_tesi/chromoHBM.R')
source('progetto_tesi/chromoHBM3C.R')
source('progetto_tesi/tesi_createHBM.R')
source('progetto_tesi/tesi_createHBM3C.R')
source('progetto_tesi/tesi_cluster_gerarchico.R')
source('progetto_tesi/confusion_matrix.R')
source('progetto_tesi/CMsup_error.R')
source("cormad.R")
source("thresholding.R")

# The initial graph
g = hierarchical_triangle_net_rec(nLevels = 4,module_size = 3,toplot=TRUE)
ADJ = get.adjacency(g,attr = "weight")

plot(g,layout = layout.reingold.tilford(g, root=2),vertex.label = NA,vertex.size = 5)
plot(g,layout = layout.reingold.tilford(g, root=2))

gen_data = gene_data_given_cov_mat(g)
X = gen_data$X
image(X)
image(gen_data$RCor)


d01 = 1-abs(gen_data$RCor)

image(as.matrix(d01))
hls1 = hclust(d = as.dist(d01),method = "ward")
plot(hls1)

listMatrixSimilarityTag = list(1-d01)#list(1-d01)
listDissimilarityMatrix = list(d01)

cluster.args <-c()
cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
#cluster.method <- function(graph) walktrap.community(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
#cluster.method <- function(graph) fastgreedy.community(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM


result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)

image(result_chromo_hbm[[1]]$hbm)
max(result_chromo_hbm[[1]]$hbm)

heatmap(result_chromo_hbm[[1]]$hbm)

result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
resultClusterGer <- cluster_gerarchico(listDissimilarityMatrix, "complete", 3, 
                                       colnames(genes), 
                                       colnames(genes) )
plot(resultClusterGer[[1]]$cluster_gerarchico)


