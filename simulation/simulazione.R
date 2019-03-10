#simulazione dati 

source('progetto_tesi/chromoHBM.R')
source('progetto_tesi/chromoHBM3C.R')
source('progetto_tesi/tesi_createHBM.R')
source('progetto_tesi/tesi_createHBM3C.R')
source('progetto_tesi/tesi_cluster_gerarchico.R')
source('progetto_tesi/confusion_matrix.R')
source('progetto_tesi/CMsup_error.R')

range01 <- function(x){(x-min(x))/(max(x)-min(x))}
sim_gene = function(nGenes,nsamples){
  M = c()
  for(i in 1:nSamples){
    M = cbind(M,rnorm(nGenes,mean = 0, sd = 1))
  }
  
  for(i in 1:nSamples){
    M = cbind(M,rnorm(nGenes,mean = 0, sd = 3))
  }
  
  for(i in 1:nSamples){
    M = cbind(M,rnorm(nGenes,mean = 4, sd =1))
  }
  
  colnames(M) = c(rep("Mean0",nSamples*2),rep("Mean10",nSamples))
  return(M)
}

sim_mir = function(nGenes,nsamples){
  M = c()
  for(i in 1:nSamples){
    M = cbind(M,rnorm(nGenes,mean = 1, sd = 1))
  }
  
  for(i in 1:nSamples){
    M = cbind(M,rnorm(nGenes,mean = 1, sd = 3))
  }
  
  for(i in 1:nSamples){
    M = cbind(M,rnorm(nGenes,mean = 3, sd =1))
  }
  
  colnames(M) = c(rep("Mean0",nSamples*2),rep("Mean10",nSamples))
  return(M)
}



nGenes = 1000
nSamples = 100
genes = sim_gene(nGenes,nSamples)
miRNA = sim_mir(nGenes,nSamples)

d1 = dist(t(noise_gene))
d2 = dist(t(miRNA))

image(as.matrix(d01))
hls1 = hclust(d = d1,method = "complete")
plot(hls1)
hls2 = hclust(d = d2,method = "complete")
plot(hls2)

d01 = range01(as.matrix(d1))
d02 = range01(as.matrix(d2))

listMatrixSimilarityTag = list(1-d01,1-d02)
listDissimilarityMatrix = list(d01,d02)

cluster.args <-c()
cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM


result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
resultClusterGer <- cluster_gerarchico(listDissimilarityMatrix, "complete", 3, 
                                       colnames(genes), 
                                       colnames(genes) )
plot(resultClusterGer[[1]]$cluster_gerarchico)
