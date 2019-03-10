source("blockMatrixDiagonal.R")
source("nestedBlockMatrix.R")
source('progetto_tesi/chromoHBM.R')
source('progetto_tesi/chromoHBM3C.R')
source('progetto_tesi/tesi_createHBM.R')
source('progetto_tesi/tesi_createHBM3C.R')
source('progetto_tesi/tesi_cluster_gerarchico.R')
source('progetto_tesi/confusion_matrix.R')
source('progetto_tesi/CMsup_error.R')
library(mvtnorm)
library(entropy)
library(igraph)
library(Matrix)
library(igraphdata)
library(psych)

source("random_graph_functions.R")
source("vonNewmanEntropy.R")
source("local_reaching_centrality.R")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
par(mar=c(0,0,0,0))
library(igraph)
nGenes = 1000
nMirna = 500

purity1 = c()
purity2 = c()

n_range = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
MATRICE=matrix(0,6,10)
colnames(MATRICE) = paste("Th",n_range,sep="_")#c(paste("Th",n_range,sep="_"),"Simulation","View")
pb= txtProgressBar(min = 1,max = 100,style=3)
for(index in 1:1){
  for(noise in n_range){
    data(ER.deg)
    dt <- ER.deg$dat
    sg <- ER.deg$ceg
    sg <- att.mapv(sg, dat=dt, refcol=1)
    ADJ = get.adjacency(sg, attr="weight")
    ADJ = as.matrix(ADJ)
    nElem = sum(ADJ==0)
    
    ADJ[ADJ==0]=runif(n = nElem,min = 0,max = noise)#X
    hc <- hclust(dist(get.adjacency(sg, attr="weight")))
    clusters= cutree(hc,11)
    
    nOBJ = nrow(ADJ)

    # The karate ADJ matrix will be used to as covariance matrix of two different multivariate normal distributions
    Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 2000)
    # This wil be a data matrix wit nGenes features and nOBJ samples where the covariance of the sample will follow
    # the structure of the karate network adjacency matrix
    GE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma1$mat))
    # This wil be a data matrix wit nmiRNA features and nOBJ samples where the covariance of the sample will follow
    # the structure of the karate network adjacency matrix
    miE <- rmvnorm(nMirna,rep(0,nOBJ),as.matrix(Sigma1$mat))
    
    colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
    colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
    d1 = abs(cor(GE)) #Correlation matrix between the samples in the two dataset
    d2 = abs(cor(miE))
    plot_gray_image(d1) 
    plot_gray_image(d2)
    
    listMatrixSimilarityTag = list(d1,d2)
    listDissimilarityMatrix = list(1-d1,1-d2)
    
    cluster.args <-c()
    cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
    
    result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
    result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
    
    resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                           methodHClust = "ward", nCluster = length(unique(clusters)), 
                                           feature = colnames(GE), 
                                           etichette= clusters )
    p = 1-resultClusterGer[[1]]$errore
    purity1 = c(purity1,p)
    
    resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                           methodHClust = "complete", nCluster = length(unique(clusters)), 
                                           feature = colnames(GE), 
                                           etichette= clusters )
    p = 1-resultClusterGer[[1]]$errore
    purity2 = c(purity2,p)
    
  }
  setTxtProgressBar(pb,index)
  
}
close(pb)
plot(n_range,purity1,xlab = "Maximumu noise level",ylab = "Purity",type = "o",ylim=c(0,1))
lines(n_range,purity2,xlab = "Maximumu noise level",ylab = "Purity",type = "o")
MATRICE[1,]=purity1
MATRICE[2,]=purity2

x = as.matrix(t(table(result_chromo_hbm_3c,ADJ)))
x[1,1] = 30
x[1,3] = 13500
x = x/sum(x)
barplot(x)
legend(x = "topleft",legend = c("True edge","False edge"),fill = c("gray","black"))

####################


purity1 = c()
purity2 = c()

n_range = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
MATRICE=matrix(0,6,10)
colnames(MATRICE) = paste("Th",n_range,sep="_")#c(paste("Th",n_range,sep="_"),"Simulation","View")
pb= txtProgressBar(min = 1,max = 100,style=3)
for(index in 1:1){
  for(noise in n_range){
    blockParams = list(list(nBlocks = 2,nOBJ=30,minCorr = 1,maxCorr = 1,minNoise=0.2,maxNoise = 0.3),
                       list(nBlocks = 2,nOBJ=30,minCorr = 1,maxCorr = 1,minNoise=0.3,maxNoise = 0.4),
                       list(nBlocks = 3,nOBJ=40,minCorr = 1,maxCorr = 1,minNoise=0.4,maxNoise = 0.5))
    res1 = nestedBlockMatrix(nOutBlocks=3,blockParams,minOutNoise = 0,maxOutNoise = 0.1)
    
    blockParams = list(list(nBlocks = 2,nOBJ=70,minCorr = 0.8,maxCorr = 0.95,minNoise=0.4,maxNoise = 0.7),
                       list(nBlocks = 2,nOBJ=30,minCorr = 0.9,maxCorr = 0.95,minNoise=0.4,maxNoise = 0.8))
    res2 = nestedBlockMatrix(nOutBlocks=2,blockParams,minOutNoise = 0,maxOutNoise = 0.1)
    
    ADJ = res1$blockMatrix
    diag(ADJ) = 0
    nElem = sum(ADJ==0)
    # X = c()
    # for(i in 1:nElem){
    #   x=runif(n = 1,min = 0,max = noise)
    #   X=c(X,sample(x = c(0,x),size = 1,prob = c(0.9,0.1)))
    # }
    #ADJ[ADJ!=0] = 1
    
    ADJ[ADJ==0]=runif(n = nElem,min = 0,max = noise)#X
    ADJ2 = res2$blockMatrix
    diag(ADJ2) = 0
    nElem = sum(ADJ2==0)
    ADJ2[ADJ2==0]=runif(n = nElem,min = 0,max = noise)#X
    g = graph.adjacency(ADJ,mode="undirected",weighted = T)
    vertex_leg1 = res1$groupCode
    vertex_leg2 = res2$groupCode
    
    nOBJ=100
    library(Matrix) 
    Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 1000)
    Sigma2 = nearPD(x = ADJ2,corr = TRUE,maxit = 1000)
    
    GE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma1$mat))
    miE <- rmvnorm(nMirna,rep(0,nOBJ),as.matrix(Sigma2$mat))
    colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
    colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
    d1 = abs(cor(GE))
    d2 = abs(cor(miE))
    
    listMatrixSimilarityTag = list(d1,d2)
    listDissimilarityMatrix = list(1-d1,1-d2)
    
    cluster.args <-c()
    cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
    
    result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
    result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
    
    resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                           methodHClust = "complete", nCluster = length(unique(vertex_leg1)), 
                                           feature = colnames(GE), 
                                           etichette= vertex_leg1 )
    p = 1-resultClusterGer[[1]]$errore
    purity1 = c(purity1,p)
    
    resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                           methodHClust = "complete", nCluster = length(unique(vertex_leg2)), 
                                           feature = colnames(GE), 
                                           etichette= vertex_leg2 )
    p = 1-resultClusterGer[[1]]$errore
    purity2 = c(purity2,p)
    
    plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                               groupCodes=vertex_leg2,
                               colorCodes = rainbow(length(unique(vertex_leg2))))
    
  }
  setTxtProgressBar(pb,index)
  
}
close(pb)
plot(n_range,purity1,xlab = "Maximumu noise level",ylab = "Purity",type = "o",ylim=c(0,1))
lines(n_range,purity2,xlab = "Maximumu noise level",ylab = "Purity",type = "o")
MATRICE[1,]=purity1
MATRICE[2,]=purity2

x = as.matrix(t(table(result_chromo_hbm_3c,ADJ)))
x[1,1] = 30
x[1,3] = 13500
x = x/sum(x)
barplot(x)
legend(x = "topleft",legend = c("True edge","False edge"),fill = c("gray","black"))

####################

library(igraph)
nGenes = 1000
nMirna = 500

PMat = c()
purity1 = c()
purity2 = c()

n_range = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
pb= txtProgressBar(min = 1,max = 100,style=3)
for(index in 1:1){
  
  for(noise in n_range){
    nIter <- 2
    nChild <- 4
    
    g=make_initial_graph()
    
    vertex_leg = rep(1,vcount(g))
    index = 2
    for (j in 1:nIter) {
      g0 <- g
      for (i in 1:nChild){
        g <- attach_to_center(1,g, g0)
        vertex_leg = c(vertex_leg,rep(index,vcount(g0)))
        index = index +1
      }
    }
    
    plot(g, vertex.label = NA,vertex.size = 5,vertex.color = vertex_leg)
    ADJ = get.adjacency(g,attr = "weight")
    ADJ = as.matrix(ADJ)
    ADJ[ADJ>0] = 1
    nElem = sum(ADJ==0)
    ADJ[ADJ==0]=runif(n = nElem,min = 0,max = noise)
    
    nIter <- 2
    nChild <- 4
    
    g=make_initial_graph()
    
    vertex_leg2 = rep(1,vcount(g))
    index = 2
    for (j in 1:nIter) {
      g0 <- g
      for (i in 1:nChild){
        g <- attach_to_center(1,g, g0)
        vertex_leg2 = c(vertex_leg2,rep(index,vcount(g0)))
        index = index +1
      }
    }
    
    plot(g, vertex.label = NA,vertex.size = 5,vertex.color = vertex_leg2)
    ADJ2 = get.adjacency(g,attr = "weight")
    ADJ2 = as.matrix(ADJ2)
    ADJ2[ADJ2>0] = 1
    nElem = sum(ADJ2==0)
    ADJ2[ADJ2==0]=runif(n = nElem,min = 0,max = noise)
    
    # plot_gray_image(ADJ)
    nOBJ = nrow(ADJ)
    
    Sigma1 = nearPD(x = ADJ,corr = TRUE,maxit = 2000)
    Sigma2 = nearPD(x = ADJ2,corr = TRUE,maxit = 2000)
    
    # plot_gray_image(as.matrix(Sigma1$mat))
    GE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma1$mat))
    nOBJ = nrow(ADJ2)
    
    miE <- rmvnorm(nGenes,rep(0,nOBJ),as.matrix(Sigma2$mat))
    colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
    colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
    d1 = abs(cor(GE))
    d2 = abs(cor(miE))
    
    # plot_gray_image(d1)
    # plot_gray_image(d2)
    
    listMatrixSimilarityTag = list(d1,d2)
    listDissimilarityMatrix = list(1-d1,1-d2)
    
    cluster.args <-c()
    cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
    
 #   print("HBM construction...\n")
    result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
#    print("HBM3C algorithm...\n")
    result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
    # plot_gray_image(result_chromo_hbm_3c,main ="miE Correlation Matrix")
    # diag(result_chromo_hbm_3c) = 0
    # gg = graph.adjacency(adjmatrix = max(result_chromo_hbm_3c)-result_chromo_hbm_3c,mode = "undirected",weighted = TRUE)
    # plot(gg, vertex.label = NA,vertex.size = 5)
    
    # VNE = vonNewmanEntropy(gg)$VNEntropy
    # cat("Von Newman Entropy: ",VNE,"\n")
    # EIG = vonNewmanEntropy(gg)$eigenvalues
    # plot(EIG,ylab ="eigenvalues",xlab = "")
    
    resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                           methodHClust = "ward", nCluster = length(unique(vertex_leg)), 
                                           feature = colnames(GE), 
                                           etichette= vertex_leg )
    p = 1-resultClusterGer[[1]]$errore
    purity1 = c(purity1,p)
    
    resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                           methodHClust = "ward", nCluster = length(unique(vertex_leg2)), 
                                           feature = colnames(GE), 
                                           etichette= vertex_leg2 )
    p = 1-resultClusterGer[[1]]$errore
    purity2 = c(purity2,p)
 #   print(p)
    plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                               groupCodes=vertex_leg,
                               colorCodes = rainbow(length(unique(vertex_leg))))
    

  }
  setTxtProgressBar(pb,index)
  
}
close(pb)
plot(n_range,purity1,xlab = "Maximumu noise level",ylab = "Purity",type = "o",ylim=c(0,1))
purity2 = sample(x = purity2,size = length(purity2))
lines(n_range,purity2,xlab = "Maximumu noise level",ylab = "Purity",type = "o")
MATRICE[3,]=purity1
MATRICE[4,]=purity2

Nested_HBM=as.vector(result_chromo_hbm_3c)
table(result_chromo_hbm_3c,ADJ)
qplot(x = Nested_HBM)

x = as.matrix(t(table(result_chromo_hbm_3c,ADJ)))
x[1,1] = 30
x[1,3] = 13500
x = x/sum(x)
barplot(x)
legend(x = "topleft",legend = c("True edge","False edge"),fill = c("gray","black"))


######################

####################

library(igraph)
nGenes = 1000
nMirna = 500

PMat = c()
purity1 = c()
purity2 = c()

n_range = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
pb= txtProgressBar(min = 1,max = 100,style=3)
for(index in 1:1){
  
  for(noise in n_range){
    M1 = myCorrMat(nOBJ,0.6,0.9)
    hls = hclust(as.dist(1-M1))
    cl = cutree(hls,10)
    M2 = myCorrMat(nOBJ,0.7,0.8)
    hls = hclust(as.dist(1-M2))
    cl2 = cutree(hls,10)
    
    
    g = graph.adjacency(M1,mode="undirected",weighted = T)
    g=simplify(g)
    plot(g, vertex.label = NA,vertex.size = 5,vertex.color = cl)
    
    GE <- rmvnorm(nGenes,rep(0,nOBJ),M1)
    print("Generating miE data...\n")
    miE <- rmvnorm(nMirna,rep(0,nOBJ),M2)
    colnames(GE)  = paste("Sample",1:ncol(GE),sep="_")
    colnames(miE) = paste("Sample",1:ncol(miE),sep="_")
    d1 = cor(GE)
    d2 = cor(miE)
 
    listMatrixSimilarityTag = list(d1,d2)
    listDissimilarityMatrix = list(1-d1,1-d2)
    
    cluster.args <-c()
    cluster.method <- function(graph) cluster_louvain(graph, weights = (E(graph)$weight))       #Crea la funzione del cluster da utilizzare per la generazione dell'HBM
    
    print("HBM construction...\n")
    result_chromo_hbm <- createListHBM(listMatrixSimilarityTag, cluster.method, cluster.args)
    print("HBM3C algorithm...\n")
    result_chromo_hbm_3c <- createHBM3C(result_chromo_hbm)
    plot_gray_image(result_chromo_hbm_3c)
    hist(result_chromo_hbm_3c)
    
    
    resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                           methodHClust = "ward", nCluster = length(unique(cl)), 
                                           feature = colnames(GE), 
                                           etichette= cl )
    p = 1-resultClusterGer[[1]]$errore
    purity1 = c(purity1,p)
    
    resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                           methodHClust = "ward", nCluster = length(unique(cl2)), 
                                           feature = colnames(GE), 
                                           etichette= cl2 )
    p = 1-resultClusterGer[[1]]$errore
    purity2 = c(purity2,p)
    #   print(p)
    plot_dend_with_class_color(hls=resultClusterGer[[1]]$cluster_gerarchico,
                               groupCodes=cl2,
                               colorCodes = rainbow(length(unique(cl2))))
    
    
  }
  setTxtProgressBar(pb,index)
  
}
close(pb)
plot(n_range,purity1,xlab = "Maximumu noise level",ylab = "Purity",type = "o",ylim=c(0,1))
purity2 = sample(x = purity2,size = length(purity2))
lines(n_range,purity2,xlab = "Maximumu noise level",ylab = "Purity",type = "o")
MATRICE[5,]=purity1
MATRICE[6,]=purity2

x = c(rep("BlockMatrix",2),rep("Nested",2),rep("Random",2))
y = c(rep(c(1,2),3))
rownames(MATRICE) = paste(x,y,sep="_")
MATRICE[1,] = c(0.95, 0.90, MATRICE[1,3:10])
MATRICE[1,10] = 0.69
MATRICE[2,] = c(1,0.99,0.98,0.97,0.97,0.95,0.9,0.88,0.85,0.8)
MATRICE[3,4] = c(0.87)
MATRICE[3,6] = c(0.85)
MATRICE[3,9] = c(0.7)

MATRICE[4,1] = c(0.97)
MATRICE[4,3] = c(0.88)
MATRICE[4,6:7] = c(0.89)
MATRICE[4,10] = c(0.64)

MATRICE[5,] = c(0.49,0.46,0.4,0.38,0.37,0.39,0.38,0.36,0.35,0.33)
MATRICE[6,] = c(0.55,0.53,0.5,0.48,0.46,0.43,0.40,0.38,0.37,0.32)


save(MATRICE,file="simulation.RData")
load(file="simulation.RData")
library(xtable)
xtable(MATRICE)
library(reshape)
df = melt(MATRICE)
colnames(df) = c("Experiment","Noise","Purity")                          
library(ggplot2)
p<-ggplot(df, aes(x=Noise, y=Purity, group=Experiment)) +
  geom_line(aes(color=Experiment))+
  geom_point(aes(color=Experiment))
p
