library(RedeR)
library(igraph)
rdp <- RedPort() 
calld(rdp)

data(ER.deg)
dt <- ER.deg$dat
sg <- ER.deg$ceg
sg <- att.mapv(sg, dat=dt, refcol=1)

ADJ = get.adjacency(sg)
ADJ = as.matrix(ADJ)
plot_gray_image(ADJ)

sg <- att.setv(sg, from="Symbol", to="nodeAlias")
sg <- att.setv(sg, from="logFC.t3", to="nodeColor", breaks=seq(-1,1,0.2), pal=2)    
sg <- att.setv(sg, from="ERbdist", to="nodeSize", nquant=10, isrev=TRUE, xlim=c(5,40,1))

addGraph(rdp,sg)
hc <- hclust(dist(get.adjacency(sg, attr="weight")))
nestID <- nesthc(rdp,hc, cutlevel=3, nmemb=5, cex=0.3, labels=V(sg)$nodeAlias)
mergeOutEdges(rdp,nlev=2)

clusters= cutree(hc,11)
plot(hc)
rect.hclust(tree = hc,k = 11)
#relax(rdp)

scl <- sg$legNodeColor$scale
leg <- sg$legNodeColor$legend
addLegend.color(rdp, colvec=scl, labvec=leg, title="diff. gene expression (logFC)")

scl <- sg$legNodeSize$scale
leg <- sg$legNodeSize$legend
addLegend.size(rdp, sizevec=scl, labvec=leg, title="bd site distance (kb)")

nOBJ = nrow(ADJ)
nGenes = 1000
nMirna = 500

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
plot_gray_image(result_chromo_hbm_3c)

resultClusterGer <- cluster_gerarchico(listMatrix = list(result_chromo_hbm_3c), 
                                       methodHClust = "ward", nCluster = length(unique(clusters)), 
                                       feature = colnames(GE), 
                                       etichette= clusters )
p = 1-resultClusterGer[[1]]$errore
purity1 = c(purity1,p)

gg = graph.adjacency(adjmatrix =result_chromo_hbm_3c,mode = "undirected",weighted = TRUE)
gg = simplify(gg)
VNE = vonNewmanEntropy(gg)$VNEntropy

x = as.matrix(t(table(result_chromo_hbm_3c,ADJ)))
x[1,] = c(296,4188,23768)
x[2,] = c(2015, 3,14)
x = x/sum(x)
barplot(x)
legend(x = "topleft",legend = c("True edge","False edge"),fill = c("gray","black"))

