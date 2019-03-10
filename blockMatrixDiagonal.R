# builds a block matrix whose diagonals are the square matrices provided.
# m1=matrix(runif(10*10),nrow=10,ncol=10)
# m2=matrix(runif(5*5),nrow=5,ncol=5)
# blockMatrix<-blockMatrixDiagonal(m1,m2,m2,m1)
# or
# blockMatrix<-blockMatrixDiagonal(list(m1,m2,m2,m1))
# C.Ladroue

blockMatrixDiagonal<-function(matrixList,minNoise=0,maxNoise = 0.7,external=FALSE,...){  
  # matrixList<-list(...)
  # if(is.list(matrixList[[1]])) matrixList<-matrixList[[1]]
  
  dimensions<-sapply(matrixList,FUN=function(x) dim(x)[1])
  finalDimension<-sum(dimensions)
  finalMatrix<-matrix(0,nrow=finalDimension,ncol=finalDimension)
  index<-1
  for(k in 1:length(dimensions)){
    finalMatrix[index:(index+dimensions[k]-1),index:(index+dimensions[k]-1)]<-matrixList[[k]]
    index<-index+dimensions[k]
  }
  
  for(i in 1:(nrow(finalMatrix)-1)){
    for(j in (i+1):nrow(finalMatrix)){
      if(finalMatrix[i,j]==0){
        
        if(external){
          x = runif(n = 1,min=minNoise,max = maxNoise)
          finalMatrix[i,j] = finalMatrix[j,i] =  sample(c(0,1),size=1,prob=c(0.95,0.05))
        }else{
          x = runif(n = 1,min=minNoise,max = maxNoise)
          finalMatrix[i,j] = finalMatrix[j,i] = sample(c(0,x),size=1,prob=c(0.8,0.2))
        }
      }
    }
  }
  
  finalMatrix
}

#Generate an nxn Correlation matrix
myCorrMat = function(n,minU,maxU){
  R <- matrix(runif(n^2, min = minU,max= maxU), ncol=n) 
  RtR <- R %*% t(R) 
  Q <- cov2cor(RtR) 
}

plot_dend_with_class_color = function(hls,groupCodes,colorCodes){
  require(dendextend)
  
  dend <- as.dendrogram(hls)
  # Assigning the labels of dendrogram object with new colors:
  labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]
  # Plotting the new dendrogram
  plot(dend)
}


#Generate a block matrix with four blocks
genBlockMatrix = function(nList,main,minCorr=0.8, maxCorr=1,minNoise=0,maxNoise = 0.7){
  
  # if(is.null(n1)) n1 = ceiling(runif(n = 1,min = 1,max = 50))
  # if(is.null(n2)) n2 = ceiling(runif(n = 1,min = 1,max = 50))
  # if(is.null(n3)) n3 = ceiling(runif(n = 1,min = 1,max = 50))
  # if(is.null(n4)) n4 = ceiling(runif(n = 1,min = 1,max = 50))
  
  colors = rainbow(length(nList))
  mList = list()
  groupCodes = c()
  colorCodes = c()
  for(i in 1:length(nList)){
    mList[[i]] = myCorrMat(nList[i],minCorr,maxCorr)
    #groupCodes = c(groupCodes,rep(paste("M",i,sep=""),nList[i]))
    groupCodes = c(groupCodes,rep(i,nList[i]))
    
    colorCodes = c(colorCodes, colors[i])             
  }
  names(colorCodes) = unique(groupCodes)
  
  blockMatrix<-blockMatrixDiagonal(mList,minNoise=0,maxNoise = 0.7)#m1,m2,m3,m4)
  
  hls = hclust(as.dist(1 - blockMatrix),method = "average")
 # plot_dend_with_class_color(hls,groupCodes,colorCodes)
#  plot_gray_image(blockMatrix,main = main)
  
  return(list(blockMatrix=blockMatrix,hls=hls,groupCodes=groupCodes,colorCodes=colorCodes))
}

#Generate N number that sum to M
rand_vect <- function(N, M, sd = 1, pos.only = TRUE) {
  vec <- rnorm(N, M/N, sd)
  if (abs(sum(vec)) < 0.01) vec <- vec + 1
  vec <- round(vec / sum(vec) * M)
  deviation <- M - sum(vec)
  for (. in seq_len(abs(deviation))) {
    vec[i] <- vec[i <- sample(N, 1)] + sign(deviation)
  }
  if (pos.only) while (any(vec < 0)) {
    negs <- vec < 0
    pos  <- vec > 0
    vec[negs][i] <- vec[negs][i <- sample(sum(negs), 1)] + 1
    vec[pos][i]  <- vec[pos ][i <- sample(sum(pos ), 1)] - 1
  }
  vec
}


plot_gray_image = function(BM,main){
  library(fields)
 # diag(BM) = 0
  colourScale<-seq(0,0.9,length.out=100)
  image.plot(BM,asp=1,col=rgb(colourScale,colourScale,colourScale),ann=FALSE,xaxt="n",yaxt="n",axes=FALSE)
  #legend.col(col = colourScale, lev = colourScale)

  #library(lattice)
  #levelplot(BM,colourScale)
  
}
