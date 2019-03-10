nestedBlockMatrix = function(nOutBlocks=2,blockParams,minOutNoise,maxOutNoise,main=""){
  
  
  BML = list()#block matrices list
  groupCodes = c()
  colorCodes = c()
  for(i in 1:nOutBlocks){
    nBlocks = blockParams[[i]]$nBlocks
    nOBJ = blockParams[[i]]$nOBJ
    minCorr = blockParams[[i]]$minCorr
    maxCorr = blockParams[[i]]$maxCorr
    minNoise = blockParams[[i]]$minNoise
    maxNoise = blockParams[[i]]$maxNoise
    
    
    nn = rand_vect(nBlocks,nOBJ)
    print("Generating first block matrix...\n")
    res1 = genBlockMatrix(nn,minCorr=minCorr, maxCorr=maxCorr,minNoise=minNoise,maxNoise = maxNoise)
    if(i ==1){
      groupCodes = c(groupCodes,res1$groupCodes)
      
    }else{
      groupCodes = c(groupCodes,res1$groupCodes + max(groupCodes))
      
    }
    colorCodes = c(colorCodes,res1$colorCodes)
    BML[[i]] = res1$blockMatrix
  }
 
  
  blockMatrix<-blockMatrixDiagonal(BML,minNoise=minOutNoise,maxNoise = maxOutNoise,external=TRUE)#m1,m2,m3,m4)
  
  hls = hclust(as.dist(1 - blockMatrix),method = "average")
#  plot_dend_with_class_color(hls,groupCodes,colorCodes)
 # plot_gray_image(blockMatrix,main = main)
  
  return(list(blockMatrix=blockMatrix,hls=hls,groupCodes=groupCodes,colorCodes=colorCodes))
}


nestedBlockMatrix3 = function(nExternalBlock = 2,nOutBlocks=2,blockParams,minOutNoise,maxOutNoise){

  BML = list()#block matrices list
  groupCodes = c()
  colorCodes = c()
  for(i in 1:nOutBlocks){
    nBlocks = blockParams[[i]]$nBlocks
    nOBJ = blockParams[[i]]$nOBJ
    minCorr = blockParams[[i]]$minCorr
    maxCorr = blockParams[[i]]$maxCorr
    minNoise = blockParams[[i]]$minNoise
    maxNoise = blockParams[[i]]$maxNoise
    
    
    nn = rand_vect(nBlocks,nOBJ)
    print("Generating first block matrix...\n")
    res1 = genBlockMatrix(nn,minCorr=minCorr, maxCorr=maxCorr,minNoise=minNoise,maxNoise = maxNoise)
    groupCodes = c(groupCodes,res1$groupCodes)
    colorCodes = c(colorCodes,res1$colorCodes)
    BML[[i]] = res1$blockMatrix
  }
  
  
  blockMatrix<-blockMatrixDiagonal(BML,minNoise=minOutNoise,maxNoise = maxOutNoise)#m1,m2,m3,m4)
  
  hls = hclust(as.dist(1 - blockMatrix),method = "average")
  plot_dend_with_class_color(hls,groupCodes,colorCodes)
  plot_gray_image(blockMatrix,main = main)
  
  return(list(blockMatrix=blockMatrix,hls=hls,groupCodes=groupCodes,colorCodes=colorCodes))
}