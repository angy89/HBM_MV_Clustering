createHBM3C<-function(listHBM){
  listHBMControllo <- c()
  listHBMControllo[[1]]<- listHBM[[1]]
  listHBMControllo[[2]]<- listHBM[[2]]
  for(i in 2:length(listHBM)){
    if(i==2){
      HBM_3C<-chromoHBM3C(listHBM[[1]]$hbm,listHBM[[2]]$hbm)
    }
    else{
      HBM_3C<-chromoHBM3C(HBM_3C,listHBM[[i]]$hbm)
    }
  }
  return(HBM_3C)
}