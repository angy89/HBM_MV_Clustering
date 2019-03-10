rmse_mat = function(x,y){
  sqrt(sum(abs(as.vector(x) - as.vector(y))^2))
}
