#Construct the subdiagnoal matrix
upper_diag<-function(p,diag_num){
  AAA<-matrix(0,nrow = p,ncol = p)
  for(j in 1:(p-1)){
    AAA[j,]<-(c(rep(0,j),rep(1,diag_num),rep(0,ifelse(p-j-diag_num>0,p-j-diag_num,0))))[1:p]
  }
  AAA_1<-AAA+t(AAA)+diag(p)
  return(AAA_1)
}
