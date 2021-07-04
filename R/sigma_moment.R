#Function for calculate the covariance matrix by moment method
sigma_momest<-function(data_use,S_depth){
  ##1. Calculate the sample covariance matrix by moment method----------------
  alpha_1<-colMeans(data_use/S_depth)
  Y_tlide<-data_use/S_depth
  sigma_me<-log(((t(Y_tlide) %*% Y_tlide)/nrow(data_use))/(matrix(alpha_1,ncol = 1) %*% matrix(alpha_1,nrow = 1)))
  diag(sigma_me)<-log(colMeans((data_use^2-data_use)/matrix(rep(S_depth^2,each=ncol(data_use)),nrow = nrow(data_use),byrow = TRUE))/alpha_1^2)
  rm(alpha_1,Y_tlide)
  ###################################

  ##2. Handle the abnormal situation of sample covariance matrix------------------
  ##2.1 non-diagnoal
  isinfinite_mat<-ifelse(is.infinite(sigma_me),1,0)
  isinfinite_mat[lower.tri(isinfinite_mat)]<-0
  min_vec<-rep(NA,dim(data_use)[2])
  for(i in 1:dim(data_use)[2]){
    min_vec[i] <- min(sigma_me[i,-i][is.finite(c(sigma_me[i,-i]))])
  }
  for(i in 1:(dim(data_use)[2]-1)){
    index_select<-which(isinfinite_mat[i,]==1)
    index_select<-index_select[which(index_select>i)]
    if(length(index_select)>0){
      for(j in 1:length(index_select)){
        isinfinite_mat[i,index_select[j]]<-min(min_vec[i],min_vec[index_select[j]])
      }
    }
  }
  isinfinite_mat<-isinfinite_mat + t(isinfinite_mat)
  diag(isinfinite_mat)<-0
  sigma_me<-ifelse(is.infinite(sigma_me),0,sigma_me)
  sigma_me<-sigma_me + isinfinite_mat
  rm(isinfinite_mat)

  ##2.2 diagnoal
  error_index<-which(!is.finite(as.vector(diag(sigma_me))))
  if(length(error_index)>0){
    nonerror_index<-setdiff(1:(dim(data_use)[2]),error_index)
    for(i in 1:length(error_index)){
      diag(sigma_me)[error_index[i]]<-max((sigma_me[error_index[i],nonerror_index])^2/(diag(sigma_me))[nonerror_index])
    }
  }
################################

  ##3. Optimatization by minimize the finite norm--------------------
  S<-Variable(dim(sigma_me)[1], dim(sigma_me)[1], PSD = TRUE)

  obj<-Minimize(max(abs(S-sigma_me)))

  prob <- CVXR::Problem(obj)

  result <- CVXR::solve(prob,solver="SCS",verbose=FALSE)
  sigma_me1 <- result$getValue(S)
  sigma_me1<-max(abs(sigma_me1-sigma_me))* diag(dim(data_use)[2]) + sigma_me1
  rm(sigma_me,S);
########################################
  return(sigma_me1)
}
