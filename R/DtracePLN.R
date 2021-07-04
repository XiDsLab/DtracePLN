
#main function DtracePLN
DtracePLN<-function(obs_mat, Sd_est = "TSS", lambda_vec = NULL, n_lambda =100,lambda_valuemax0 = 10){

  ##names
  Cellname<-rownames(obs_mat)
  Genename<-colnames(obs_mat)
  ##1. Determine the S_depth--------------------
  if(length(Sd_est) == 1){
    S_depth<-compute_offset(obs_mat,offset = Sd_est)
  }else{
    S_depth <- Sd_est
  }
  message("1. library size determination has been completed.")
#################

  ##2. Estimate the covariance matrix by moment method------------------------
  time_sigmamoment1<-Sys.time()
  cov_Y1<-sigma_momest(data_use = obs_mat,S_depth = S_depth)
  time_sigmamoment2<-Sys.time()
  time_sigmamoment<-as.double(difftime(time_sigmamoment2,time_sigmamoment1,units = "hours"))
  message("2. The estimation of the covariance matrix by moment method has been completed.")
###################################

  ##3. Determine the lambda_vec--------------------------
  eigres<-eigen(cov_Y1)
  eigvalue_vec<-ifelse(eigres$values<0,0,eigres$values)
  obs_mat_s<-t(eigres$vectors%*%diag(sqrt(eigvalue_vec))*sqrt(ncol(obs_mat)-1))
  rm(eigres);

  if(is.null(lambda_vec)){
    #find the suitable lambda_max and lambda_min
    lambda_valuemax<-lambda_valuemax0
    con_1<-TRUE
    while (con_1) {
      aa<-as.matrix((EQUAL(obs_mat_s,lambda = lambda_valuemax))$Omega[[1]])
      nonzero_num <- length(which(aa[upper.tri(aa)]!=0))
      if(nonzero_num==0){
        lambda_valuemax<-lambda_valuemax*0.8
      }else{
        lambda_valuemax<-lambda_valuemax * 1.2
        con_1<-FALSE
      }
      rm(aa);
    }
    lambda_valuemin<-lambda_valuemax
    con_2<-TRUE
    while (con_2) {
      aa<-as.matrix((EQUAL(obs_mat_s,lambda = lambda_valuemin))$Omega[[1]])
      zero_num <- length(which(aa[upper.tri(aa)]==0))
      if(zero_num!=0){
        lambda_valuemin<-lambda_valuemin*0.8
      }else{
        con_2<-FALSE
      }
      rm(aa)
    }
    # lambda_vec <- exp(seq(from = log(lambda_valuemin),to=log(lambda_valuemax),length.out = n_lambda))
    lambda_vec <- seq(from = lambda_valuemin,to=lambda_valuemax,length.out = n_lambda)
  }
  message("3. Lambda_vec determination has been completed.")
######################################

  ##4. Run EQUAL for sparse precision matrix estimation-----------------------
  time_EQUAL1<-Sys.time()
  EQUAL_res<-EQUAL(obs_mat_s,lambda = lambda_vec)
  time_EQUAL2<-Sys.time()
  time_EQUAL<-as.double(difftime(time_EQUAL2,time_EQUAL1,units = "hours"))
  rm(obs_mat_s)
  lambda_vec<-EQUAL_res$lambda
  message("4. The estimation of precision matrix has been completed.")
##########################################

  ##5. Hyper-parameter selection by BIC criterion and AIC criterion----------------------
  ##5.1 BIC criterion
  BIC_vec<-c()
  for (l in 1:length(lambda_vec)) {
    BIC_value<-nrow(obs_mat) * sqrt(sum(((EQUAL_res$Omega[[l]] %*% cov_Y1 + cov_Y1 %*% EQUAL_res$Omega[[l]])/2 - diag(nrow(cov_Y1)))^2)) + log(nrow(obs_mat)) *
      (sum(ifelse(EQUAL_res$Omega[[l]]==0,0,1)))/2
    BIC_vec<-c(BIC_vec,BIC_value)
  }
  lambda_chooseB<-lambda_vec[which.min(BIC_vec)]
  Omega_chooseB<-EQUAL_res$Omega[[which.min(BIC_vec)]]

  ##5.2 AIC criterion
  AIC_vec<-c()
  for (l in 1:length(lambda_vec)) {
    AIC_value<-nrow(obs_mat) * sqrt(sum(((EQUAL_res$Omega[[l]] %*% cov_Y1 + cov_Y1 %*% EQUAL_res$Omega[[l]])/2 - diag(nrow(cov_Y1)))^2)) + 2 *
      (sum(ifelse(EQUAL_res$Omega[[l]]==0,0,1)))/2
    AIC_vec<-c(AIC_vec,AIC_value)
  }
  lambda_chooseA<-lambda_vec[which.min(AIC_vec)]
  Omega_chooseA<-EQUAL_res$Omega[[which.min(AIC_vec)]]
  message("5. Hyper-parameter selection by BIC criterion and AIC criterion have been completed.")
############################################

  ##names
  for(l in 1:length(EQUAL_res$lambda)){
    rownames(EQUAL_res$Omega[[l]])<-Genename
    colnames(EQUAL_res$Omega[[l]])<-Genename
  }
  names(S_depth)<-Cellname
  rownames(Omega_chooseB)<-Genename
  colnames(Omega_chooseB)<-Genename
  rownames(Omega_chooseA)<-Genename
  colnames(Omega_chooseA)<-Genename

  # rownames(Omega_est)<-Omega_est = EQUAL_res$Omega
  DtracePLN_res<-list(Omega_est = EQUAL_res$Omega, lambda_vec = EQUAL_res$lambda, S_depth = S_depth,
                      BIC_vec = BIC_vec, Omega_chooseB = Omega_chooseB,
                      AIC_vec = AIC_vec, Omega_chooseA = Omega_chooseA,
                      time_sigmamoment = time_sigmamoment,time_EQUAL = time_EQUAL)
  return(DtracePLN_res)
}

