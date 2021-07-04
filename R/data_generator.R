#Function for generating the synthetic data
data_generator<-function(n = 100, p = 50, densy_degree = 0.05, sd_ls = 0.1, mean_ls = log(10),value_nondiagonal = 0.3){
  ##
  ## Step 1: generate the graph and the corresponding precision matrix
  G_1<-matrix(rep(0,p^2),nrow = p)
  G_1[upper.tri(G_1)][sample(1:(p*(p-1)/2),floor(densy_degree*(p*(p-1)/2))+1,replace = FALSE)]<-1
  index_1<-G_1[upper.tri(G_1)]
  G_1<-t(G_1)+G_1
  diag(G_1)<-rep(1,p)
  bg_mat<-diag(p)
  bg_mat[upper.tri(bg_mat)][index_1==1]<-value_nondiagonal
  bg_mat<-bg_mat+t(bg_mat)-diag(p)
  #generate a precision matrix
  scalar_para1<-0.3
  scalar_para2<-0.1
  Pre_mat<-G_1*bg_mat+diag(rep(abs(min(eigen(G_1*bg_mat)$values))+scalar_para2,p))

  ##
  ## Step 2: generate the library size
  S_depth<-exp(rnorm(n,mean = mean_ls,sd = sd_ls))
  S_depth<-ifelse(S_depth==0,0.5,S_depth)

  ##
  ## Step 3: generator the observed count data
  log_a<-MASS::mvrnorm(n=n, c(rep(-2,floor(p*0.6)),rep(-0.5,(p - floor(p*0.6)))), solve(Pre_mat))
  log_a<-ifelse(log_a>50,50,log_a)
  a<-exp(log_a)
  b<-a * matrix(rep(S_depth,each=p),nrow = n,byrow = TRUE)
  Y<-matrix(NA,nrow = n,ncol = p)
  for(i in 1:n){
    for(j in 1:p){
      Y[i,j]<-rpois(1,b[i,j])
    }
  }

  ##names
  colnames(Y)<-paste("Gene",1:ncol(Y),sep = "")
  rownames(Y)<-paste("Cell",1:nrow(Y),sep = "")
  rownames(Pre_mat)<-paste("Gene",1:ncol(Y),sep = "")
  colnames(Pre_mat)<-paste("Gene",1:ncol(Y),sep = "")
  names(S_depth)<-paste("Cell",1:nrow(Y),sep = "")
  ##
  ## Return the result
  return(list(obs_mat = Y, S_depth = S_depth, Pre_mat = Pre_mat))
  ##
}
