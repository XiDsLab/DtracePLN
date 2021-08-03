\name{PLNet}
\alias{PLNet}
\title{Estimate the sparse precision matrix for PLN model.}
\usage{
PLNet(
  obs_mat,
  Sd_est = "TSS",
  lambda_vec = NULL,
  n_lambda = 100,
  lambda_valuemax0 = 10,
  penalize.diagonal = TRUE
)
}
\description{
Estimate the precision matrix for PLN model by estimating the covariance matrix using moment method and estimating the precision matrix using the dtrace loss.
}

\arguments{

\item{obs_mat}{An abundance count matrix with Cell rows and Gene columns.}

\item{Sd_est}{Optional. Normalisation scheme used to compute scaling factors used as offset during PLN inference. Available schemes are "TSS" (Total Sum Scaling, default), "CSS" (Cumulative Sum Scaling, used in metagenomeSeq), "GMPR" (Geometric Mean of Pairwise Ratio, introduced in Chen et al., 2018) or "none". Alternatively the user can supply its own vector or matrix of offsets (see note for specification of the user-supplied offsets).}

\item{lambda_vec}{Optional. The user can supply its own vector.}

\item{n_lambda}{Optional. When the lambda_vec is empty, the package will provide a suitable lambda_vec with length n_lambda.}

\item{lambda_valuemax0}{Optional. When the lambda_vec is NULL, the package will provide a suitable lambda_vec searching from lambda_valuemax0.}

\item{penalize.diagonal}{Should diagonal of inverse covariance be penalized? Dafault TRUE.}

}

\value{

A list which contain the following result:

\item{Omega_est}{A list which contain the precision matrix for different lambda value.}

\item{lambda_vec}{The lambda_vec used.}

\item{S_depth}{The vector of estimated library size.}

\item{BIC_vec}{The vector of BIC value.}

\item{Omega_chooseB}{The estimated precision matrix choosed by the BIC criterion.}

\item{AIC_vec}{The vector of AIC value.}

\item{Omega_chooseA}{The estimated precision matrix choosed by the AIC criterion.}

\item{time_sigmamoment}{The time for calculating the covariance matrix for PLN using moment method.}

\item{time_EQUAL}{The time for calculating the path of precision matrix using EQUAL package.}

}

\references{
Zhang T, Zou H. Sparse precision matrix estimation via lasso penalized D-trace loss[J]. Biometrika, 2014, 101(1): 103-120.
}
\examples{

##
library(MASS)
library(CVXR)
library(EQUAL)

##
# Run PLNet
PLNet_res<-PLNet::PLNet(obs_mat = PLNet::obs_mat,
Sd_est = "TSS",lambda_vec = NULL , n_lambda =100 , lambda_valuemax0 = 10,penalize.diagonal = TRUE)


}
