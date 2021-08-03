\name{sigma_momest}
\alias{sigma_momest}
\title{Estimate the covariance matrix using moment method for PLN model}
\usage{
sigma_momest(
  data_use,
  S_depth
)
}
\description{
Estimate the covariance matrix using moment method for PLN model.
}

\arguments{

\item{data_use}{The matrix of count data.}

\item{S_depth}{A Vector, which contains the library size estimation.}

}

\value{

The estimated covariance matrix for PLN model by using moment method.

}

\examples{
# Library package
library(MASS)
library(CVXR)
library(EQUAL)

##
# Generator the syntheic data
data_generator_res<-PLNet::data_generator(n = 100, p = 50, densy_degree = 0.05,
sd_ls = 0.1, mean_ls = log(10),value_nondiagonal = 0.3)

##
S_depth<-PLNet::compute_offset(counts = data_generator_res$obs_mat,offset = "TSS")
cov_Y1<-PLNet::sigma_momest(data_use = data_generator_res$obs_mat,S_depth = S_depth)

}
