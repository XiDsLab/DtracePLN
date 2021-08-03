\name{data_generator}
\alias{data_generator}

\title{
Syntheic data generator
}
\description{
Generate the syntheic data from multivariate poisson log-normal distribution.
}
\usage{
data_generator(n = 100, p = 50, densy_degree = 0.05, sd_ls = 0.1, mean_ls = log(10), value_nondiagonal = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{Sample size.}
  \item{p}{Dimension.}
  \item{densy_degree}{The proportion of edge connected in network.}
  \item{sd_ls}{The standard deviation of library size.}
  \item{mean_ls}{The expectation of library size.}
  \item{value_nondiagonal}{The valus of non-diagonal element in precision matrix.}
}
\value{
A list which contain the following result:

\item{obs_mat}{The matrix of syntheic count data.}

\item{S_depth}{The vector of syntheic library size.}

\item{Pre_mat}{The matrix of syntheic precision matrix.}

}


\examples{
# Library package
library(MASS)

##
# Generator the syntheic data
data_generator_res<-PLNet::data_generator(n = 100, p = 50, densy_degree = 0.05,
sd_ls = 0.1, mean_ls = log(10),value_nondiagonal = 0.3)
}

