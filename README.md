# DtracePLN: Dtrace method for poisson log-normal model

> DtracePLN is a package for calculating the precision matrix of poission log-normal model using moment method and dtrace loss. We first calculate the covariance matrix by moment method and then calculate the precision matrix by dtrace loss using R package EQUAL.

## Installation
**DtracePLN** is available on [Github](https://github.com/git-tjj/DtracePLN).

### R Package installation
- Before installing **DtracePLN** package, please install **cescwang85-EQUAL-da9cbbb.tar.gz** on [Github](https://github.com/git-tjj/DtracePLN).
- For the development version, use the github install
```{r package github, eval = FALSE}
remotes::install_github("https://github.com/git-tjj/DtracePLN")
```

## Usage and main functions
The package comes with a single-cell RNA sequence data about Interferon β-Stimulated PBMCs to present the functionality of main function.

### DtracePLN
```{r load DtracePLN, eval = FALSE}
library(DtracePLN)

## simulation data for testing
data("obs_mat")
DtracePLN_res_sim<-DtracePLN::DtracePLN(obs_mat = DtracePLN::obs_mat,
Sd_est = "GMPR",lambda_vec = NULL , n_lambda = 50 , lambda_valuemax0 = 10)

## real data analysis
data("realdatamatrix")
DtracePLN_res_real<-DtracePLN::DtracePLN(obs_mat = DtracePLN::realdatamatrix,
Sd_est = "GMPR",lambda_vec = NULL , n_lambda = 50 , lambda_valuemax0 = 10)
```

### data_generator
```{r, warning = FALSE}
library(MASS)
data_generator_res<-DtracePLN::data_generator(n = 100, p = 50, densy_degree = 0.05,
sd_ls = 0.1, mean_ls = log(10),value_nondiagonal = 0.3)
```


## Reference

Please cite our work using the following references:
