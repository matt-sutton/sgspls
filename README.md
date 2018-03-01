Sparse group-subgroup PLS
=========================

`sgspls` is an R package that provides an implementation of an alternating convex search algorithm for computing partial least squares solutions with multiple grouping levels. Most of the linear algebra is written in C++ using [Rcpp](http://www.rcpp.org) and the [Armadillo](http://arma.sourceforge.net) C++ library. 

A preliminary paper has been submitted describing the statistical properties of the sparse group-subgroup PLS with grouping structure and will be made avaliable when the manuscript is accepted. 

Installation
------------

Use [devtools](https://github.com/hadley/devtools) to install:

```R
library(devtools)
install_github("sgspls", "matt-sutton")
```

Example Usage (regression)
--------------------------

```R
 library(sgspls)

 set.seed(1)
 n = 50; p = 500; 
 size.groups = 30; size.subgroups = 5
 groupX <- ceiling(1:p / size.groups)
 subgroupX <- ceiling(1:p / size.subgroups)
 
 X = matrix(rnorm(n * p), ncol = p, nrow = n)
 
 beta <- rep(0,p)
 bSG <- -2:2; b0 <- rep(0,length(bSG))
 betaG <- c(bSG, b0, bSG, b0, bSG, b0)
 beta[1:size.groups] <- betaG
 
 y = X %*% beta + 0.1*rnorm(n)
 
 #--------------------------------------#
 #-- Set up a basic model to tune --#
 
 parameters <- list(X=X, Y=y, groupX=groupX, subgroupX=subgroupX)
 
 #---------------------------------------------#
 #-- Tune over 1 to 2 groups and multiple    --#
 #-- sparsity levels for the first component --#
 
 cv_pls_comp1 <- tune_sgspls(parameters = parameters, group_seq = 1:2, scale_resp = F)
 
 #-- MSEP is on the original scale for the response --#
 cv_pls_comp1$results_tuning
 cv_pls_comp1$best
 
 # Use the optimal fit for the first component and tune the second component
 cv_pls_comp2 <- tune_sgspls(parameters = cv_pls_comp1$parameters, group_seq = 1:2, scale_resp = F)
 cv_pls_comp2
 
 # Use the optimal fit for the second component and tune the third component
 cv_pls_comp3 <- tune_sgspls(parameters = cv_pls_comp2$parameters, group_seq = 1:2, scale_resp = F)
 
 #--------------------------------------#
 #------ Inspect tuning function  ------#
 
 # Look at the plot of validation curves
 plot(cv_pls_comp3)
 
 # get the optimal values
 cv_pls_comp3
 
 # Fit the pls model with these parameters
 model <- do.call(sgspls, args = cv_pls_comp3$parameters)
 
 # See model details
 model
 
 # get the regression coefficients
 model_coef <- coef(model, type = "coefficients", comps = 3)
 
 cbind(beta, model_coef$B)

```

Example Usage (canonical)
--------------------------

```R
set.seed(1)
n = 50; p = 500; q = 100
ncomp <- 3

Scores <- MASS::mvrnorm(n, rep(0, ncomp), Sigma = diag(1, ncomp, ncomp))
Tm <- Scores
Sm <- Scores + MASS::mvrnorm(n, rep(0, ncomp), Sigma = diag(0.1, ncomp, ncomp))


subgroupC <- -2:2; p_sg <- length(subgroupC) 
groupC <- c(subgroupC, rep(0,p_sg), subgroupC, rep(0,p_sg), subgroupC, rep(0,p_sg))
p_g <- length(groupC) 

groupX <- ceiling(1:p / p_g)
subgroupX <- ceiling(1:p / p_sg)


Cm <- matrix(0,nrow = p, ncol = ncomp)
Cm[groupX == 1, 1] <- groupC
Cm[groupX == 2, 2] <- groupC
Cm[groupX == 3, 3] <- groupC
Dm <- matrix(runif(q*ncomp), ncol = ncomp)

X <- Tm%*%t(Cm) + MASS::mvrnorm(n, rep(0, p), Sigma = diag(0.1,p,p))
Y <- Sm%*%t(Dm) + MASS::mvrnorm(n, rep(0, q), Sigma = diag(0.1,q,q))

model <- sgspls(X, Y, ncomp = 5, mode = "regression", keepX = 17,
                groupX = groupX, subgroupX = subgroupX,
                indiv_sparsity_x = 0.8, subgroup_sparsity_x = 0.15)

#-- Scree type plot for the sgspls canonical method --#
plot(model)

#-- Sparse Version --#
model <- sgspls(X, Y, ncomp = 3, mode = "regression", keepX = 1,
                groupX = groupX, subgroupX = subgroupX,
                indiv_sparsity_x = 0.8, subgroup_sparsity_x = 0.15)

# check the recovery
model$weights$X

```
