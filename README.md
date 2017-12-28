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

Example Usage
-------------

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

