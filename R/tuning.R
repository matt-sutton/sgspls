#' Calculate percentage of variance explained (PVE)
#'
#' Function to evaluate the PEV and recommend a number of components to select.
#'
#' @param object Object of class inheriting from \code{"sgspls"}. 
#' The function will retrieve some key parameters stored in that object.
#' 
#' @export
#' @return \code{PVE} returns a vector of percentage of variance explained by the scores in the PLS method.
#' See references for details.
#' 
#' @seealso Tuning functions \code{\link{calc_pve}},
#' \code{\link{tune_sgspls}}, \code{\link{tune_groups}}. 
#' Model performance and estimation  \code{\link{predict}}, \code{\link{perf.sgspls}}, \code{\link{coeff.sgspls}} 
#' @references Liquet Benoit, Lafaye de Micheaux, Boris Hejblum, Rodolphe
#'   Thiebaut. A group and Sparse Group Partial Least Square approach applied in
#'   Genomics context. \emph{Submitted}.
#' 
#' @examples
#'
#' set.seed(1)
#' n = 50; p = 500; 
#' 
#' X = matrix(rnorm(n * p), ncol = p, nrow = n)
#' 
#' beta <- rep(0,p)
#' beta[1:20] <- 1:20
#' 
#' y = X %*% beta + 0.1*rnorm(n)
#' 
#' model <- sgspls(X, y, ncomp = 3, mode = "regression")
#' calc_pve(model)
#' 

calc_pve <- function(object){
  
  #--- Take the predictions from the PLS object --#
  Ypred <- predict(object, X)
  
  PVE <- matrix(
    1 - apply(Ypred, 3, function(Xi) sum(apply(Y - Xi,2,var))/sum(apply(Y,2,var))), 
    ncol = object$parameters$ncomp 
    )
  
  #--- Get the pve ---#
  if (max(PVE) > 0.98){
    ncomp <- which(PVE >0.95)[1]
  } else{
    ncomp <- which.min(PVE)[1]
  }
  cat("Recommend to use", ncomp, "components.\n\n")
  colnames(PVE) <- paste0("Comp ", 1:length(PVE))
  rownames(PVE) <- "PVE"
  
  return(PVE)
}



generate_sparsities <- function(){
  
  #--- Original code to generate the sparsities ---# 
  
  # # Choose the sparsities for individual level
  # sparsities = seq(0.01,0.99,0.01)
  # grp = c(seq(0.05,0.95,0.1))
  # 
  # # third column to ensure individual sparsity
  # sparsities.s <- do.call(expand.grid,list(grp,sparsities,sparsities)) 
  # sparsities <- NULL
  # ind <- which(sparsities.s[,1] %in% grp )
  # sparsities.s <- sparsities.s[ind,]
  # for(sp in grp){
  #
  #   #-- take only the valid sparsities --#
  #   ind <- which(apply(sparsities.s,1,sum) == 1 & sparsities.s[,1] == sp )
  #
  #   #-- take the last, first and middle ones --#
  #   ind <- c(tail(ind,1),ind[floor(length(ind)/2)],head(ind,1))
  #
  #   sparsities <- rbind(sparsities,sparsities.s[ind,])
  # }
  # sparsities <- as.matrix(unique(sparsities))
  # rownames(sparsities) <- NULL
  # colnames(sparsities) <- c("Group", "Subgroup","Individual")
  
  return(
    structure(c(0.05, 0.05, 0.05, 0.15, 0.15, 0.15, 0.25, 0.25, 0.25, 
                0.35, 0.35, 0.35, 0.45, 0.45, 0.45, 0.55, 0.55, 0.55, 0.75, 0.75, 
                0.75, 0.85, 0.85, 0.85, 0.95, 0.95, 0.95, 0.01, 0.48, 0.94, 0.01, 
                0.43, 0.84, 0.01, 0.38, 0.74, 0.01, 0.33, 0.64, 0.01, 0.28, 0.54, 
                0.01, 0.23, 0.44, 0.03, 0.13, 0.22, 0.01, 0.08, 0.14, 0.01, 0.03, 
                0.04, 0.94, 0.47, 0.01, 0.84, 0.42, 0.01, 0.74, 0.37, 0.01, 0.64, 
                0.32, 0.01, 0.54, 0.27, 0.01, 0.44, 0.22, 0.01, 0.22, 0.12, 0.03, 
                0.14, 0.07, 0.01, 0.04, 0.02, 0.01), .Dim = c(27L, 3L), .Dimnames = list(
                  NULL, c("Group", "Subgroup", "Individual")))
  )
  
}


#' Compute cross-validated mean squared prediction error for sgspls regression
#'
#' Tuning function for finding the number of groups, and sparsities to select for an sgspls object.
#' Offers a sequential way to find optimal sparsities, and number of groups for either block.
#' 
#' @param parameters List of parameters to use in the current PLS object (see examples below).
#' @param sparsities Matrix of sparsities, with columns corresponding to group, subgroup and individual 
#' sparsity levels to tune over. If it is NULL then a preselected set of sparsity levels is used.
#' @param group_seq a vector containing the number of groups to tune over.
#' @param block A string either "X" or "Y" to indicate which block to tune parameters over.
#' @param folds The number of folds to use in cross validation.
#' @param progressBar Logical, indicating if a progress bar is shown.
#' @param setseed False, or integer for replicating tuning parameters.
#' @param scale_resp Logical, the MSEP is standardised across responses (see perf function for details).
#' 
#' @export
#' @return \code{tune_sgspls} returns a list with class \code{cv.sgspls} with measures:
#' \item{result_tuning}{A matrix containing the tuning parameters and MSEP values.}
#' \item{best}{A vector containing the optimal tuning parameters. }
#' \item{parameters}{A list of the parameters for a sgspls object.}
#' \item{tuning_sparsities}{A matrix of group, subgroup and individual sparsities tuned over.}
#' \item{folds}{Number of folds used in cross validation.}
#' \item{min_cv}{Minimum MSEP score.}
#' \item{group_seq}{Groups tuned over in cross validation.}
#' 
#' @seealso \code{\link{sgspls}} Tuning functions \code{\link{calc_pve}}, \code{\link{tune_groups}}. 
#' Model performance and estimation  \code{\link{predict}}, \code{\link{perf.sgspls}}, \code{\link{coeff.sgspls}} 
#' 
#' @references Liquet Benoit, Lafaye de Micheaux, Boris Hejblum, Rodolphe
#'   Thiebaut. A group and Sparse Group Partial Least Square approach applied in
#'   Genomics context. \emph{Submitted}.
#' 
#' @examples
#'
#'  set.seed(1)
#'  n = 50; p = 500; 
#'  size.groups = 30; size.subgroups = 5
#'  groupX <- ceiling(1:p / size.groups)
#'  subgroupX <- ceiling(1:p / size.subgroups)
#'  
#'  X = matrix(rnorm(n * p), ncol = p, nrow = n)
#'  
#'  beta <- rep(0,p)
#'  bSG <- -2:2; b0 <- rep(0,length(bSG))
#'  betaG <- c(bSG, b0, bSG, b0, bSG, b0)
#'  beta[1:size.groups] <- betaG
#'  
#'  y = X %*% beta + 0.1*rnorm(n)
#'  
#'  #--------------------------------------#
#'  #-- Set up a basic model to tune --#
#'  
#'  parameters <- list(X=X, Y=y, groupX=groupX, subgroupX=subgroupX)
#'  
#'  #---------------------------------------------#
#'  #-- Tune over 1 to 2 groups and multiple    --#
#'  #-- sparsity levels for the first component --#
#'  
#'  cv_pls_comp1 <- tune_sgspls(parameters = parameters, group_seq = 1:2, scale_resp = F)
#'  
#'  #-- MSEP is on the original scale for the response --#
#'  cv_pls_comp1$results_tuning
#'  cv_pls_comp1$best
#'  
#'  # Use the optimal fit for the first component and tune the second component
#'  cv_pls_comp2 <- tune_sgspls(parameters = cv_pls_comp1$parameters, group_seq = 1:2, scale_resp = F)
#'  cv_pls_comp2
#'  
#'  # Use the optimal fit for the second component and tune the third component
#'  cv_pls_comp3 <- tune_sgspls(parameters = cv_pls_comp2$parameters, group_seq = 1:2, scale_resp = F)
#'  
#'  #--------------------------------------#
#'  #------ Inspect tuning function  ------#
#'  
#'  # Look at the plot of validation curves
#'  plot(cv_pls_comp3)
#'  
#'  # get the optimal values
#'  cv_pls_comp3
#'  
#'  # Fit the pls model with these parameters
#'  model <- do.call(sgspls, args = cv_pls_comp3$parameters)
#'  
#'  # See model details
#'  model
#'  
#'  # get the regression coefficients
#'  model_coef <- coef(model, type = "coefficients", comps = 3)
#'  
#'  cbind(beta, model_coef$B)
#'  
#'  
tune_sgspls <-
  function(parameters, sparsities=NULL, group_seq=NULL, block = "X",
           folds= 10, progressBar=TRUE, setseed = 1, scale_resp = TRUE) {

    #----------------------#
    #--- Get sparsities ---#
    
    if(is.null(sparsities)){
      sparsities <- generate_sparsities()
    }
    
    #--- Progress Bar ---#
    if (progressBar == TRUE)
      pb <- txtProgressBar(style = 3)

    #--- Initalise ---#
    ngroups <- length(group_seq)
    nsparsities <- dim(sparsities)[1]
    results_tuning <- NULL
    
    if(is.null(parameters$ncomp)) {
      parameters$ncomp <- 1
    }
    
    if(!(class(parameters) %in% "cv.sgspls")){
      ncomp <- 1:parameters$ncomp
    } else {
      ncomp <- parameters$ncomp + 1 
    }
    
    #-------------------------------------------------#
    #--- Loop through sparsities combinations on j ---#
    
    for( j in 1:nsparsities){
      
      #-- Update the sparsity ---#
      
      if(block == "X" || block == "x")
      {
        parameters$indiv_sparsity_x[ncomp] <- sparsities[j,3]
        parameters$subgroup_sparsity_x[ncomp] <- sparsities[j,2]
      }
      
      if(block == "Y" || block == "y")
      {
        parameters$indiv_sparsity_y[ncomp] <- sparsities[j,3]
        parameters$subgroup_sparsity_y[ncomp] <- sparsities[j,2]
      }
      
      #--- Loop through groups ---#
      object_groups <- tune_groups(parameters, group_seq=group_seq, setseed = setseed,
                                   block = block, scale_resp = scale_resp)
      
      results <- cbind( group_seq, sparsities[j,3], sparsities[j,2], t(object_groups$MSEP[ncomp, ,drop=F]) )

      results_tuning <- rbind(results_tuning, results)
      if (progressBar == TRUE) {
        setTxtProgressBar(pb, j/(nsparsities))
      }
    }
    
    #------------------------------------#
    #-- create table of tuning results --#
    
    colnames(results_tuning) <- c("Groups", "Individual", "Subgroup", paste0("MSEP comp",ncomp))
    MSEPopt <- min(results_tuning[,-c(1:3)])
    best <- results_tuning[which(results_tuning == MSEPopt, arr.ind = T)[1],]
    keep <- best[1];
    indiv_sparsity <- best[2];
    subgroup_sparsity <- best[3]

    #-- Return parameters with optimal fit --#
    parameters <- object_groups$parameters
    parameters$ncomp <- max(ncomp)
    
    if(block == "X" || block == "x"){
      parameters$keepX[ncomp] <- keep
      parameters$indiv_sparsity_x[ncomp] <- indiv_sparsity
      parameters$subgroup_sparsity_x[ncomp] <- subgroup_sparsity
    }
    
    if(block == "Y" || block == "y"){
      parameters$keepY[ncomp] <- keep
      parameters$indiv_sparsity_y[ncomp] <- indiv_sparsity
      parameters$subgroup_sparsity_y[ncomp] <- subgroup_sparsity
    } 
    
    #-- alphas to match the paper --#
    rownames(results_tuning) <- paste("alpha",rep(1:length(sparsities[,1]), each = length(group_seq)),sep = "_")

    result <- list(results_tuning = results_tuning,
                   best = best,
                   parameters = parameters,
                   tuning_sparsities = sparsities,
                   folds = folds,
                   min_cv = min(results_tuning[,4]),
                   group_seq = group_seq,
                   folds = object_groups$folds)

    class(result) = c("cv.sgspls")
    return(result)
    }


#' Compute cross-validated mean squared prediction error for sgspls regression
#' 
#' Tuning function for finding the number of groups to select for an sgspls object.
#' Function is consistent in the sparsity parameters but may be tuned over the number of 
#' groups.
#'
#' @param parameters List of parameters to use in the current PLS object (see examples below).
#' @param group_seq a vector containing the number of groups to tune over.
#' @param block A string either "X" or "Y" to indicate which block to tune parameters over.
#' @param folds The number of folds to use in cross validation.
#' @param progressBar Logical, indicating if a progress bar is shown.
#' @param setseed False, or integer for replicating tuning parameters.
#' @param scale_resp Logical, the MSEP is standardised across responses (see perf function for details).
#' 
#' @export
#' @return \code{tune_groups} returns a list with class \code{cv.sgspls} with measures:
#' \item{result_tuning}{A matrix containing the tuning parameters and MSEP values.}
#' \item{best}{A vector containing the optimal tuning parameters. }
#' \item{parameters}{A list of the parameters for a sgspls object.}
#' \item{tuning_sparsities}{A matrix of group, subgroup and individual sparsities tuned over.}
#' \item{folds}{Number of folds used in cross validation.}
#' \item{min_cv}{Minimum MSEP score.}
#' \item{group_seq}{Groups tuned over in cross validation.}
#' 
#' @seealso \code{\link{sgspls}} Tuning functions \code{\link{calc_pve}}, \code{\link{tune_sgspls}}. 
#' Model performance and estimation  \code{\link{predict}}, \code{\link{perf.sgspls}}, \code{\link{coeff.sgspls}} 
#' @references Liquet Benoit, Lafaye de Micheaux, Boris Hejblum, Rodolphe
#'   Thiebaut. A group and Sparse Group Partial Least Square approach applied in
#'   Genomics context. \emph{Submitted}.
#' 
#' @examples
#'
#'  set.seed(1)
#'  n = 50; p = 500; 
#'  size.groups = 30; size.subgroups = 5
#'  groupX <- ceiling(1:p / size.groups)
#'  subgroupX <- ceiling(1:p / size.subgroups)
#'  
#'  X = matrix(rnorm(n * p), ncol = p, nrow = n)
#'  
#'  beta <- rep(0,p)
#'  bSG <- -2:2; b0 <- rep(0,length(bSG))
#'  betaG <- c(bSG, b0, bSG, b0, bSG, b0)
#'  beta[1:size.groups] <- betaG
#'  
#'  y = X %*% beta + 0.1*rnorm(n)
#'  
#'  #--------------------------------------#
#'  #-- Set up a basic model to tune --#
#'  
#'  parameters <- list(X=X, Y=y, groupX=groupX, subgroupX=subgroupX)
#'  
#'  #---------------------------------------------#
#'  #-- Tune over 1 to 2 groups and multiple    --#
#'  #-- sparsity levels for the first component --#
#'  
#'  cv_pls_comp1 <- tune_groups(parameters = parameters, group_seq = 1:2, scale_resp = F)
#'  
#'  #-- MSEP is on the original scale for the response --#
#'  cv_pls_comp1$MSEP
#'  
#'  # Use the optimal fit for the first component and tune the second component
#'  cv_pls_comp2 <- tune_groups(parameters = cv_pls_comp1$parameters, group_seq = 1:2, scale_resp = F)
#'  cv_pls_comp2$MSEP
#'  
#'  # Use the optimal fit for the second component and tune the third component
#'  cv_pls_comp3 <- tune_groups(parameters = cv_pls_comp2$parameters, group_seq = 1:2, scale_resp = F)
#'  cv_pls_comp3$MSEP
#'  
#'  #--------------------------------------#
#'  #------ Inspect tuning function  ------#
#'  
#'  # Fit the pls model with these parameters
#'  model <- do.call(sgspls, args = cv_pls_comp3$parameters)
#'  
#'  # See model details
#'  model
#'  
#'  # get the regression coefficients
#'  model_coef <- coef(model, type = "coefficients", comps = 3)
#'  
#'  cbind(beta, model_coef$B)
#'  
#'  
tune_groups <-
  function(parameters, group_seq = NULL, block = "X", folds = 10, scale_resp = TRUE, setseed = 1) {
    
    #-- Get the length of the groups to tune over ---#
    ngroups = length(group_seq); 
    
    if(is.null(group_seq))
      stop("Enter a sequence of groups to tune over")
    
    if(is.null(parameters$ncomp)) {
      parameters$ncomp <- 1
      }
    
    if(!(class(parameters) %in% "cv.sgspls")){
      ncomp <- 1:parameters$ncomp
    } else {
      ncomp <- parameters$ncomp + 1 
    }
    
    parameters$ncomp <- max(ncomp)
    
    MSEP <- matrix(0, ncol = ngroups, nrow = max(ncomp))
    
    #-- Tune over the group sequence --#
    for(i in 1:ngroups){
      
      #-- Update the group number ---#
      
      if(block == "X" || block == "x")
        parameters$keepX[ncomp] = group_seq[i]
      
      if(block == "Y" || block == "y")
        parameters$keepY[ncomp] = group_seq[i]
      
      
      object <- invisible(do.call(sgspls, args = parameters))
      
      #--- Estimate the performance ---#
      object_perf <- perf(object, folds = folds, progressBar = FALSE, 
                          scale_resp = scale_resp, setseed = setseed)
      
      MSEP[ , i ] <- apply(object_perf$MSEP, 1, sum)
      
    }
    
    #-- Find optimal tuning parameters --#
    mMSEP <- min(MSEP)
    
    colnames(MSEP) <- paste0(group_seq)
    rownames(MSEP) <- paste0("ncomp ", 1:max(ncomp))
    
    mMSEPcol <- apply(MSEP, 2, min)
    mMSEProw <- apply(MSEP, 1, min)
    
    ncomp_opt <- min( (1:max(ncomp))[mMSEProw == mMSEP])
    keep_opt <- min( group_seq[mMSEPcol == mMSEP])
    
    #-- Return parameters with optimal fit --#
    parameters <- object$parameters
    parameters$ncomp <- ncomp_opt
    
    if(block == "X" || block == "x")
      parameters$keepX[ncomp] <- keep_opt
    
    if(block == "Y" || block == "y")
      parameters$keepY[ncomp] <- keep_opt
    
    
    res <- list(
      MSEP = MSEP, ncomp_opt = ncomp_opt, keep_opt = keep_opt, parameters = parameters, folds = object_perf$folds
    )
    return(res)
  }


