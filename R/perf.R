
perf = function(object, ...) UseMethod("perf")


#' Performance evaluation of sgsPLS objects
#'
#' Function to evaluate the performance of the fitted the PLS models using various criteria. 
#' Evaluation is made for each component in the PLS object and can only be evaluated for regression PLS.
#'
#' @param object Object of class inheriting from \code{"sgspls"}. 
#' The function will retrieve some key parameters stored in that object.
#' @param validation What kind of cross validation to use, matching one of \code{"Mfold"} or 
#' \code{"loo"} (leave one out). Default is \code{"Mfold"}.
#' @param folds Number of folds to use in the cross validation.
#' @param BIC Return the BIC type criterion for large sample sizes (see paper for details). 
#' Note that this is not an actual BIC criterion.
#' @param progressBar Logical to show a progress bar while computing the performances
#' @param scale_resp Logical to scale the responses. This is useful if comparing the fit on 
#' multiple responses.
#' 
#' @export
#' @return \code{perf} returns a list that contains the following performance measures: 
#' 
#' \item{MSEP}{A matrix of Mean Square Error of Prediction (MSEP) estimates by cross validation. The penalty is defined
#' as: \deqn{MSEP = 1/n \sum \sum (f_k(x_i) - y_i)^2} see the references for details. }
#' \item{PRESS0}{A vector of cross validated Predictive Residual Sum of Squares (PRESS) values. 
#' Each column corresponds to a response. Matches the PLS package.}
#' \item{PRESS}{A matrix of cross validated Predictive Residual Sum of Squares (PRESS) values. 
#' Each row contains the values for a different component and each column corresponds to a response.} 
#' \item{R2}{a matrix of \eqn{R^2} values of the \eqn{Y}-variables with \code{ncomp} components} 
#' \item{BIC}{A BIC type criterion for large samples (see references for details). 
#' Note that this is not an actual BIC criterion.} 
#' \item{cvPred}{an array with the cross-validated predictions.}
#' \item{folds}{A list of the folds used in the cross validation.}
#' 
#' @references Mevik, Bjørn-Helge, and Henrik René Cederkvist. 2004. 
#'   Mean Squared Error of Prediction (MSEP) Estimates for Principal Component Regression (PCR) 
#'   and Partial Least Squares Regression (PLSR). 
#'   \emph{Journal of Chemometrics} \bold{18} (9). John Wiley & Sons, Ltd.:422–29.
#'   
#' @seealso Tuning functions \code{\link[sgspls]{calc_pve}},
#' \code{\link[sgspls]{tune_sgspls}}, \code{\link[sgspls]{tune_groups}}. 
#' Model performance and estimation  \code{\link[sgspls]{predict.sgspls}}, \code{\link[sgspls]{perf.sgspls}}, \code{\link[sgspls]{coef.sgspls}}
#' 
#' @examples
#'
#' set.seed(1)
#' n = 50; p = 500; 
#' 
#' size.groups = 30; size.subgroups = 5
#' groupX <- ceiling(1:p / size.groups)
#' subgroupX <- ceiling(1:p / size.subgroups)
#' 
#' X = matrix(rnorm(n * p), ncol = p, nrow = n)
#' 
#' beta <- rep(0,p)
#' bSG <- -2:2; b0 <- rep(0,length(bSG))
#' betaG <- c(bSG, b0, bSG, b0, bSG, b0)
#' beta[1:size.groups] <- betaG
#' 
#' y = X %*% beta + 0.1*rnorm(n)
#' 
#' model <- sgspls(X, y, ncomp = 3, mode = "regression", keepX = 1,
#'                 groupX = groupX, subgroupX = subgroupX,
#'                 indiv_sparsity_x = 0.8, subgroup_sparsity_x = 0.15)
#'
#' model_perf <- perf(model, folds = 5)
#'
#' # Check model performance
#' model_perf$MSEP
#' 

perf.sgspls <- 
  function (object, validation = c("Mfold","loo"), folds = 10, BIC = FALSE,
            progressBar = TRUE, setseed = FALSE, scale_resp = FALSE) {
    
  if(setseed) set.seed(setseed)
  pls_parameters = object$parameters

  X = pls_parameters$X
  Y = pls_parameters$Y
  
  ncomp = pls_parameters$ncomp
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  
  validation = match.arg(validation)
  
  #---------------------#
  #-- check perf args --#
  
  if(is.null(folds))
    stop("Enter a number of folds")
  if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix for validation.")
  if (object$parameters$mode == "canonical")
    stop("mode should be set to regression")
  if (any(is.na(X)) || any(is.na(Y)))
    stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
  
  
  #---------------------#
  #-- define folds -----#
  
  if (validation == "Mfold") 
    {
      if (is.list(folds)) 
        {
        
          if (length(folds) < 2 | length(folds) > n)
          stop("Invalid number of folds.")
        
          if (length(unique(unlist(folds))) != n)
          stop("Invalid folds. The total number of samples in folds must be equal to ",n,".")
        
          if (length(unique(unlist(folds))) != n)
          stop("Invalid folds. Repeated samples in folds.")
        
          M = length(folds)
          
        } else {
          
          if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
            {
              stop("Invalid number of folds.")
            } else {
              M = round(folds)
              folds = split(sample(1:n), rep(1:M, length = n))
            }
        }
    } else {
      folds = split(1:n, rep(1:n, length = n))
      M = n
      }
  
  #-- set-up progress bar --#
  
  if ( progressBar == TRUE )
  {
    pb <- txtProgressBar(style = 3)
  }
  
  #-- initialise -----#

  RSS <- rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
  PRESS <- Q2 <- MSEP <- R2 <- matrix(nrow = ncomp, ncol = q)
  BIC_mat <- matrix(nrow = ncomp, ncol = 1)
  
  Ypred <- RSS_indiv <- array(NA, c(n, q, ncomp))
  
  #--------------------------------#
  #-- loop for cross validation ---#
  for (i in 1:M) {
    
    if (progressBar == TRUE)
    {
      setTxtProgressBar(pb, i/M)
    }
    
    omit = folds[[i]]
    X_train = X[-omit, ]
    Y_train = Y[-omit, ]
    X_test = matrix(X[omit, ], nrow = length(omit))
    Y_test = matrix(Y[omit, ], nrow = length(omit))
    
    pls_parameters$X = X_train
    pls_parameters$Y = Y_train
    
    #-- fit the PLS method using the cv dataset ---#
    spls_res = do.call(sgspls, args = pls_parameters)
    Y_hat =  predict(spls_res, X_test)
    
    #-- loop through h components ---#
    for (h in 1:ncomp) 
      {
        Ypred[omit, , h] = Y_hat[, , h]
        RSS_indiv[omit, , h] <- (Y_test - Y_hat[, , h])^2
    }
  }
  
  #-------------------------------------#
  #-- calculate performance measures ---#
  
  PRESS0 <- apply(Y, 2, var) * n^2/(n-1)
  nonzX <- NULL
  Y_hat <- predict(object, X)
  
  for (h in 1:ncomp) 
    {
    
    #-- calculate MSEP ---#
    MSEP[ h, ] = apply(as.matrix(RSS_indiv[, , h]), 2, mean, na.rm = TRUE)
    
    #-- calculate R2 ---#
    R2[ h, ] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
    
    #-- calculate Q2 ---#
    PRESS[ h, ] <- colSums(as.matrix(RSS_indiv[, , h]), na.rm = TRUE)
    
    #-- BIC type statistic --#
    if ( BIC ){
      
      nonzX <- c(nonzX, which(abs(object$weights$X[,h])>0))
      kappa <- length(unique(nonzX))
      rss_bic <- sum((Y - Y_hat[,,h])^2)
      
      BIC_mat[ h ] = (n*q)*log(rss_bic/(n*q)) + (2+kappa)*log(n*q)  + (n*q) + n*q*log(2*pi) 
      ## additional terms to match stats::BIC call
      }
  }
  
  #-- Get multiple responses on the same scale (MSEP and PRESS) --#
  if(scale_resp)
    MSEP <- scale(MSEP, center = F, scale = diag(var(Y)))
  
  #-- Add dimnames --#
  rownames(MSEP) <- rownames(R2) <- paste("ncomp", c(1:ncomp), sep = " ")
  colnames(MSEP) <- colnames(R2) <- colnames(Y)

  #-- Progress bar update --#  
  if (progressBar == TRUE)
    cat("\n")
  
  res <- list(
    MSEP = MSEP, PRESS0 = PRESS0, PRESS = PRESS, R2 = R2,
    BIC = BIC_mat, cvPred = Ypred, folds = folds
  )
  
  return(res)
  }

