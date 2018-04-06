#' Predict Method for sgspls
#'
#' Predicted values based on sparse group subgroup PLS. New responses are predicted using a fitted model and a new matrix of observations.
#'
#' @param object Object of class inheriting from \code{"sgspls"}. 
#' @param newdata Data matrix in which to look for for explanatory variables to be used for prediction.
#' @param ... Not currently used.
#' 
#' @export
#' @return \code{perf} returns a list that contains the following performance measures: 
#' \code{predict} function produces predicted values, obtained by evaluating the sparse group subgroup PLS. 
#' The prediction values are calculated based on the regression coefficients of \code{object$Y} onto \code{object$variates$X}.
#' 

predict.sgspls <-
  function(object, newdata,  ...)  {
    
    # predict function for sgspls

    newdata <-  as.matrix(newdata)
    nobs <- nrow(newdata)
    p <- ncol(newdata)
    nresp <- ncol(object$parameters$Y)
    npred <- ncol(object$parameters$X)
    ncomp <- object$parameters$ncomp
    
    #-- validation des arguments --#
    if (missing(newdata)){
      stop("No new data available.")
    }
    
    if (length(dim(newdata)) == 0) {
      if (length(newdata) != p)
        stop("'newdata' must be a numeric matrix with ncol = ", p,
             " or a vector of length = ", p, ".")
      dim(newdata) = c(1, p)
    }
    

    B <- array(0, dim = c(npred, nresp, ncomp))
    B_coef <- coef(object, type = "coefficients")
    B <- B_coef$B
    B0 <- B_coef$B0
    
    pred <- array(dim = c(nobs, nresp, ncomp))
    
    for( i in 1:ncomp ){
      pred[,,i] <- newdata%*%B[ , , i] + matrix(rep(B0[i,], each = nobs), nrow = nobs)
    }

    return(pred)
  }

