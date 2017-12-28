
# predict function for sgspls

predict.sgspls <-
  function(object, newdata,  ...)  {

    newdata <-  as.matrix(newdata)
    nobs <- nrow(newdata)
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

