## Calculate the PVE


## Calculate the number of components to look at:
calc.num.components <- function(X,Y){
  plsmodel <- sgspls(X,Y, ncomp = 15, mode = "regression")
  predmodel <- predict(plsmodel, X)
  PVE = 1 - apply(predmodel$predict, 3, function(x) sum(apply(Y - x,2,var))/sum(apply(Y,2,var)))
  if (max(PVE) > 0.98){
    ncomp <- which(PVE >0.95)[1]
  } else{
    ncomp <- which.min(PVE)[1]
  }
  return(list(ncomp, PVE))
}

per.var.sgspls <- function(model){
  predmodel <- predict(model, model$X)
  PVE = 1 - apply(predmodel$predict, 3, function(x) sum(apply(model$Y - x,2,var))/sum(apply(model$Y,2,var)))
}


sgspls.tune.groups <-
  function(pls.obj,  indiv.sparsity.x=0, subgroup.sparsity.x=0, indiv.sparsity.y=0, subgroup.sparsity.y=0, keepY = NULL,
           group.seq = NULL, newdata = NULL) {

    # Calculate the alphas corresponding to the sparsities
    alpha1.x <- 1 - rowSums(cbind(indiv.sparsity.x, subgroup.sparsity.x))
    alpha1.x <- pmax(alpha1.x, 1e-6)

    alpha2.x <- checkalphas(subgroup.sparsity.x)
    alpha3.x <- checkalphas(indiv.sparsity.x)
    alpha1.y <- 1 - rowSums(cbind(indiv.sparsity.y, subgroup.sparsity.y))
    alpha1.y <- pmax(alpha1.y, 1e-6)

    alpha2.y <- checkalphas(subgroup.sparsity.y)
    alpha3.y <- checkalphas(indiv.sparsity.y)

    # Get the parameters for the current pls object
    if (is.null(pls.obj$X)) stop("Enter a X matrix")
    if (is.null(pls.obj$Y)) stop("Enter a Y matrix")
    X <- pls.obj$X
    Y <- pls.obj$Y
    Y <- matrix(Y, ncol = max(ncol(Y), 1))

    parameters = list()
    parameters$X <- X
    parameters$Y <- Y
    if(is.null(pls.obj$ncomp)) pls.obj$ncomp = 0
    if(is.null(pls.obj$mode)) pls.obj$mode = "regression"
    if(is.null(pls.obj$groupY)) pls.obj$groupY = rep(1,ncol(Y))
    if(is.null(pls.obj$groupX)) pls.obj$groupX = rep(1,ncol(X))
    if(is.null(pls.obj$subgroupY)) pls.obj$subgroupY = rep(1,ncol(Y))
    if(is.null(pls.obj$subgroupX)) pls.obj$subgroupX = rep(1,ncol(X))

    # Fill any missing pls parameters
    if(pls.obj$ncomp == 0) {
      parameters <- pls.obj
      parameters$scale.x <- if(is.null(pls.obj$scale.x))  T  else  pls.obj$scale.x
      parameters$scale.y <- if(is.null(pls.obj$scale.y))  T  else  pls.obj$scale.y
      Xh <- scale(X, center = T, scale = parameters$scale.x);
      Yh <- scale(Y, center = T, scale = parameters$scale.y);
    }  else {
      pls.obj <- do.call(sgspls, args = pls.obj)
      parameters <- pls.obj$parameters
      Xh <- as.matrix(parameters$X.residuals)
      Yh <- as.matrix(parameters$Y.residuals)
    }

    # Set the arguments for finding the pls weights
    p = ncol(X); q = ncol(Y); n = nrow(X)
    M = t(Xh) %*% Yh
    groupX = parameters$groupX; subgroupX = parameters$subgroupX;
    groupY = parameters$groupY; subgroupY = parameters$subgroupY;
    svd.M = svd(M, nu = 1, nv = 1)

    if (is.null(keepY)) keepY <- length(unique(groupY))

    load.args = c(list(M = M, svd.M = svd.M, groupX = groupX, subgroupX = subgroupX, keepY = keepY,
                       alpha1.y = alpha1.y, alpha2.y = alpha2.y, groupY=groupY, subgroupY=subgroupY,
                       alpha1.x = alpha1.x, alpha2.x = alpha2.x))

    ngroups = length(group.seq); nonzero = numeric(ngroups);
    if(is.null(newdata)) newdata <- X
    predict.Y = array( 0, c(nrow(newdata), q, ngroups));

    for(i in 1:ngroups){
      load.args$keepX = group.seq[i]
      loadings <- invisible(do.call(group.sparse.subgroup.penalty, args = load.args))
      # speedup since SVD does not need to be recalculated
      u <- loadings$u
      v <- loadings$v
      nonzero[i] = sum(abs(u)>0)

      # calculate the Beta estimate
      variates <- loadings <- list()
      variates$X <- cbind( pls.obj$variates$X, Xh%*%u )
      mat.c <- cbind( pls.obj$mat.c, t(Xh)%*%matrix(Xh%*%u,ncol=1)/e.norm(Xh%*%u)^2 )
      loadings$X <- cbind( pls.obj$loadings$X, u )
      loadings$Y <- cbind( pls.obj$loadings$Y, v )

      sgsplspredict <- list(X = X, Y = Y, variates = variates, loadings = loadings, mat.c = mat.c, ncomp = parameters$ncomp + 1 )
      sgsplspredict$parameters$scale.x = parameters$scale.x; sgsplspredict$parameters$scale.y = parameters$scale.y
      class(sgsplspredict)<- c("sgspls")

      sgsplspredict <- predict(sgsplspredict, newdata = newdata)

      predict.Y[ , , i] <- sgsplspredict$predict[ , , parameters$ncomp + 1 ]
    }

    tuningparameters <- list(indiv.sparsity.x=indiv.sparsity.x, subgroup.sparsity.x=subgroup.sparsity.x,
                             indiv.sparsity.y=indiv.sparsity.y, subgroup.sparsity.y=subgroup.sparsity.y, keepY = keepY)

    result = list(group.seq = group.seq, nonzero=nonzero, predict = predict.Y, parameters = parameters, tuningparameters = tuningparameters)
    class(result) <- c("sgspls.groups")
    return(result)
  }

sgspls.tune <-
  function(pls.obj, sparsities=NULL, group.seq=NULL, indiv.sparsity.y=0, subgroup.sparsity.y=0, keepY = NULL,
           folds= 10, nrepeats = 1, progressBar=T, setseed = 1) {

    if(is.null(sparsities)){
      sparsities = seq(0.01,0.99,0.01)
      grp = c(seq(0.05,0.95,0.1))

      sparsities.s <- do.call(expand.grid,list(grp,sparsities,sparsities)) # third column to ensure individual sparsity
      sparsities <- NULL
      ind <- which(sparsities.s[,1] %in% grp )
      sparsities.s <- sparsities.s[ind,]
      for(sp in grp){
        ind <- which(apply(sparsities.s,1,sum) == 1 & sparsities.s[,1] == sp )
        ind <- c(tail(ind,1),ind[floor(length(ind)/2)],head(ind,1))
        sparsities <- rbind(sparsities,sparsities.s[ind,])
      }
      sparsities <- unique(sparsities)
      rownames(sparsities) <- NULL
      colnames(sparsities) <- c("Group", "Subgroup","Individual")
    }

    if (progressBar == TRUE)
      pb <- txtProgressBar(style = 3)

    ngroups <- length(group.seq)
    nsparsities <- dim(sparsities)[1]
    results.tuning <- NULL

    for( j in 1:nsparsities){
      sgspls <- sgspls.tune.groups(pls.obj = pls.obj, indiv.sparsity.x = sparsities[j,3], subgroup.sparsity.x = sparsities[j,2], group.seq=group.seq,
                                   indiv.sparsity.y=indiv.sparsity.y, subgroup.sparsity.y=subgroup.sparsity.y, keepY = keepY)

      sgspls.cv <- rperf(sgspls, folds = folds, nrepeats = nrepeats, progressBar = F,choicesetseed = setseed)

      results <- cbind(group.seq, sparsities[j,3], sparsities[j,2], colSums(sgspls.cv$MSEP))

      results.tuning <- rbind(results.tuning, results)
      if (progressBar == TRUE) {
        setTxtProgressBar(pb, j/(nsparsities))
      }
    }
    colnames(results.tuning) <- c("Groups", "Individual", "Subgroup", "MSEP")
    best <- results.tuning[which.min(results.tuning[,4]),]
    keepX <- best[1];
    indiv.sparsity.x <- best[2];
    subgroup.sparsity.x <- best[3]

    # update parameters
    parameters <- sgspls$parameters
    parameters$ncomp <- parameters$ncomp + 1
    parameters$keepX <- c(parameters$keepX, keepX)
    parameters$indiv.sparsity.x <- c(parameters$indiv.sparsity.x, indiv.sparsity.x)
    parameters$subgroup.sparsity.x <- c(parameters$subgroup.sparsity.x, subgroup.sparsity.x)

    parameters$keepY <- c(parameters$keepY, keepY)
    parameters$indiv.sparsity.y <- c(parameters$indiv.sparsity.y, indiv.sparsity.y)
    parameters$subgroup.sparsity.y <- c(parameters$subgroup.sparsity.y, subgroup.sparsity.y)

    rownames(results.tuning) <- paste("alpha",rep(1:length(sparsities[,1]), each = length(group.seq)),sep = "_")

    result <- list(results.tuning = results.tuning,
                   keepX = keepX,
                   indiv.sparsity = indiv.sparsity.x,
                   subgroup.sparsity=subgroup.sparsity.x,
                   parameters = parameters,
                   tuning.parameters = sparsities,
                   folds = folds,
                   min.cv = min(results.tuning[,4]),
                   group.seq = group.seq)

    class(result) = c("cv.sgspls")
    return(result)
    }

