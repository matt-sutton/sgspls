
perf = function(object, ...) UseMethod("perf")

## A quick function to calulate repeated cross validation for the perf function
rperf <- function(object, folds = 10, progressBar = T, nrepeats = 10, choicesetseed = 1){
 if("sgspls.groups" %in% class(object)){
    nrows = ncol(object$parameters$Y)
    ncols = length(object$group.seq)
    parameters = object$group.seq
  } else {
    nrows = ncol(object$Y)
    ncols = object$ncomp
    parameters = 1:ncols
  }

  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)

  choicesetseed <- choicesetseed

  mseps = matrix(0, nrow = nrows, ncol = ncols)
  n <- nrow(object$parameters$Y)
  press.mat <- array(NA, c(n*nrepeats, nrows, ncols))
  for ( r in 1:nrepeats){

    temp.valid <- perf(object, validation = "Mfold", criterion = "MSEP", folds = folds, setseed = choicesetseed, progressBar = F)
    mseps <- mseps + temp.valid$MSEP
    if(!is.null(temp.valid$press.mat)) press.mat[(1+(r-1)*n):(r*n), ,] <- temp.valid$press.mat
    choicesetseed <- choicesetseed + 1

    if (progressBar == TRUE)
      setTxtProgressBar(pb, r/nrepeats)
  }
  mseps <- scale(mseps,center = F,scale = rep(nrepeats,ncols))
  attr(mseps, "scaled:scale") = NULL
  sMSEP = colSums(mseps)
  opt.paramater = parameters[which.min(sMSEP)]
  return(list(MSEP=mseps, press.mat=press.mat, opt.paramater = opt.paramater, press.mat=press.mat, bic = temp.valid$BIC))
}


## Perf for sgspls goes across the number of components
perf.sgspls <- function (object, criterion = c("all", "MSEP", "R2", "Q2"), validation = c("Mfold",
                                                                           "loo"), folds = 10, progressBar = TRUE, setseed = 1, ...) {
  set.seed(setseed)
  X = object$X#scale(object$X, center = T, scale = object$parameters$scale.x)
  Y = object$Y #scale(object$Y, center = T, scale = object$parameters$scale.y)

  ncomp = object$ncomp

  pls.parameters = object$parameters

  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)

  res = list()
  validation = match.arg(validation)
  featuresX = featuresY = list()
  for (k in 1:ncomp) {
    featuresX[[k]] = featuresY[[k]] = NA
  }
  if(is.null(folds))
    stop("Enter a number of folds")
  if (length(dim(X)) != 2)
    stop("'X' must be a numeric matrix for validation.")
  if (object$call$mode == "canonical")
    stop("mode should be set to regression")
  if (any(criterion == "Q2") & ncomp == 1)
    stop("'ncomp' must be > 1 for Q2 criterion.")
  if (any(is.na(X)) || any(is.na(Y)))
    stop("Missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.")
  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n)
        stop("Invalid number of folds.")
      if (length(unique(unlist(folds))) != n)
        stop("Invalid folds.")
      M = length(folds)
    }
    else {
      if (is.null(folds) || !is.numeric(folds) || folds <
          2 || folds > n)
        stop("Invalid number of folds.")
      else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n))
      }
    }
  }
  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
  }
  RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
  RSS.indiv = array(NA, c(n, q, ncomp + 1))
  PRESS.inside = Q2.inside = matrix(nrow = ncomp, ncol = q)
  if (any(criterion %in% c("all", "MSEP", "R2", "Q2"))) {
    press.mat = Ypred = array(NA, c(n, q, ncomp))
    MSEP = R2 = matrix(NA, nrow = q, ncol = ncomp)
    rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
    dimnames(press.mat)[[2]] = as.list(colnames(Y))
    stop.user = FALSE
    if (progressBar == TRUE)
      pb <- txtProgressBar(style = 3)
    for (i in 1:M) {
      if (progressBar == TRUE)
        setTxtProgressBar(pb, i/M)
      omit = folds[[i]]
      if (length(omit) == 1)
        stop.user = TRUE
      X.train = X[-omit, ]
      Y.train = Y[-omit, ]
      X.test = matrix(X[omit, ], nrow = length(omit))
      Y.test = matrix(Y[omit, ], nrow = length(omit))

      pls.parameters$X = X.train
      pls.parameters$Y = Y.train

      spls.res = do.call(sgspls, args = pls.parameters)

      Y.hat =  predict(spls.res, X.test)$predict
      for (h in 1:ncomp) {
        Ypred[omit, , h] = Y.hat[, , h]
        press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
        RSS.indiv[omit, , h + 1] = (Y.test - Y.hat[,
                                                   , h])^2
      }
    }
    if (stop.user == TRUE & validation == "Mfold")
      stop("The folds value was set too high to perform cross validation. Choose validation = \"loo\" or set folds to a lower value")
    for (h in 1:ncomp) {
      MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2,
                        mean, na.rm = TRUE)
      R2[, h] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
      if (q > 1) {
        RSS[h + 1, ] = t(apply(RSS.indiv[, , h + 1],
                               2, sum))
        PRESS.inside[h, ] = colSums(press.mat[, , h],
                                    na.rm = TRUE)
      }
      else {
        RSS[h + 1, q] = sum(RSS.indiv[, q, h + 1])
        PRESS.inside[h, q] = sum(press.mat[, q, h], na.rm = TRUE)
      }
      Q2.inside[h, ] = 1 - PRESS.inside[h, ]/RSS[h, ]
    }
    colnames(MSEP) = colnames(R2) = rownames(Q2.inside) = paste("ncomp",
                                                                c(1:ncomp), sep = " ")
    rownames(MSEP) = rownames(R2) = colnames(Q2.inside) = colnames(Y)
    if (q == 1)
      rownames(MSEP) = rownames(R2) = ""
    if (ncomp > 1) {
      if (q > 1) {
        Q2.total = 1 - rowSums(PRESS.inside, na.rm = TRUE)/rowSums(RSS[-(ncomp +
                                                                           1), ], na.rm = TRUE)
      }
      else {
        Q2.total = t(1 - PRESS.inside/RSS[-(ncomp + 1),
                                          ])
      }
    }
    else {
      Q2.total = NA
    }
    names(Q2.total) = paste("comp", 1:ncomp, sep = " ")
  }


  ### BIC type statistic

  if (any(criterion %in% c("all", "BIC"))) {
    # Calculations for scaled versions
    #     means.X = attr(X, "scaled:center")
    #     means.Y = attr(Y, "scaled:center")
    #     sigma.X = attr(X, "scaled:scale")
    #     sigma.Y = attr(Y, "scaled:scale")
    #     if(is.null(sigma.X)) sigma.X <- 1
    #     if(is.null(sigma.Y)) sigma.Y <- 1
    #
    #     X.org <- scale(X,center=FALSE,scale=1/sigma.X)
    #     X.org <- scale(X.org,center=-means.X,scale=FALSE)
    #     Y.org <- scale(Y,center=FALSE,scale=1/sigma.Y)
    #     Y.org <- scale(Y.org,center=-means.Y,scale=FALSE)

    Y.hat <- predict(object, X)$predict

    BIC <- rep(0, ncomp)

    # nonzY <-  list(which(abs(object$loadings$Y)>0))
    # if(ncomp > 1){
    #   nonzX <-  apply(object$loadings$X, 2, function(x) which(abs(x)>0))
    #   nonzX <- Reduce(union, nonzX, accumulate = T)
    #   nonzY <- apply(object$loadings$Y, 2, function(x) which(abs(x)>0))
    #   nonzY <- Reduce(union, nonzY, accumulate = T)
    # }
    nonzX <- NULL

    for (h in 1:ncomp) {
      nonzX <- c(nonzX, which(abs(object$loadings$X[,h])>0))
      nonzX <- unique(nonzX)
      kappa <- length(nonzX)
      rss.bic <- sum((Y - Y.hat[,,h])^2)

      BIC[ h ] = (n*q)*log(rss.bic/(n*q)) + (2+kappa)*log(n*q)  + (n*q) + n*q*log(2*pi) ## additional terms to match stats::BIC call
    }
  }

  if (progressBar == TRUE)
    cat("\n")
  list.featuresX = list.featuresY = list()
  for (k in 1:ncomp) {
    remove.naX = which(is.na(featuresX[[k]]))
    remove.naY = which(is.na(featuresY[[k]]))
    list.featuresX[[k]] = sort(summary(as.factor(featuresX[[k]][-remove.naX]))/M,
                               decreasing = TRUE)
    list.featuresY[[k]] = sort(summary(as.factor(featuresY[[k]][-remove.naY]))/M,
                               decreasing = TRUE)
  }
  features.finalX = features.finalY = list()
#   if(ncomp > 1){
#     for (k in 1:ncomp) {
#       if(packageVersion("mixOmics") >= 5.2){
#         features.finalX[[k]] = selectVar(convert.package(object), comp = k)$X$value
#         features.finalY[[k]] = selectVar(convert.package(object), comp = k)$Y$value
#       } else {
#         features.finalX[[k]] = selectVar(convert.package(object), comp = k)$value.X
#         features.finalY[[k]] = selectVar(convert.package(object), comp = k)$value.Y
#       }
#     }
#     names(features.finalX) = names(features.finalY) = names(list.featuresX) = names(list.featuresX) = paste("comp",
#                                                                                                             1:ncomp)
#   }

  if (any(criterion %in% c("all", "MSEP"))) {
    res$MSEP = MSEP
    res$press.mat = press.mat
    res$RSS.indiv = RSS.indiv
    res$PRESS.inside = PRESS.inside
    res$RSS = RSS
    }
  if (any(criterion %in% c("all", "R2"))) {
    res$press.mat = press.mat
    res$RSS.indiv = RSS.indiv
    res$PRESS.inside = PRESS.inside
    res$RSS = RSS
    res$R2 = R2
  }
  if (any(criterion %in% c("all", "Q2"))) {
    res$press.mat = press.mat
    res$RSS.indiv = RSS.indiv
    res$PRESS.inside = PRESS.inside
    res$RSS = RSS
    res$Q2 = t(Q2.inside)
    res$Q2.total = Q2.total
  }
  if (any(criterion %in% c("all", "BIC")))
    res$BIC = BIC
  res$features$stable.X = list.featuresX
  res$features$stable.Y = list.featuresY
  res$features$final.X = features.finalX
  res$features$final.Y = features.finalY
  res$segments = folds
  method = "pls.mthd"
  class(res) = c("perf", method)
  return(invisible(res))
}


# perf for testing across multiple tuning parameters (groups) for a fixed fold.

perf.sgspls.groups <- function (object, validation = c("Mfold", "loo"), folds = 10, progressBar = TRUE, setseed = 1, ...) {
  set.seed(setseed)

  group.seq = object$group.seq
  nonzero = object$nonzero
  ngroups = length(group.seq)
  parameters <- object$parameters
  X = parameters$X#scale(parameters$X, center = T, scale = object$parameters$scale.x)
  Y = parameters$Y #scale(parameters$Y, center = T, scale = object$parameters$scale.y)

  object <- object$tuningparameters
  object$group.seq <- group.seq
  object$pls.obj <- parameters

  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)

  res = list()
  validation = match.arg(validation)

  if(is.null(folds))
    {stop("Enter a number of folds")}
  if (length(dim(X)) != 2)
    {stop("'X' must be a numeric matrix for validation.")}
  if (validation == "Mfold") {
    if (is.list(folds)) {
      if (length(folds) < 2 | length(folds) > n) {
        stop("Invalid number of folds.")
      }
      if (length(unique(unlist(folds))) != n) {
        stop("Invalid folds.")
      }
      M = length(folds)
    }    else {
      if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n) {
        stop("Invalid number of folds.")
      }  else {
        M = round(folds)
        folds = split(sample(1:n), rep(1:M, length = n))
      }
    }
  }  else {
    folds = split(1:n, rep(1:n, length = n))
    M = n
  }
  PRESS.inside = matrix(nrow = ngroups, ncol = q)
  press.mat = Ypred = array(NA, c(n, q, ngroups))
  MSEP = R2 = matrix(NA, nrow = q, ncol = ngroups)
  rownames(MSEP) = rownames(R2) = colnames(Y)
  dimnames(press.mat)[[2]] = as.list(colnames(Y))
  stop.user = FALSE
  if (progressBar == TRUE)
    pb <- txtProgressBar(style = 3)
  #cat(M)
  for (i in 1:M) {
    if (progressBar == TRUE)
      setTxtProgressBar(pb, i/M)
    omit = folds[[i]]
    if (length(omit) == 1)
      stop.user = TRUE
    X.train = X[-omit, ]
    Y.train = Y[-omit, ]
    X.test = matrix(X[omit, ], nrow = length(omit))
    Y.test = matrix(Y[omit, ], nrow = length(omit))

    object$pls.obj$X = X.train
    object$pls.obj$Y = Y.train
    object$newdata = X.test

    spls.res = do.call(sgspls.tune.groups, args = object)
    Y.hat = array(spls.res$predict, dim = c(nrow(X.test), q, ngroups))

    for (h in 1:ngroups) {
        Ypred[omit, , h] = Y.hat[, , h]
        press.mat[omit, , h] = (Y.test - Y.hat[, , h])^2
        }
  }

  if (stop.user == TRUE & validation == "Mfold") {
    stop("The folds value was set too high to perform cross validation. Choose validation = \"loo\" or set folds to a lower value")
  }
  for (h in 1:ngroups) {
      MSEP[, h] = apply(as.matrix(press.mat[, , h]), 2, mean, na.rm = TRUE)
      R2[, h] = (diag(cor(Y, Ypred[, , h], use = "pairwise")))^2
      if (q > 1) {
        PRESS.inside[h, ] = colSums(press.mat[, , h],
                                    na.rm = TRUE)
      }
      else {
        PRESS.inside[h, q] = sum(press.mat[, q, h], na.rm = TRUE)
      }
    }
    colnames(MSEP) = colnames(R2) = paste(group.seq)
    rownames(MSEP) = rownames(R2) = colnames(Y)
    if (q == 1){
      rownames(MSEP) = rownames(R2) = ""
    }
    if (progressBar == TRUE)
      cat("\n")


    res$group.seq = group.seq
    res$nonzero = nonzero
    res$keepX = group.seq[which.min(colSums(MSEP))]
    res$MSEP = MSEP
    res$R2 = R2
    res$press.mat = press.mat
    res$PRESS.inside = PRESS.inside
    method = "pls.mthd"
    class(res) = c("perf", method)
    return(invisible(res))
    }

