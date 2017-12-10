#' Sparse Group Subgroup PLS
#'
#' Fit a PLS model to two blocks of data via the sparse group subgroup Partial
#' Least Squares (sgspls) algorithm. The sgspls algorithm enables selection of
#' variables at the group, subgroup and single feature levels.
#'
#' @param X A matrix of regressors (n x p). By default the matrix will be
#'   centered to have mean zero.
#' @param Y A matrix of continuous responses (n x q). By default the matrix will
#'   be centered to have mean zero.
#' @param ncomp The number of components to include in the model.
#' @param mode A character string. What type of PLS algorithm to use, matching
#'   one of "regression", "canonical". See Details.
#' @param tol A positive real tolerance for the PLS algorithm.
#' @param groupX A vector describing the group details of the X variable. (see
#'   example in Details).
#' @param groupY A vector describing the group details of the Y variable. (see
#'   example in Details).
#' @param subgroupX A vector describing the subgroup details of the X block (see
#'   example in Details).
#' @param subgroupY A vector describing the subgroup details of the Y block (see
#'   example in Details).
#' @param keepX Numeric vector of length \code{ncomp}, the number of groups to
#'   select in \eqn{X}-loadings. Default selects all groups.
#' @param keepY Numeric vector of length \code{ncomp}, the number of groups to
#'   select in \eqn{Y}-loadings. Default selects all groups.
#' @param indiv.sparsity.x Individual sparisty parameter (value between 0 and 1)
#'   related to the sparisty within subgroups for the \eqn{X} block.
#' @param indiv.sparsity.y Individual sparisty parameter (value between 0 and 1)
#'   related to the sparisty within subgroups for the \eqn{Y} block.
#' @param subgroup.sparsity.x Sub-group sparisty parameter (value between 0 and
#'   1) related to the number of subgroups selected for the PLS \eqn{X} weights.
#' @param subgroup.sparsity.y Sub-group sparisty parameter (value between 0 and
#'   1) related to the number of subgroups selected for the PLS \eqn{Y} weights.
#' @param scale.x Scale predictors by their standard deviation.
#' @param scale.y Scale responses by their standard deviation.
#'
#' @export
#' @return \code{sgspls} returns an object of class \code{"sgspls"}, a list that
#' contains the following components: \item{X}{The centered (and standardized)
#' original predictor matrix.} \item{Y}{The centered (and standardized) original
#' predictor matrix.} \item{ncomp}{The number of components included in the
#' model.} \item{mode}{The type of PLS algorithm.} \item{keepX}{Number of
#' \eqn{X} groups kept in the model on each component.} \item{keepY}{Number of
#' \eqn{Y} groups kept in the model on each component.} \item{mat.c}{Matrix of
#' coefficients to be used internally by \code{predict}.} \item{variates}{List
#' containing the variates.} \item{loadings}{List containing the estimated
#' loadings for the \eqn{X} and \eqn{Y} variates.} \item{names}{List containing
#' the names to be used for individuals and variables.} \item{tol}{The tolerance
#' used in the iterative algorithm, used for subsequent S3 methods.}
#' \item{max.iter}{The maximum number of iterations, used for subsequent S3
#' methods.} \item{iter}{Vector containing the number of iterations for
#' convergence in each component.} \item{parameters}{List containing the
#' parameters of the model that was fitted.}
#'
#' @author Matthew Sutton \email{m5.sutton@hdr.qut.edu.au}
#' @references Liquet Benoit, Lafaye de Micheaux, Boris Hejblum, Rodolphe
#'   Thiebaut. A group and Sparse Group Partial Least Square approach applied in
#'   Genomics context. \emph{Submitted}.
#'
#'   L\^e Cao, K.-A., Martin, P.G.P., Robert-Grani\'e, C. and Besse, P. (2009).
#'   Sparse canonical methods for biological data integration: application to a
#'   cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.
#'
#'   Shen, H. and Huang, J. Z. (2008). Sparse principal component analysis via
#'   regularized  low rank matrix approximation. \emph{Journal of Multivariate
#'   Analysis} \bold{99}, 1015-1034.
#'
#'   Tenenhaus, M. (1998). \emph{La r\'egression PLS: th\'eorie et pratique}.
#'   Paris: Editions Technic.
#'
#' @seealso Tuning functions \code{\link{sgspls.tune.groups}},
#' \code{\link{sgspls.tune}}, \code{\link{predict}}, \code{\link{perf}} and
#' functions from \code{mixOmics} package: \code{summary}, \code{plotIndiv},
#' \code{plotVar}, \code{plot3dIndiv}, \code{plot3dVar}.
#'
#' @examples
#'
#' n = 100; p = 200; size.groups = 25; size.subgroups = 5
#' groupX <- ceiling(1:p / size.groups)
#' subgroupX <- ceiling(1:p / size.subgroups)
#' X = matrix(rnorm(n * p), ncol = p, nrow = n)
#'
#' beta <- rep(0,p)
#' bSG <- -2:2; b0 <- rep(0,length(bSG))
#' betaG <- c(bSG, b0, bSG, b0, b0)
#' beta[1:size.groups] <- betaG
#'
#' y = X %*% beta + 0.1*rnorm(n)
#'
#' model <- sgspls(X, y, ncomp = 3, mode = "regression", keepX = 1,
#'                        groupX = groupX, subgroupX = subgroupX,indiv.sparsity.x = 0.8,
#'                        subgroup.sparsity.x = 0.15, penalty = "sgsgPLS")
#' predict(model, X)
#' coef(model)
#'
#' \dontrun{
#' cbind(model.sgsplsR$B.hat[,,3], beta)[1:30,]
#' }


sgspls <-
  function(X, Y, ncomp = 2, mode = c("regression", "canonical"),max.iter = 500,  keepX = NA,
           keepY = NA, groupX = rep(1,ncol(X)), groupY = rep(1,ncol(Y)),
           subgroupX = rep(1,ncol(X)), subgroupY = rep(1,ncol(Y)),
           indiv.sparsity.x = rep(0,ncomp), subgroup.sparsity.x = rep(0,ncomp),
           indiv.sparsity.y = rep(0,ncomp), subgroup.sparsity.y = rep(0,ncomp),
           tol = 1e-12, scale.x = T, scale.y = T, newtry = F,
           ...)
  {

    # Repeat parameters if left missing
    ngY <- length(unique(groupY))
    ngX <- length(unique(groupX))
    if(is.na(keepY[1])) keepY <- ngY
    if(is.na(keepX[1])) keepX <- ngX
    keepX = add.penalty.param(keepX,ncomp)
    keepY = add.penalty.param(keepY,ncomp)
    indiv.sparsity.x = add.penalty.param(indiv.sparsity.x,ncomp)
    subgroup.sparsity.x = add.penalty.param(subgroup.sparsity.x,ncomp)
    indiv.sparsity.y = add.penalty.param(indiv.sparsity.y,ncomp)
    subgroup.sparsity.y = add.penalty.param(subgroup.sparsity.y,ncomp)

    # Set the alphas for the sparsity values
    # There must be positive penalty on the group so that the keepX and keepY parameters work
    alpha1.x <- 1 - rowSums(cbind(matrix(indiv.sparsity.x, ncol = 1), matrix(subgroup.sparsity.x, ncol = 1)))

    if(any(alpha1.x == 0)){
      ind <- which(alpha1.x == 0)
      alpha1.x[ind] <- 1e-06
      for (i in ind){
        if(indiv.sparsity.x[ind] < subgroup.sparsity.x[ind]) {
          subgroup.sparsity.x[ind] = subgroup.sparsity.x[ind] - 1e-06
        } else{
          indiv.sparsity.x[ind] = indiv.sparsity.x[ind] - 1e-06
        }
      }
    }
    alpha2.x <- subgroup.sparsity.x
    alpha3.x <- indiv.sparsity.x

    alpha1.y <- 1 - rowSums(cbind(indiv.sparsity.y, subgroup.sparsity.y))
    alpha1.y <- pmax(alpha1.y, 1e-6)

    if(any(alpha1.y == 0)){
      ind <- which(alpha1.x == 0)
      alpha1.x[ind] <- 1e-06
      for (i in ind){
        if(indiv.sparsity.x[ind] < subgroup.sparsity.x[ind]) {
          subgroup.sparsity.x[ind] = subgroup.sparsity.x[ind] - 1e-06
        } else{
          indiv.sparsity.x[ind] = indiv.sparsity.x[ind] - 1e-06
        }
      }
    }

    alpha2.y <- subgroup.sparsity.y
    alpha3.y <- indiv.sparsity.y

    if(any(round(rowSums(cbind(alpha1.x, alpha2.x, alpha3.x)),6) != 1)) stop("Sparsities don't sum to 1 on X")
    if(any(round(rowSums(cbind(alpha1.y, alpha2.y, alpha3.y)),6) != 1)) stop("Sparsities don't sum to 1 on Y")

    ## Specialist method to preform sgspls with lambda and alpha's specified on X
    # additional Args #
    lambda.x = rep(0,ncomp)
    lambda.y = rep(0,ncomp)

    add.args = list(...)
    if(!is.null(add.args$lambda.x)){
      lambda.x <- add.args$lambda.x
      alpha1.x <- add.args$alpha1.x
      alpha2.x <- add.args$alpha2.x
      alpha3.x <- add.args$alpha3.x
    }

    # Check the data #
    X = as.matrix(X)
    Y = as.matrix(Y)
    p = ncol(X); q = ncol(Y); n = nrow(X);

#    if (p == 1 || q < 1) stop("\nBlock X must contain more than one column")
#    if (nrow(Y) != n)  stop("\nX and Y have different number of rows")
    if (any(!is.finite(Y))) stop("\ninfinite, NA or NaN values in X or Y")

    if (is.null(colnames(X))) colnames(X) = paste(rep("X",p), 1:p, sep="")
    if (is.null(rownames(X))) rownames(X) = rep(1:n)

    if (is.null(colnames(Y))) colnames(Y) = paste(rep("Y",q), 1:q, sep="")
    if (is.null(rownames(Y))) rownames(Y) = rep(1:n)

    # Scaleing similar to the pls package #
    Y_h <- scale(Y, center = T, scale = scale.y)
    X_h <- scale(X, center = T, scale = scale.x)

    mat.t <- matrix(nrow=n, ncol = ncomp) # X variate
    mat.u <- matrix(nrow=n, ncol = ncomp) # Y variate
    load.u <- matrix(nrow = p, ncol = ncomp)
    load.v <- matrix(nrow = q, ncol = ncomp)
    mat.c <-matrix(nrow = p, ncol = ncomp)
    mat.d <- matrix(nrow = q, ncol = ncomp)
    mat.e <- matrix(nrow = q, ncol = ncomp)

    # Algorithm #
    for(h in 1:ncomp){
      M = t(X_h)%*%Y_h
      svd.M <- try(svd(M,nu = 1,nv = 1))

      if(class(svd.M) == "try-error"){
        require(irlba,quietly = T)
        irlba(M,nu = 1,nv = 1)
      }

      if (svd.M$d[1] < .Machine$double.eps) {
        cat("Too many PLS components, stoped at ",h-1," components",svd.M$d[1])
        ncomp = h-1;
        scores = scores[,1:(h-1)]; weights = weights[,1:(h-1)]; loadings = loadings[,1:(h-1)]
        path.coeff = path.coeff[1:(h-1)]; comps = comps[1:(h-1)]
        break
      }

      # Calaulate the loadings
      res.load <- group.sparse.subgroup.penalty(M = M,svd.M = svd.M,keepX = keepX[h],keepY = keepY[h],groupX = groupX,
                                                groupY = groupY,subgroupX = subgroupX,subgroupY = subgroupY,alpha1.x = alpha1.x[h],
                                                alpha2.x = alpha2.x[h],alpha1.y = alpha1.y[h],alpha2.y = alpha2.y[h],tol, max.iter,
                                                lambda.x = lambda.x[h], lambda.y = lambda.y[h],newtry=newtry)
      load.u[,h] <- res.load$u
      load.v[,h] <- res.load$v

      # calculate loadings, scores and deflate X & Y
      res.deflat <- deflate.pls(X=X_h,Y=Y_h,res.load$u,res.load$v,mode=mode)
      mat.t[,h] <- X_h%*%matrix(res.load$u,ncol = 1)
      mat.u[,h] <- Y_h%*%matrix(res.load$v,ncol = 1)
      X_h <- res.deflat$X_h
      Y_h <- res.deflat$Y_h
      mat.c[,h] <- res.deflat$c
      if (mode=="regression") mat.d[,h] <- res.deflat$d else mat.e[,h] <- res.deflat$e
    }

    cl = match.call()
    if (is.null(keepY)){
      keepY <- rep(length(unique(groupY)),ncomp)
    }

    # Save parameters for quick code in perf and tunning functions
    parameters <- list(X=X,
                       Y=Y,
                       ncomp=ncomp,
                       max.iter = max.iter,
                       mode=mode,
                       keepX=keepX,
                       keepY=keepY,
                       groupX = groupX, groupY = groupY,
                       subgroupX = subgroupX, subgroupY = subgroupY,
                       indiv.sparsity.x = indiv.sparsity.x, subgroup.sparsity.x = subgroup.sparsity.x,
                       indiv.sparsity.y = indiv.sparsity.y, subgroup.sparsity.y = subgroup.sparsity.y,
                       tol = tol, scale.x = scale.x, scale.y = scale.y,
                       X.residuals = X_h, Y.residuals = Y_h)

    result <- list(call = cl,X=X,Y=Y, ncomp=ncomp,mode=mode,keepX=keepX,keepY=keepY,mat.c=mat.c,mat.d=mat.d,mat.e=mat.e,loadings = list(X = load.u, Y = load.v),variates = list(X = mat.t, Y = mat.u),
                   names = list(X = colnames(X),Y = colnames(Y), indiv = rownames(X)),tol=tol,max.iter=max.iter,parameters=parameters)
    class(result) = c("sgspls")
    return(invisible(result))
  }
