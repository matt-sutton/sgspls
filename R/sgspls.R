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
#' @param indiv_sparsity_x Individual sparisty parameter (value between 0 and 1)
#'   related to the sparisty within subgroups for the \eqn{X} block.
#' @param indiv_sparsity_y Individual sparisty parameter (value between 0 and 1)
#'   related to the sparisty within subgroups for the \eqn{Y} block.
#' @param subgroup_sparsity_x Sub-group sparisty parameter (value between 0 and
#'   1) related to the number of subgroups selected for the PLS \eqn{X} weights.
#' @param subgroup_sparsity_y Sub-group sparisty parameter (value between 0 and
#'   1) related to the number of subgroups selected for the PLS \eqn{Y} weights.
#' @param scale.x Scale predictors by their standard deviation.
#' @param scale.y Scale responses by their standard deviation.
#' @param ... additional arguments for low level functionality.
#' @param max.iter How many iterations should be performed? Default is 500.
#'
#' @export
#' @return \code{sgspls} returns an object of class \code{"sgspls"}, a list that
#' contains the following components: 
#' \item{weights}{a list containing the X and Y pls weights.} 
#' \item{scores}{a list containing the X and Y pls scores.} 
#' \item{names}{a list containing the X and Y names.} 
#' \item{parameters}{a list containing the parameters of the model that was fitted.}
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
#' @seealso Tuning functions \code{\link[sgspls]{calc_pve}},
#' \code{\link[sgspls]{tune_sgspls}}, \code{\link[sgspls]{tune_groups}}. 
#' Model performance and estimation  \code{\link[sgspls]{predict.sgspls}}, \code{\link[sgspls]{perf.sgspls}}, \code{\link[sgspls]{coef.sgspls}} 
#'
#' @examples
#'
#' set.seed(1)
#' n = 50; p = 510; 
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
#' reg_coef <- coef(model, type = "coefficients")
#'
#' # Check model fit
#' cbind(reg_coef$B[ , , 3], beta)
#'
#' \dontrun{
#' cbind(model.sgsplsR$B.hat[,,3], beta)[1:30,]
#' }


sgspls <-
  function(X, Y, ncomp = 2, mode = "regression",  
           keepX = NA, keepY = NA, max.iter = 500,
           tol = 1e-12, scale.x = T, scale.y = F,
           groupX = rep(1,ncol(X)), groupY = rep(1,ncol(Y)),
           subgroupX = rep(1,ncol(X)), subgroupY = rep(1,ncol(Y)),
           indiv_sparsity_x = rep(0,ncomp), subgroup_sparsity_x = rep(0,ncomp),
           indiv_sparsity_y = rep(0,ncomp), subgroup_sparsity_y = rep(0,ncomp),
           ...)
  {
    
    # Check the data #
    X = as.matrix(X)
    Y = as.matrix(Y)
    p = ncol(X); q = ncol(Y); n = nrow(X);
    
    if (any(!is.finite(Y))) stop("\ninfinite, NA or NaN values in X or Y")
    
    if (is.null(colnames(X))) colnames(X) = paste(rep("X",p), 1:p, sep="")
    if (is.null(rownames(X))) rownames(X) = rep(1:n)
    
    if (is.null(colnames(Y))) colnames(Y) = paste(rep("Y",q), 1:q, sep="")
    if (is.null(rownames(Y))) rownames(Y) = rep(1:n)

    # Add group information
    ngX <- length(unique(groupX))
    ngY <- length(unique(groupY))
    
    # Fix number of selected groups
    if(is.na(keepY[1])) keepY <- ngY
    if(is.na(keepX[1])) keepX <- ngX
    
    # Repeat arguments to match ncomp length
    keepX <- rep_param(keepX, ncomp)
    keepY <- rep_param(keepY, ncomp)
    
    indiv_sparsity_x <- rep_param(indiv_sparsity_x,ncomp)
    subgroup_sparsity_x <- rep_param(subgroup_sparsity_x,ncomp)
    
    indiv_sparsity_y <- rep_param(indiv_sparsity_y,ncomp)
    subgroup_sparsity_y <- rep_param(subgroup_sparsity_y,ncomp)

    # Set the alphas for the sparsity values
    # There must be positive penalty on the group so that the keepX and keepY parameters work
    alpha1.x <- 1 - rowSums(cbind(matrix(indiv_sparsity_x, ncol = 1), matrix(subgroup_sparsity_x, ncol = 1)))

    if(any(alpha1.x == 0)){
      ind <- which(alpha1.x == 0)
      alpha1.x[ind] <- 1e-06
      for (i in ind){
        if(indiv_sparsity_x[ind] < subgroup_sparsity_x[ind]) {
          subgroup_sparsity_x[ind] = subgroup_sparsity_x[ind] - 1e-06
        } else{
          indiv_sparsity_x[ind] = indiv_sparsity_x[ind] - 1e-06
        }
      }
    }
    alpha2.x <- subgroup_sparsity_x
    alpha3.x <- indiv_sparsity_x

    alpha1.y <- 1 - rowSums(cbind(indiv_sparsity_y, subgroup_sparsity_y))
    alpha1.y <- pmax(alpha1.y, 1e-6)

    if(any(alpha1.y == 0)){
      ind <- which(alpha1.x == 0)
      alpha1.x[ind] <- 1e-06
      for (i in ind){
        if(indiv_sparsity_x[ind] < subgroup_sparsity_x[ind]) {
          subgroup_sparsity_x[ind] = subgroup_sparsity_x[ind] - 1e-06
        } else{
          indiv_sparsity_x[ind] = indiv_sparsity_x[ind] - 1e-06
        }
      }
    }

    alpha2.y <- subgroup_sparsity_y
    alpha3.y <- indiv_sparsity_y

    if(any(round(rowSums(cbind(alpha1.x, alpha2.x, alpha3.x)),6) != 1)) stop("Sparsities don't sum to 1 on X")
    if(any(round(rowSums(cbind(alpha1.y, alpha2.y, alpha3.y)),6) != 1)) stop("Sparsities don't sum to 1 on Y")

    ## Specialist method to preform sgspls with lambda and alpha's specified on X
    # additional Args #
    lambda.x = rep(0,ncomp)
    lambda.y = rep(0,ncomp)
    singular_vals <- rep(0, ncomp)

    add.args = list(...)
    if(!is.null(add.args$lambda.x)){
      lambda.x <- add.args$lambda.x
      alpha1.x <- add.args$alpha1.x
      alpha2.x <- add.args$alpha2.x
      alpha3.x <- add.args$alpha3.x
    }

    # Scaleing similar to the pls package #
    Y_h <- scale(Y, center = T, scale = scale.y)
    X_h <- scale(X, center = T, scale = scale.x)

    x_scores <- matrix(nrow=n, ncol = ncomp) # X variate
    y_scores <- matrix(nrow=n, ncol = ncomp) # Y variate

    x_weights <- matrix(nrow=p, ncol = ncomp) 
    y_weights <- matrix(nrow=q, ncol = ncomp) 

    # Algorithm #
    for(h in 1:ncomp){
      M = t(X_h)%*%Y_h
      svd.M <- try(svd(M,nu = 1,nv = 1))

      if(class(svd.M) == "try-error"){
        irlba::irlba(M,nu = 1,nv = 1)
      }

      if (svd.M$d[1] < .Machine$double.eps) {
        cat("Too many PLS components, stoped at ",h-1," components",svd.M$d[1])
        ncomp = h-1;
        scores = scores[,1:(h-1)]; weights = weights[,1:(h-1)]; loadings = loadings[,1:(h-1)]
        path.coeff = path.coeff[1:(h-1)]; comps = comps[1:(h-1)]
        break
      }

      # Calaulate the weights
      weights <- cal_weights(M = M,svd.M = svd.M,keepX = keepX[h],keepY = keepY[h],groupX = groupX,
                                                groupY = groupY,subgroupX = subgroupX,subgroupY = subgroupY,alpha1.x = alpha1.x[h],
                                                alpha2.x = alpha2.x[h],alpha1.y = alpha1.y[h],alpha2.y = alpha2.y[h],tol, max.iter,
                                                lambda.x = lambda.x[h], lambda.y = lambda.y[h])
      x_weights[,h] <- weights$x
      y_weights[,h] <- weights$y

      # calculate loadings, scores and deflate X & Y
      x_scores[,h] <- x_score <- X_h %*% weights$x
      y_scores[,h]  <- y_score <- Y_h %*% weights$y
      
      # Singular value type value
      singular_vals[h] <- weights$d
      
      ### Deflation step
      proj_x <- diag(1,n,n) - tcrossprod(x_score)/drop(crossprod(x_score))
      X_h <- proj_x %*% X_h
      
      if (mode=="regression") {
        Y_h <- proj_x %*% Y_h
        }
      else 
        {
          proj_y <- diag(1,n,n) - tcrossprod(y_score)/drop(crossprod(y_score))
          Y_h <- proj_y %*% Y_h
        }
    }

    cl = match.call()

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
                       indiv_sparsity_x = indiv_sparsity_x, subgroup_sparsity_x = subgroup_sparsity_x,
                       indiv_sparsity_y = indiv_sparsity_y, subgroup_sparsity_y = subgroup_sparsity_y,
                       tol = tol, scale.x = scale.x, scale.y = scale.y,
                       X.residuals = X_h, Y.residuals = Y_h)

    result <- list(call = cl, weights = list(X = x_weights, Y = y_weights), scores = list(X = x_scores, Y = y_scores),
                   names = list(X = colnames(X),Y = colnames(Y), indiv = rownames(X)), singular_vals = singular_vals, parameters=parameters)
    class(result) = c("sgspls")
    return(invisible(result))
  }
