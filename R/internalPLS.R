#' Internal sgspls functions
#'
#' Internal sgspls functions
#'
#' @details These are not intended for use by the user.
#'   \code{softhresh(Mx,lambda)} provides soft thresholding on vector Mx with
#'   parameter lambda, \code{checkalphas} checks if the alpha values entered are
#'   appropriate.
#'
#' @author Matthew Sutton
#' @name sgspls-internal
NULL

#' @rdname sgspls-internal
cal_weights <- function(M,svd.M,keepX=NA,keepY=NA,groupX=NA,groupY=NA,subgroupX=NA,subgroupY=NA,
                                          alpha1.x=0,alpha2.x=0,alpha1.y=0,alpha2.y=0,tol=1e-06,max.iter=500,
                                          lambda.x=0, lambda.y=0){

  u <- uold <- matrix(svd.M$u)
  v <- vold <- matrix(svd.M$v)

  # Change group and subgroups for speed if no penalisation is given
  if(alpha1.x == 0) groupX = rep(1,length(u))
  if(alpha2.x == 0) subgroupX = rep(1,length(u))
  if(alpha1.y == 0) groupY = rep(1,length(v))
  if(alpha2.y == 0) subgroupY = rep(1,length(v))

  ngx = length(unique(groupX)); ngy = length(unique(groupY))
  gx = ngx - keepX; gy = ngy - keepY;
  GX = GY = NULL
  Mv = M %*% vold
  Mu = t(M) %*% uold

  # A good initiall guess for the value setting the weights to 0
  lamb.upper.x = 2*max(abs(Mv))/abs(1-alpha1.x - alpha2.x) + 1
  lamb.upper.x = min(lamb.upper.x, 10^5)
  
  lamb.upper.y = 2*max(abs(Mu))/abs(1-alpha1.y - alpha2.y) + 1
  lamb.upper.y = min(lamb.upper.x, 10^5)
  
  gXunique = unique(groupX)
  gYunique = unique(groupY)

  # determine if penalisation is required for X and Y blocks
  penaliseX = (alpha1.x < 1 || (keepX != ngx) || lambda.x >0)
  penaliseY = (alpha1.y < 1 || (keepY != ngy) || lambda.y >0)


  lambda1 <- lambda.x
  lambda2 <- lambda.y
  ittr = 0

  repeat{
    ittr = ittr + 1
    # Update X-weights (u)
    GX = NULL
    Mv = M%*%vold
    # find lambda corresponding to X group number
    if(penaliseX){
      for(k in 1:ngx){
        ind = which(groupX == gXunique[k])
        GX = c(GX, uniroot(lambdazerosubgroup, lower=0, upper=lamb.upper.x, Mx=Mv[ind], subgroupX=subgroupX[ind], alpha1=alpha1.x,alpha2=alpha2.x,tol = tol)$root)
      }
      lambda.x <- lambda1 <- max(sort(GX)[gx], 0, na.rm = T)
    }
    # Updata u (calls Cpp funciton)
    u = if (penaliseX) updataU(Mx = Mv,groupX,subgroupX,lambda1,alpha1.x,alpha2.x ) else Mv
    u = scale_vec(u, scale = 0.5)

    # Update Y-weights (v)
    GY = NULL
    Mu = t(M) %*% uold

    # find lambda corresponding to Y group number
    if(penaliseY){
      for(ell in 1:ngy){
        ind = which(groupY == gYunique[ell])
        GY = c(GY, uniroot(lambdazerosubgroup, lower=0, upper=lamb.upper.y, Mx=Mu[ind], subgroupX=subgroupY[ind], alpha1=alpha1.y,alpha2=alpha2.y, tol = tol)$root)
      }
      lambda.y <- lambda2 <- max(sort(GY)[gy], 0, na.rm = T)
    }
    # Updata v (calls Cpp funciton)
    v = if (penaliseY) updataU(Mx = Mu,groupY,subgroupY,lambda2,alpha1.y,alpha2.y ) else Mu
    v = scale_vec(v, scale = 0.5)
    # Check convergence
    if (crossprod(u-uold)+crossprod(v-vold) < tol || ittr >= max.iter || all(u==0)&all(v==0)) {break}
    uold = u; vold = v
  }
  return(list(x=u,y=v,lambda.x = lambda1, lambda.y = lambda2))
}


#' @rdname sgspls-internal
checkalphas <- function(alpha)  return(pmin(pmax(alpha,0),0.9999999)) # alpha kept between 0 and 0.99999


#' Get the alpha values corresponding to the required sparsity
#'
#' Checks the sparsity levels are correct and returns 
#' corresponding alpha values. Alpha values should add to one.
#' The values correspond to the percentage of penalty put on each 
#' of the penalised terms in the sgspls penalty.
#' Current implementation forces some penalisation on the group component.
#'
#' @param  indiv_sparsity matrix of values between 0 and 1 for the individual sparsity penalty per component
#' @param  subgroup_sparsity matrix of values between 0 and 1 for the subgroup sparsity penalty per component
#' 
#' @export
#' 
#' @examples
#'
#' # Return the alphas used for the required sparsities
#' get_alphas(indiv_sparsity = c(0.4,0), subgroup_sparsity = c(0.2,1))
#' 


get_alphas <- function(indiv_sparsity, subgroup_sparsity){
  
  alpha1 <- 1 - rowSums(cbind(indiv_sparsity, subgroup_sparsity))
  alpha2 <- subgroup_sparsity
  alpha3 <- indiv_sparsity
  
  #-- Find which components have no group penalty --#
  nogroup <- which(alpha1 < 1e-6)
  
  #-- Set group penalty --#
  for(i in nogroup)
    {
    alpha1[i] <-  1e-6
    alpha2[i] <- alpha2[i] - alpha2[i]*1e-6
    alpha3[i] <- alpha3[i] - alpha3[i]*1e-6
  }
  
  res <- list(alpha1 = alpha1, 
              alpha2 = alpha2, 
              alpha3 = alpha3)
  return(res)
}


