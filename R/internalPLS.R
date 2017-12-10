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
group.sparse.subgroup.penalty <- function(M,svd.M,keepX=NA,keepY=NA,groupX=NA,groupY=NA,subgroupX=NA,subgroupY=NA,
                                          alpha1.x=0,alpha2.x=0,alpha1.y=0,alpha2.y=0,tol=1e-06,max.iter=500,
                                          lambda.x=0, lambda.y=0, newtry = T){
  #
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
  Mv = M%*%vold
  Mu = t(M)%*%uold

  # A good initiall guess for the value setting the weights to 0
  lamb.upper.x = min( 4*max(abs(Mv)+100)/abs(1-alpha1.x - alpha2.x), 10^5)
  lamb.upper.y = min( 4*max(abs(Mu)+100)/abs(1-alpha1.y - alpha2.y), 10^5)

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
    if(penaliseX & lambda.x == 0){
      for(k in 1:ngx){
        ind = which(groupX == gXunique[k])
        GX = c(GX, uniroot(lambdazerosubgroup, lower=0, upper=lamb.upper.x, Mx=Mv[ind], subgroupX=subgroupX[ind], alpha1=alpha1.x,alpha2=alpha2.x,tol = tol)$root)
      }
      lu <- min(sort(GX)[gx+1], max(GX), na.rm = T)
      ll <- max(sort(GX)[gx], 0, na.rm = T)

      if(newtry){
        lambda.x <- lambda1 <- ll #+ abs(lu-ll)*(1-alpha1.x) #ll
        #cat("a",alpha1.x,"l",lambda.x,"\n")
      } else{
        lambda1 <- ll
      }
      #cat(ll,lu,lambda1,"\n")
    }
    # Updata u (calls Cpp funciton)
    u = if (penaliseX) updataU(Mx = Mv,groupX,subgroupX,lambda1,alpha1.x,alpha2.x ) else Mv
    #cat("u",u,"\n")
    #cat("alpha3", 1 - alpha1.x - alpha2.x)
    u = u/e.norm(u)

    # Update Y-weights (v)
    GY = NULL
    Mu = t(M)%*%uold

    # find lambda corresponding to Y group number
    if(penaliseY & lambda.y == 0){
      for(ell in 1:ngy){
        ind = which(groupY == gYunique[ell])
        GY = c(GY, uniroot(lambdazerosubgroup, lower=0, upper=lamb.upper.y, Mx=Mu[ind], subgroupX=subgroupY[ind], alpha1=alpha1.y,alpha2=alpha2.y, tol = tol)$root)
      }
      lu <- min(sort(GY)[gy+1], max(GY), na.rm = T)
      ll <- max(sort(GY)[gy], 0, na.rm = T)

      lambda1 = ll
    }
    # Updata v (calls Cpp funciton)
    v = if (penaliseY) updataU(Mx = Mu,groupY,subgroupY,lambda2,alpha1.y,alpha2.y ) else Mu
    v = v/e.norm(v)
    # Check convergence
    if (e.norm(u-uold,exact=T)/e.norm(u)+e.norm(v-vold,exact=T)/e.norm(v) < tol || ittr >= max.iter || all(u==0)&all(v==0)) {break}
    uold = u; vold = v
  }
  return(list(u=u,v=v,lambda.x = lambda1, lambda.y = lambda2))
  #return(list(u=u,v=v))
}

#' @rdname sgspls-internal
deflate.pls <- function(X,Y,u,v,mode){
  ### Step d
  xi.h <- X%*% matrix(u,ncol=1)/e.norm(u)^2
  w.h  <- Y%*% matrix(v,ncol=1)/e.norm(v)^2

  ### Step e
  c.h <- t(X)%*%matrix(xi.h,ncol=1)/e.norm(xi.h)^2

  d.rh <- t(Y)%*%matrix(xi.h,ncol=1)/(sum(xi.h*xi.h) + (sum(xi.h)==0)*1e-16)

  d.h <- t(Y)%*%matrix(w.h,ncol=1)/(sum(w.h*w.h) + (sum(w.h)==0)*1e-16)

  ###Step f and g
  X_h <- X - xi.h%*%t(c.h)
  if (mode=="regression") Y_h <- Y - xi.h%*%t(d.rh) else Y_h <- Y - w.h%*%t(d.h)

  res <- list(X_h=X_h,Y_h=Y_h,c=c.h,d=d.rh,e=d.h)
  return(res)
}

#' @rdname sgspls-internal
checkalphas <- function(alpha)  return(pmin(pmax(alpha,0),0.9999999)) # alpha kept between 0 and 0.99999
