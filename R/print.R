# print fit
print.sgspls <-
  function( x, ... )
  {
    ncomp <- x$ncomp
    mode <- x$mode
    keepX <- x$keepX
    keepY <- x$keepY
    p <- ncol(x$X)
    q <- ncol(x$Y)
    groupX <- x$parameters$groupX
    subgroupX <- x$parameters$subgroupX
    groupY <- x$parameters$groupY
    subgroupY <- x$parameters$subgroupY
    indiv.sparsity.x <- x$parameters$indiv.sparsity.x
    indiv.sparsity.y <- x$parameters$indiv.sparsity.y
    subgroup.sparsity.x <- x$parameters$subgroup.sparsity.x
    subgroup.sparsity.y <- x$parameters$subgroup.sparsity.y
    
    # cat( "\nSparse Group Sub-group Partial Least Squares \n" )
    # cat( "----\n\n")
    cat( paste0("\n sgspls with a '",mode,"' mode with ",ncomp," components. \n\n"))
    # cat(" You entered data X of dimensions:", n, p, "\n")
    # cat(" You entered data Y of dimensions:", n, q, "\n\n")
    
    cat( paste0(" Selected ",sum(rowSums(abs(x$loadings$X)>0)>0)," variables among ",p," variables on the X block. \n"))
    cat( " Selected variables per component:",colSums(abs(x$loadings$X)>0),"\n")
    cat( " Selected Groups:",unique(unlist(apply(x$loadings$X,2,function(x)groupX[which(abs(x)>0)]))),"\n\n")
    cat( paste0(" Selected ",sum(rowSums(abs(x$loadings$Y)>0)>0)," variables among ",q," variables on the Y block.\n") )
    cat( " Selected variables per component:",colSums(abs(x$loadings$Y)>0),"\n")
    cat( " Selected Groups:",unique(unlist(apply(x$loadings$Y,2,function(x)groupY[which(abs(x)>0)]))),"\n\n")
    
    cat( " Model Parameters:\n\n") 
    if(keepX[1]==length(unique(groupX)) && indiv.sparsity.x[1] -subgroup.sparsity.x[1] == 0){
      cat(" No X block penalisation")
    } else{
    cat( "  keepX ", " group ", " subgroup ", " individual \n" )
    for(i in 1:ncomp){
      cat("   ",keepX[i],"    ", 1-indiv.sparsity.x[i] -subgroup.sparsity.x[i], "    ",subgroup.sparsity.x[i], "      ",indiv.sparsity.x[i]," \n" )
    }
    }
    cat("\n")
    if(is.null(keepY))    {keepY <- length(unique(groupY))}
    if(keepY[1]==length(unique(groupY)) && indiv.sparsity.y[1] -subgroup.sparsity.y[1] == 0){
      cat(" No Y block penalisation")
    } else{
    
    cat( "  keepY ", " group ", " subgroup ", " individual \n" )
    for(i in 1:ncomp){
      cat("   ",keepY[i],"    ",1-indiv.sparsity.y[i] -subgroup.sparsity.y[i], "    ",subgroup.sparsity.y[i], "      ",indiv.sparsity.y[i]," \n" )
    }
    }
    cat("\n\n")
    
    cat(" Available components: \n", 
        "-------------------- \n")
    
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
  }

# print Cross validation
print.cv.sgspls <-
  function( x, ... )
  {
    cat("\n ========================================== \n")
    cat("\n",x$folds,"fold Cross Validation for sgspls \n\n")
    
    cat( " Optimal Parameters:\n\n") 
    cat( "  keepX ", " group ", " subgroup ", " individual \n" )
    cat("   ",x$keepX,"    ", 1 - x$indiv.sparsity -x$subgroup.sparsity, "    ",x$subgroup.sparsity, "      ",x$indiv.sparsity," \n" )
    cat("\n")
    cat( " Optimal MSEP:",x$min.cv,"\n\n")
    cat(" Available plotting: see plot(object) \n")
    cat("\n ========================================== \n")
  }