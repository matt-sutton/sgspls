print.sgspls <-
  function( x, ... )
  {
    parameters <- x$parameters
    ncomp <- parameters$ncomp
    mode <- parameters$mode
    keepX <- parameters$keepX
    keepY <- parameters$keepY
    p <- ncol(parameters$X)
    q <- ncol(parameters$Y)
    
    groupX <- parameters$groupX
    subgroupX <- parameters$subgroupX
    groupY <- parameters$groupY
    subgroupY <- parameters$subgroupY
    
    indiv_sparsity_x <- parameters$indiv_sparsity_x
    indiv_sparsity_y <- parameters$indiv_sparsity_y
    subgroup_sparsity_x <- parameters$subgroup_sparsity_x
    subgroup_sparsity_y <- parameters$subgroup_sparsity_y
    
    # cat( "\nSparse Group Sub-group Partial Least Squares \n" )
    # cat( "----\n\n")
    cat( paste0("\n sgspls with a '",mode,"' mode with ",ncomp," components. \n\n"))
    # cat(" You entered data X of dimensions:", n, p, "\n")
    # cat(" You entered data Y of dimensions:", n, q, "\n\n")
    
    cat( paste0(" Selected ",sum(rowSums(abs(x$weights$X)>0)>0)," variables among ",p," variables on the X block. \n"))
    cat( " Selected variables per component:",colSums(abs(x$weights$X)>0),"\n")
    cat( " Selected Groups:",unique(unlist(apply(x$weights$X,2,function(x)groupX[which(abs(x)>0)]))),"\n\n")
    
    cat( paste0(" Selected ",sum(rowSums(abs(x$weights$Y)>0)>0)," variables among ",q," variables on the Y block.\n") )
    cat( " Selected variables per component:",colSums(abs(x$weights$Y)>0),"\n")
    cat( " Selected Groups:",unique(unlist(apply(x$weights$Y,2,function(x)groupY[which(abs(x)>0)]))),"\n\n")
    
    cat( " Model Parameters:\n\n") 
    if(keepX[1]==length(unique(groupX)) && indiv_sparsity_x[1] -subgroup_sparsity_x[1] == 0){
      cat(" No X block penalisation")
    } else{
    cat( "  keepX ", " group ", " subgroup ", " individual \n" )
    for(i in 1:ncomp){
      cat("   ",keepX[i],"    ", 1-indiv_sparsity_x[i] -subgroup_sparsity_x[i], "    ",subgroup_sparsity_x[i], "      ",indiv_sparsity_x[i]," \n" )
    }
    }
    cat("\n")
    if(is.null(keepY))    {keepY <- length(unique(groupY))}
    if(keepY[1]==length(unique(groupY)) && indiv_sparsity_y[1] -subgroup_sparsity_y[1] == 0){
      cat(" No Y block penalisation")
    } else{
    
    cat( "  keepY ", " group ", " subgroup ", " individual \n" )
    for(i in 1:ncomp){
      cat("   ",keepY[i],"    ",1-indiv_sparsity_y[i] -subgroup_sparsity_y[i], "    ",subgroup_sparsity_y[i], "      ",indiv_sparsity_y[i]," \n" )
    }
    }
    cat("\n\n")
    
    cat(" Available components: \n", 
        "-------------------- \n")
    
    cat(" weights: see object$weights \n")
    cat(" scores: see object$scores \n")
    cat(" variable names: see object$names \n")
  }

print.cv.sgspls <-
  function( x, ... )
  {
    cat("\n ========================================== \n")
    cat("\n",x$folds,"fold Cross Validation for sgspls \n\n")
    
    cat( " Optimal Parameters:\n\n") 
    cat( "  keepX ", " group ", " subgroup ", " individual \n" )
    cat("   ",x$best[1],"    ", 1 - x$best[2] - x$best[3], "    ",x$best[3], "      ",x$best[2]," \n" )
    cat("\n")
    cat( " Optimal MSEP:",x$min_cv,"\n\n")
    cat(" Available plotting: see plot(object) \n")
    cat("\n ========================================== \n")
  }