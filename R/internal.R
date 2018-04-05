#' Find selected Modules, Genes and Times
#'
#' This function finds the selected groups, subgroups and individual predictors
#' from a sgspls method.
#'
#' @param model object of class inhereting from "sgspls".
#' @param module A p-vector indicating group membership for each covariate in
#'   the X-block
#' @param gene A p-vector indicating gene membership for each covariate in
#'   the X-block
#' @param time A p-vector indicating time membership for each covariate in
#'   the X-block
#'
#' @export
#' @return Returns a list with the following selected parameter information:
#'
#' \item{select.table.X}{A table detailing the number of times each gene has been selected at each timepoint and the number of consistently selected genes.}
#' \item{summary.table}{A table sumarising the number of modules, genes, timepoints and covariates selected.}
#' \item{tab.gene.X}{Lists the number of timepoint that each gene in a given module and component occurs.}
#' \item{tab.gene.time.X}{Table of the selected genes against time points that they occur.}
#' \item{consistent.genes.X}{Returns the genes that occur in every time point.}
#' \item{select.gene.X}{Returns the genes that are selected at least once for a given component.}
#' \item{select.gene.X.total}{Returns the genes that are selected at least once across any component.}
#' \item{selected.table.gene.X}{Returns the total number of genes selected at each component.}
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
#' # Beta contains 1 active module with 3 active genes
#' # index for modules, genes and times
#' mod_index <- groupX
#' gene_index <- subgroupX
#' time_index <- rep(rep(1:5,6), 17)
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
#' # See the estimated regression coefficient
#' cbind(reg_coef$B[,,3], beta, mod_index, gene_index, time_index)[1:30,]
#' selectedVar <- select.sgspls(model, module = mod_index, gene = gene_index, time = time_index)
#' # show number of selected genes for component 1
#' selectedVar$select.table.X$comp1
#' # show number of modules, genes, times and total variables selected
#' selectedVar$summary.table
#' # Show when genes were selected from given module
#'

select.sgspls <- function (model, module, gene, time) {
  ##-- Selected variable metrics --##
  module <- as.factor(paste0("M",module))
  gene <- paste0("G", gene)
  
  ncomp <- model$parameters$ncomp
  
  result <- vector("list", length = ncomp)
  res.select <- vector("list", length = ncomp)
  tab.gene <- vector("list", length = ncomp)
  tab.gene.time <- vector("list", length = ncomp)
  consistent.genes <- vector("list", length = ncomp)
  names(result) <- names(res.select) <- names(tab.gene) <- names(tab.gene.time) <- names(consistent.genes) <- paste0("comp",1:ncomp)
  
  result.summary <- matrix(0, nrow = ncomp, ncol = 4)
  colnames(result.summary) <- c("# Modules", "# Genes", "# Times", " Total ")
  rownames(result.summary) <- paste0("comp ",1:ncomp)
  
  n.modules <- length(unique(module))
  n.times <- length(unique(time))
  size.group <- diag(table(module[time==1], module[time==1]))
  
  tab.gene.select <- matrix(0, nrow = n.modules, ncol = 1 + ncomp)
  tab.gene.select[,1] <- size.group
  colnames(tab.gene.select) <- c("size.group", paste0("comp ",1:ncomp))
  rownames(tab.gene.select) <- unique(module)
  
  for (h in 1:ncomp) {
    set.ind.zero <- which(abs(model$weights$X[, h]) > 0)
    # summary stats
    ns.mod <- length(unique(module[set.ind.zero]))
    ns.gene <- length(unique(gene[set.ind.zero]))
    ns.time <- length(unique(time[set.ind.zero]))
    ns <- length(set.ind.zero)
    result.summary[h,] <- c(ns.mod, ns.gene, ns.time, ns)
    
    names(set.ind.zero) <- model$names$X[set.ind.zero]
    res.select[[h]] <- set.ind.zero
    res.gene.time <- res.gene <- vector("list", length = n.modules)
    names(res.gene.time )<-names(res.gene) <- unique(module)
    consistent <- res <- NULL
    for (mod in 1:n.modules) {
      
      inds <- intersect(set.ind.zero,which(as.numeric(module)==mod))
      res.gene.time[[mod]] <- smods <- table(as.character(gene[inds]), time[inds])
      res.gene[[mod]] <- sort(rowSums(smods),decreasing = T)
      consistent <- c(consistent,length(which(rowSums(smods)==n.times)))
      temp <- rep(0,n.times)
      temp[as.numeric(names(colSums(smods)))] <- colSums(smods)
      res <- rbind(res, temp)
      
    }
    tab.gene.time[[h]] <- res.gene.time
    tab.gene[[h]] <- res.gene
    colnames(res) <- paste0("T",1:ncol(res))
    result[[h]] <- cbind(size.group, consistent, res)
    tab.gene.select[,h+1] <- unlist(lapply(res.gene, function(x) length(x)))
    consistent.genes[[h]] <- unlist(lapply(res.gene, function(x) names(which(x==6))))
    cons.names <- NULL
    #     for(i in 1:n.modules){
    #       cons.names <- c(cons.names, rep(unique(as.character(module))[i], consistent[i]))
    #     }
    #     names(consistent.genes[[h]]) <- cons.names
  }
  select.x <- lapply(res.select, function(x) unique(as.character(gene[x])))
  ind.total.x <- sort(unique(as.character(unlist(select.x))))
  return(list(select.table.X = result, summary.table = result.summary))
}


# Duplicate penalty parameter to match the number of components
rep_param <- function(arg, ncomp){
  if(ncomp != length(arg)){
    arg = rep(arg, length.out = ncomp)
  }
  return(arg)
}

#' Compute scaled vector multiplicaiton 
#'
#' Return multiplation of a matrix by scaled vector
#' return 0 vector if scale is zero (avoid division by zero).
#'
#' @param M n x p matrix of numeric variables
#' @param u p vector of numeric variables
#'
#' @return Output will be a numeric matrix or vector.
#'

scale_vec <- function(u, scale = 1){
  p <- length(u)
  scale_factor <- drop(crossprod(u)^scale)
  
  res <- if (scale_factor > 10e-8)  u/scale_factor else matrix(rep(0, p), nrow = p)
  return(res)
}

# scale and center a matrix
pls.scale <- function(X, scale.x){
  meanx <- apply(X, MARGIN = 2, mean)
  normx <- rep(1,ncol(X))
  if(scale.x) {
    normx <- apply(X, 2, sd)
    if (any(normx < .Machine$double.eps)) {
      stop("Some of the columns of the predictor matrix have zero variance.")
    }
  }
  X = scale(X,meanx,normx)
  attr(X, "scaled:scaled") <- normx
  return(X)
}


#' Return the nonzero columns and rows
#'
#' Quick function for finding the nonzero columns of a matrix.
#'
#' @param Z A matrix of variables, where some columns or rows are entierly zero.
#' @param thresh A threshold value to round values off at.
#' @param zero.col A vector of columns to take for the removal of the matrix.
#'
#' @examples
#' a <- matrix(c(rep(0,3),1,0,1,rep(0,3)),3,3)
#' a
#' nonZero(a) #removes all rows and columns with nonzero elements
#' nonZero(a,zero.col = 2) # removes only nonzero rows and columns from row 2.

nonZero <- function(Z, thresh = 1e-06,zero.col = NULL){
  if(is.null(zero.col)) {
    Z.nz = if(is.null(dim(Z))) Z[which(abs(Z)>thresh)] else Z[rowSums(abs(Z))>thresh,colSums(abs(Z))>thresh]
  }
  else {
    Z.nz = if(length(zero.col)==1) Z[which(abs(Z[,zero.col])>thresh),] else Z[rowSums(abs(Z[,zero.col]))>thresh,colSums(abs(Z[,zero.col]))>thresh]
  }
  return(Z.nz)
}

#' Plot the cross validation curves of a sgspls.cv object
#'
#' Produces a plot of the cross validation curves for a "sgspls.cv" object.
#'
#' @param obj fitted sgspls.cv object
#'
#' @seealso sgspls.tune, perf
#'
#'
plot.cv.sgspls <- function(obj, verbose = T){
  if(verbose){
    library(knitr)
    cat("\n ========================================== \n")
    cat("\n",obj$folds,"fold Cross Validation for sgspls \n\n")
    nalpha <- length(obj$tuning_sparsities[,1])

    cat( " Tuning Parameters:")
    table <- data.frame(Legend = paste("alpha",1:nalpha,sep = "_"),obj$tuning_sparsities)
    print(kable(table, format = "pandoc"))
    cat("\n ========================================== \n")
  }

  legendVal <- parse(text = paste("alpha[",1:nalpha,"]",sep=""))
  group_seq <- obj$group_seq
  cv.scores <- obj$results_tuning[,4]
  cv.scores <- matrix(cv.scores, nrow = length(group_seq), ncol = nalpha)
  minalpha <- which(cv.scores == min(cv.scores), arr.ind = T)[2]
  alphawidth <- rep(1,nalpha)
  alphawidth[minalpha] <- 3

  # get nice Y-range
  placement = max(cv.scores[1,]) > max(cv.scores[nrow(cv.scores),])

  if(placement) {
    maxy<- max(cv.scores[nrow(cv.scores),]) + max(0.8*(max(cv.scores[1,])-max(cv.scores[nrow(cv.scores),])), 0.2)
    miny <- min(cv.scores)
  } else{
    miny <- min(cv.scores[nrow(cv.scores),]) - max(min(cv.scores[nrow(cv.scores),]) - min(cv.scores[1,]), 0.2)
    maxy <- max(c(cv.scores[nrow(cv.scores),],cv.scores[1,]))
  }

  placement <- if(placement) "topright" else "bottomright"
  matplot(group_seq,cv.scores,type="l",ylab="MSEP",ylim = c(miny,maxy),
          pch=1,col=rainbow(nalpha),lty=1:nalpha,lwd=1.5*alphawidth,xlab="Number of groups")

  # How many colums to use <10 = 1, <20 = 2, <30 = 3, otherwise 4.
  ncols <- which(table(cut(nalpha, breaks = c(0,10,20,30,100)))>0)

  legend(placement,legendVal,ncol=ncols,col = rainbow(nalpha),lty=1:nalpha,lwd=alphawidth)

}

plot.sgspls <- function(obj, verbose = T){
  if(obj$parameters$mode == "regression"){
    pve <- calc_pve(obj)
    plot(drop(pve), xlab = "Number of components", ylab = "Percentage Variance Explained", type = "l")
  } else{
    singular_values <- obj$singular_vals/obj$singular_vals[1] 
    plot(singular_values, xlab = "Number of components", ylab = "Normalised singular values", type = "l")
  }
}


sim_regression <- function(n = 100, q = 1, coef_subgroup = rep(1, 5), nonzero_subgroups = 3, 
                           zero_subgroups = 2, nonzero_groups = 5, zero_groups = 5, sigma = 0.5) {
  
  generate_group <- function(coef_g, nonzero_groups, zero_groups){
    
    n_groups <- nonzero_groups + zero_groups
    coef <- matrix(0, nrow = length(coef_g)*n_groups)
    
    groupInd <- ceiling(1:length(coef) / length(coef_g))
    
    sel_groups <- sample(1:n_groups, replace = F, size = nonzero_groups)
    coef[which(groupInd %in% sel_groups)] <- coef_g
    return(coef)
  }
  
  B <- NULL
  for(resp in 1:q){
    coef_group <- generate_group(coef_subgroup, nonzero_subgroups, zero_subgroups)
    coef <- generate_group(coef_group, nonzero_groups, zero_groups)
    B <- cbind(B, coef)
  }
  p <- length(coef)
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  E <- MASS::mvrnorm(n, rep(0,q), Sigma = diag(sigma, q, q))
  Y <- X%*%B + E
  
  groupX <- ceiling(1:p / length(coef_group))
  subgroupX <- ceiling(1:p / length(coef_subgroup))
  
  return(list(X=X, Y=Y, B = B, groupX = groupX, subgroupX = subgroupX))
}


