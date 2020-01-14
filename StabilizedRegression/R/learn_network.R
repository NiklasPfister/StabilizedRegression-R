##' Learn a network model for a collection of variables.
##'
##' Uses StabilizedRegression, Lasso or correlation to construct a
##' node-wise network between all variables in X.
##' @title Learn network model
##' @param X data matrix. Numeric matrix of size n times d, where
##'   columns correspond to individual variables.
##' @param A stabilizing variable. Numeric vector of length n which
##'   can be interpreted as a factor.
##' @param method specifies which method to use. "SR" for Stabilized
##'   Regression (both standard and predictive version), "SRstab" for
##'   only the standard version of SR, "SRpred" for only the
##'   predictive version of SR, "OLS" for linear OLS regression,
##'   "lasso" for Lasso and "correlation" for correlation test.
##' @param resampling_method specifies which resampling method to
##'   use. Should be one of "none", "stability_selection" or
##'   "permutation".
##' @param numB number of resamples to use.
##' @param cutoff tuning parameter used in stability selection to
##'   determine which sets count as selected.
##' @param pars list of additional parameters passed to SR
##'   regression. See \link{StabilizedRegression} for more details.
##' @param verbose 0 for no output, 1 for text output and 2 for text
##'   and diagnostic plots.
##' @param cores number of cores to use in resampling step.
##' 
##' @return A list consisting of the following
##'   elements
##'
##' \item{Amat}{adjacency matrix, where Amat[i,j] is a score (depending on the resampling_method) for the edge from i to j. For "stability_selection" scores correspond to selection probabilities, for "permutation" scores correspond to permutation p-values and for "none" scores correspond to variable importance of the method.}
##' \item{p}{Total number of potential edges which can be used to compute upper bound on false discovery rate (only computed if resampling_method == "stability_selection").}
##' \item{qest}{Average number of selected edges in stability selection, which can be used to compute upper bound on false discovery rate (only computed if resampling_method == "stability_selection").}
##'
##' If method=="SR" result is a list with two entries SRstab and
##' SRpred each consisting of a list of the form described above.
##' 
##' @export
##'
##' @import stats glmnet parallel
##'
##' @author Niklas Pfister
##'
##' @examples
##' ## Example
##' set.seed(1)
##' X1 <- rnorm(200)
##' X2 <- X1 + rnorm(200)
##' X3 <- 0.5 * X1 + X2 + 0.2 * c(rnorm(100), rnorm(100)+20)
##'
##' X <- cbind(X1, X2, X3)
##' A <- as.factor(rep(c(0, 1), each=100))
##'
##' network <- learn_network(X, A, method="SR", resampling_method="none")
##'
##' print(network[[1]]$Amat)
##' print(network[[2]]$Amat)


learn_network <- function(X, A=NA,
                          method="correlation",
                          resampling_method="stability_selection",
                          numB=100,
                          cutoff=0,
                          pars=list(m=ncol(X),
                                    B=NA,
                                    alpha_stab=0.05,
                                    alpha_pred=0.05,
                                    size_weight="linear",
                                    use_resampling=FALSE,
                                    prescreen_size=nrow(X)-1,
                                    prescreen_type="correlation",
                                    stab_test="exact",
                                    pred_score="mse",
                                    variable_importance="scaled_coefficient"),
                          verbose=0,
                          cores=1){


  d <- ncol(X)
  n <- nrow(X)


  if(method == "SRstab"){
    pars$compute_predictive_model <- FALSE
    adjacency_matrix <- function(ind_a, ind_b){
      Amat <- matrix(0, d, d)
      for(j in 1:d){
        tmp <- StabilizedRegression(X[ind_b, -j, drop=FALSE], X[ind_a,j], A[ind_a],
                                    pars,
                                    verbose)
        Amat[-j, j] <- tmp$variable_importance
      }
      return(Amat)
    }
  }
  else if(method == "SRpred"){
    pars$compute_predictive_model <- TRUE
    adjacency_matrix <- function(ind_a, ind_b){
      Amat <- matrix(0, d, d)
      for(j in 1:d){
        tmp <- StabilizedRegression(X[ind_b, -j, drop=FALSE], X[ind_a,j], A[ind_a],
                                    pars,
                                    verbose)
        Amat[-j, j] <- tmp$variable_importance_pred
      }
      return(Amat)
    }
  }
  else if(method == "SR"){
    pars$compute_predictive_model <- TRUE
    adjacency_matrix <- function(ind_a, ind_b){
      Amat <- matrix(0, d, d)
      Amat_pred <- matrix(0, d, d)
      for(j in 1:d){
        tmp <- StabilizedRegression(X[ind_b, -j, drop=FALSE], X[ind_a,j], A[ind_a],
                                    pars,
                                    verbose)
        Amat[-j, j] <- tmp$variable_importance
        Amat_pred[-j, j] <- tmp$variable_importance_pred
      }
      return(list(Amat=Amat,
                  Amat_pred=Amat_pred))
    }
  }
  else if(method == "correlation"){
    adjacency_matrix <- function(ind_a, ind_b){
      Amat <- matrix(1, d, d)
      combs <- combn(1:d, 2, simplify=FALSE)
      for(i in 1:length(combs)){
        Amat[combs[[i]][1], combs[[i]][2]] <- cor.test(X[ind_a, combs[[i]][1]],
                                                       X[ind_b, combs[[i]][2]])$p.value
        Amat[combs[[i]][2], combs[[i]][1]] <- Amat[combs[[i]][1], combs[[i]][2]]
      }
      return(1-Amat)
    }
  }
  else if(method == "OLS"){
    adjacency_matrix <- function(ind_a, ind_b){
      Amat <- matrix(1, d, d)
      for(j in 1:d){
        Amat[-j, j] <- unname(
          summary(lm(X[ind_a, j] ~ X[ind_b, -j, drop=FALSE]))$coefficients[,"Pr(>|t|)"][-1])
      }
      return(1-Amat)
    }
  }
  else if(method == "lasso"){
    adjacency_matrix <- function(ind_a, ind_b){
      Amat <- matrix(1, d, d)
      for(j in 1:d){
        fit <- glmnet(X[ind_b, -j, drop=FALSE], X[ind_a, j])
        sel_matrix <- (fit$beta != 0)
        first_entrance <- apply(sel_matrix, 1, which.max)
        first_entrance[which(apply(sel_matrix, 1, sum) == 0)] <- Inf
        Amat[-j, j] <- 1/first_entrance
      }
      return(Amat)
    }      
  }
  else{
    stop("Specified method does not exist.")
  }

  ## Learn network
  if(method == "SR"){
    # Apply resampling
    if(resampling_method == "stability_selection"){
      Alist <- mclapply(1:numB,
                        function(x){
                          resample_ind <- sample(1:n, floor(n/2), replace=TRUE)
                          tmpmat <- adjacency_matrix(resample_ind, resample_ind) 
                          mat1 <- tmpmat$Amat_pred
                          mat2 <- tmpmat$Amat
                          return(list(mat1, mat2))
                        }, mc.cores=cores)
      Alist <- lapply(Alist, function(x) list(x[[1]] > cutoff,
                                              x[[2]] > cutoff))
      Amat_pred <- Reduce('+', lapply(Alist, function(x) x[[1]]))/numB
      Amat <- Reduce('+', lapply(Alist, function(x) x[[2]]))/numB
      p <- choose(d, 2)*2
      qest_pred <- mean(sapply(lapply(Alist, function(x) x[[1]]), function(x) sum(x)))
      qest <- mean(sapply(lapply(Alist, function(x) x[[2]]), function(x) sum(x)))
    }
    else if(resampling_method == "permutation"){
      Astat <- adjacency_matrix(1:n, 1:n)
      Alist <- mclapply(1:numB,
                        function(x){
                          resample_ind <- sample(1:n)
                          tmpmat <- adjacency_matrix(1:n, resample_ind)
                          mat1 <- tmpmat$Amat_pred
                          mat2 <- tmpmat$Amat
                          return(list(mat1, mat2))
                        }, mc.cores=cores)
      Amat_pred <- (Reduce('+',
                           lapply(Alist,
                                  function(x) x[[1]] > Astat$Amat_pred))+1)/(numB+1)
      Amat <- (Reduce('+',
                      lapply(Alist,
                             function(x) x[[2]] > Astat$Amat))+1)/(numB+1)
      diag(Amat) <- 0
      diag(Amat_pred) <- 0
      p <- NA
      qest_pred <- NA
      qest <- NA
    }
    else if(resampling_method == "none"){
      tmp <- adjacency_matrix(1:n, 1:n)
      Amat <- tmp$Amat
      Amat_pred <- tmp$Amat_pred
      p <- NA
      qest_pred <- NA
      qest <- NA
    }
    else{
      stop("Specified resampling_method does not exist.")
    }
    # Collect results
    reslist <- list(SRpred=list(Amat=Amat_pred,
                                p=p,
                                qest=qest_pred),
                    SRstab=list(Amat=Amat,
                                p=p,
                                qest=qest))
  }
  else{
    # Apply resampling
    if(resampling_method == "stability_selection"){
      Alist <- mclapply(1:numB,
                        function(x){
                          resample_ind <- sample(1:n, floor(n/2), replace=TRUE)
                          return(adjacency_matrix(resample_ind, resample_ind))
                        }, mc.cores=cores)
      Alist <- lapply(Alist, function(x) x > cutoff)
      Amat <- Reduce('+', Alist)/numB
      p <- choose(d, 2)*2
      qest <- mean(sapply(Alist, function(x) sum(x)))      
    }
    else if(resampling_method == "permutation"){
      Astat <- adjacency_matrix(1:n, 1:n)
      Alist <- mclapply(1:numB,
                        function(x){
                          resample_ind <- sample(1:n)
                          return(adjacency_matrix(1:n, resample_ind))
                        }, mc.cores=cores)
      Amat <- (Reduce('+', lapply(Alist, function(x) x > Astat))+1)/(numB+1)
      diag(Amat) <- 0
      p <- NA
      qest <- NA
    }
    else if(resampling_method == "none"){
      Amat <- adjacency_matrix(1:n, 1:n)
      p <- NA
      qest <- NA
    }
    else{
      stop("Specified resampling_method does not exist.")
    }
    # Collect results
    reslist <- list(Amat=Amat,
                    p=p,
                    qest=qest)
  }
  

  return(reslist)
}
