###
# Score functions
###

##' @import corpcor


bootstrap_mse <- function(Y, X, A, num, pred_score){
  if(pred_score == "mse"){
    n <- length(Y)
    mse <- function(ind){
      beta <- matrix(coefficients(lm.fit(X[ind,,drop=FALSE], Y[ind])), ncol=1)
      is_na <- is.na(beta)
      return(mean((X[,!is_na, drop=FALSE] %*% beta[!is_na] - Y)^2))
    }
    predmat <- sapply(1:num, function(x) mse(sample(1:n, n, replace=TRUE)))
    return(predmat)
  }
  else if(pred_score == "mse_env"){
    Alist <- lapply(unique(A), function(a) which(A == a))
    predmat <- matrix(NA, num, length(Alist))
    mse_env <- function(ind, Xin, Yin){
      beta <- matrix(coefficients(lm.fit(Xin[ind,,drop=FALSE], Yin[ind])), ncol=1)
      is_na <- is.na(beta)
      return(mean((Xin[,!is_na, drop=FALSE] %*% beta[!is_na] - Yin)^2))
    }
    for(i in 1:length(Alist)){
      Xtmp <- X[Alist[[i]],,drop=FALSE]
      Ytmp <- Y[Alist[[i]]]
      n <- length(Ytmp)
      predmat[,i] <- sapply(1:num, function(x)
        mse_env(sample(1:n, n, replace=TRUE), Xtmp, Ytmp))
    }
    return(apply(predmat, 1, min))
  }
  else if(pred_score == "expvar_env"){
    Alist <- lapply(unique(A), function(a) which(A == a))
    predmat <- matrix(NA, num, length(Alist))
    mse_env <- function(ind, Xin, Yin){
      beta <- matrix(coefficients(lm.fit(Xin[ind,,drop=FALSE], Yin[ind])), ncol=1)
      is_na <- is.na(beta)
      return(mean((Xin[,!is_na, drop=FALSE] %*% beta[!is_na] - Yin)^2)/mean((Yin-mean(Yin))^2))
    }
    for(i in 1:length(Alist)){
      Xtmp <- X[Alist[[i]],,drop=FALSE]
      Ytmp <- Y[Alist[[i]]]
      n <- length(Ytmp)
      predmat[,i] <- sapply(1:num, function(x)
        mse_env(sample(1:n, n, replace=TRUE), Xtmp, Ytmp))
    }
    return(apply(predmat, 1, min))
  }
}



VariableImportance <- function(X, Y, estimator, B=2, classification=FALSE){
  d <- ncol(X)
  n <- nrow(X)

  variable_importance <- rep(0, d)
  if(classification){
    metric <- function(Yhat){
      return(mean(Y == Yhat))
    }
  }
  else{
    varY <- mean((Y - mean(Y))^2)
    metric <- function(Yhat){
      return(1-mean((Y - Yhat)^2)/varY)
    }
  }
  baseline <- metric(predict(estimator, X))
  for(j in 1:d){
    Xtmp <- X
    tmp_imp <- rep(0, B)
    for(i in 1:B){
      Xtmp[,j] <- Xtmp[sample(1:n),j]
      tmp_imp[i] <- baseline-metric(predict(estimator, Xtmp))
    }
    variable_importance[j] <- mean(tmp_imp)
  }  
  return(variable_importance)
}


deconfounding_correlation <- function(Xmat, Y){
  Xmat <- cbind(Y, Xmat)
  n <- nrow(Xmat)
  d <- ncol(Xmat)
  Xmat <- (Xmat-matrix(colMeans(Xmat), n, d, byrow=T))/matrix(apply(Xmat, 2, sd),
                                                              n, d, byrow=T)
  ## subsample if n too large
  B <- 20
  corrmat <- matrix(NA, B, d)
  num_var <- sqrt(d)
  if(num_var < n){
    for(i in 1:B){
      svd_fit <- corpcor::fast.svd(Xmat[sample(1:n, num_var, replace=TRUE),,drop=FALSE])
      corrmat[i,] <- (svd_fit$v %*% t(svd_fit$v))[1,]
    }
    corrvec <- colMeans(corrmat)
  }
  else{
    svd_fit <- corpcor::fast.svd(Xmat)
    corrvec <- (svd_fit$v %*% t(svd_fit$v))[1,]
  }
  return(corrvec)
}
