##' R6 Class Representing a Linear Regression
##'
##' @description An R6-class for linear regression that is used
##'   within the StabilizedRegression framework.
##'
##'   Currently this is the only regression procedure that has been
##'   implemented. In order to extend the StabilizedRegression
##'   framework to a different regression procedure a custom R6-class
##'   with the same structure as this function can be written and used
##'   within StabilizedRegression.
##'
##' @details Constructer method initializes a linear regression object
##'   specifying on which subset of variables \code{S} to fit the
##'   regression and which type of stability test and prediction score
##'   to compute. The methods \code{fit()} and \code{predict()} can be
##'   applied to the object to fit and predict, respectively.
##'
##' @export
##' 
##' @import stats utils R6
##'
##' @author Niklas Pfister


linear_regressor <- R6Class("linear_regressor",
                            public = list(
                              #' @field estimator Numeric vector of
                              #' regression coefficients.
                              estimator = numeric(),
                              #' @field S Numeric vector specifying
                              #' the subset of variables to perform
                              #' regression on.
                              S = numeric(),
                              #' @field scores Numeric vector of
                              #' fitted stability and prediction scores.
                              scores = numeric(),
                              #' @field pars List specifying the
                              #' stability test via \code{test} and
                              #' prediction score via \code{pred_score}.
                              pars = list(),
                              #' @description
                              #' Create a new linear_regression object.
                              #' @param S Subset of variables.
                              #' @param pars Parameters.
                              #' @return A new `linear_regression` object.
                              initialize = function(S=numeric(),
                                                    pars=list(test="mean",
                                                              pred_score=c("mse",
                                                                           "mse"))){
                                self$S=S
                                if(length(pars$test) == 0){
                                  pars$test <- "mean"
                                  pars$pred_score=c("mse", "mse")
                                }
                                self$pars=pars
                              },
                              #' @description
                              #' Fit a 'linear_regression' object on data
                              #' and computes the stability and prediction scores.                              
                              #' @param X Predictor matrix.
                              #' @param Y response vector.
                              #' @param A environemnt indicator.
                              #' @param extra not required (placeholder)
                              #' @return A fitted `linear_regression` object.
                              fit = function(X, Y, A, extra=NA){
                                ## Pooled estimator
                                X <- cbind(rep(1, nrow(X)), X[, self$S, drop=FALSE])
                                Alist <- lapply(unique(A), function(a) which(A == a))
                                res <- getpval(Y, X, Alist, test=self$pars$test,
                                               pred_score=self$pars$pred_score,
                                               extra=extra)
                                self$estimator <- matrix(res$coefficients, ncol=1)
                                # make sure no sigularities occured
                                if(sum(is.na(self$estimator))>0){
                                  nonna <- !is.na(self$estimator)
                                  self$S <- self$S[nonna[-1]]
                                  self$estimator <- self$estimator[nonna,,drop=FALSE]
                                }
                                self$scores <- c(res$pval, res$pred_score)
                              },
                              #' @description
                              #' Predict using a fitted 'linear_regression' object.
                              #' @param X Predictor matrix on which to predict response.
                              #' @return Numeric vector of predicted response.
                              predict = function(X){
                                X <- cbind(rep(1, nrow(X)), X[, self$S, drop=FALSE])
                                invisible(X %*% self$estimator)
                              })
                            )

getpval <- function(Y, X, Alist, maxNoObs=1000,
                    test="exact", pred_score=c("mse", "mse"),
                    extra=NA){

  linm <- lm(Y~-1+X)
  coefficients <- coefficients(linm)
  pred <- fitted.values(linm)
  mse <- mean(residuals(linm)^2)
  predscore <- vector("numeric", 2)
  for(i in 1:2){
    if(pred_score[i] == "mse"){
      predscore[i] <- mse
    }
    else if(pred_score[i] == "mse_env"){
      predscore[i] <- min(sapply(Alist, function(Aind)
        mean(residuals(lm.fit(X[Aind,,drop=FALSE], Y[Aind]))^2)))
    }
    else if(pred_score[i] == "expvar_env"){
      predscore[i] <- min(sapply(Alist, function(Aind)
        mean(residuals(
          lm.fit(X[Aind,,drop=FALSE],
                 Y[Aind]))^2)/mean((Y[Aind]-mean(Y[Aind]))^2)))
    }
    else if(pred_score[i] == "aic"){
      predscore[i] <- extractAIC(linm, k=2)[2]
    }
    else if(pred_score[i] == "bic"){
      predscore[i] <- extractAIC(linm, k=log(length(Y)))[2]
    }
  }
  
  K <- length(Alist)
  n <- nrow(X)
  if(test == "exact"){
    pvalvec <- numeric(length(Alist))
    for (ki in 1:length(Alist)){
      nk <- length(Alist[[ki]])
      nko <- n-nk
      p <- ncol(X)
      linm <- lm.fit(X[-Alist[[ki]], ,drop=FALSE], Y[-Alist[[ki]]])
      is_ok <- !is.na(linm$coefficients)
      pred <- as.numeric(X[Alist[[ki]], is_ok ,drop=FALSE] %*% linm$coefficients[is_ok])
      diff <- Y[Alist[[ki]]] - pred
      
      selobs <-  if(nk>maxNoObs)  sample(Alist[[ki]], maxNoObs) else Alist[[ki]]
      if(nk>maxNoObs){
        diff <- diff[ Alist[[ki]] %in% selobs]
      }
      nk <- length(selobs)
      COV <- diag(length(diff)) + X[selobs, is_ok, drop=FALSE] %*% solve(
        t(X[-Alist[[ki]], is_ok, drop=FALSE]) %*% X[-Alist[[ki]], is_ok, drop=FALSE],
        t(X[selobs, is_ok, drop=FALSE]))
      stat <- (t(diff)%*% solve(COV, diff)) / (nk * var(residuals(linm))*nko/(nko-p))
      pval <- 1-pf(stat, nk, n-nk - ncol(X))
      
      pvalvec[ki] <- pval
    }
    if(length(Alist) == 2){
      pval <- min(pvalvec)
    }
    else{
      pval <- min(pvalvec) * length(Alist)
    }
    pval <- min(1, pval)
  }
  else if(test == "mean_sres"){
    n <- length(Y)
    R.scaled <- vector("list", length(Alist))
    Projections <- vector("list", length(Alist))
    for(i in 1:K){
      fit <- lm.fit(X[Alist[[i]],,drop=FALSE], Y[Alist[[i]]])
      is_ok <- !is.na(fit$coefficients)
      Rsmall <- qr.R(fit$qr)[is_ok,,drop=FALSE][,is_ok,drop=FALSE]
      Projections[[i]] <- (X[,is_ok,drop=FALSE] %*%
                             backsolve(Rsmall,
                                       t(qr.Q(fit$qr)[,is_ok,drop=FALSE])))
      R.scaled[[i]] <- Y - Projections[[i]] %*% Y[Alist[[i]]]
      R.scaled[[i]] <- R.scaled[[i]]/sd(R.scaled[[i]])
    }
    pairwise <- t(combn(K, 2))
    resample_fn <- function(res){
      meanvec <- sapply(Alist, function(Aind) mean(res[Aind]))
      statfun_mean <- function(i, j){
        T0 <- abs(meanvec[i] - meanvec[j])
        return(T0)
      }
      return(sum(mapply(statfun_mean, pairwise[, 1], pairwise[, 2])))
    }
    stab_score <- sum(sapply(R.scaled, function(res) resample_fn(res)))
    stab_scores_null <- sapply(1:100,
                               function(i)
                                 sum(
                                   sapply(
                                     1:length(Alist),
                                     function(k){
                                       res <- rnorm(n)
                                       res <- res - Projections[[k]] %*% res[Alist[[k]]]
                                       return(resample_fn(res))
                                     })))
    pval <- (sum(stab_scores_null>=stab_score)+1)/(length(stab_scores_null)+1)
  }
  else if(test == "meanvar_sres"){
    n <- length(Y)
    R.scaled <- vector("list", length(Alist))
    Projections <- vector("list", length(Alist))
    for(i in 1:K){
      fit <- lm.fit(X[Alist[[i]],,drop=FALSE], Y[Alist[[i]]])
      is_ok <- !is.na(fit$coefficients)
      Rsmall <- qr.R(fit$qr)[is_ok,,drop=FALSE][,is_ok,drop=FALSE]
      Projections[[i]] <- (X[,is_ok,drop=FALSE] %*%
                             backsolve(Rsmall,
                                       t(qr.Q(fit$qr)[,is_ok,drop=FALSE])))
      R.scaled[[i]] <- Y - Projections[[i]] %*% Y[Alist[[i]]]
      R.scaled[[i]] <- R.scaled[[i]]/sd(R.scaled[[i]])
    }
    pairwise <- t(combn(K, 2))
    Aweights <- sapply(Alist, length)/n
    resample_fn <- function(res){
      momvec <- sapply(Alist, function(Aind) c(mean(res[Aind]), var(res[Aind])))
      statfun <- function(i, j){
        T0 <- abs(momvec[,i] - momvec[,j])*(Aweights[i]*Aweights[j])
        return(T0)
      }
      return(rowSums(mapply(statfun, pairwise[, 1], pairwise[, 2])))
    }
    stab_scores <- colSums(t(sapply(R.scaled, function(res) resample_fn(res)))*Aweights)
    stab_scores_null <- sapply(1:100,
                               function(i)
                                 colSums(t(sapply(
                                   1:length(Alist),
                                   function(k){
                                     res <- rnorm(n)
                                     res <- res - Projections[[k]] %*% res[Alist[[k]]]
                                     return(resample_fn(res))
                                   }))*Aweights))
    pval <- min(c(1,
    (rowSums(stab_scores_null >= stab_scores)+1)/(ncol(stab_scores_null)+1)*2))
  }
  else if(test == "none"){
    pval <- 1
  }
  else{
    stop(paste("Stability test", test, "does not exist."))
  }

  return(list(pval=pval,
              coefficients=coefficients,
              pred_score=predscore))
}
