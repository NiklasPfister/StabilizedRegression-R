##' StabilizedRegression based on linear OLS
##'
##' Performs a linear regression of a response \code{Y} on a set of
##' predictors \code{X} while ensuring stability across different
##' values of a stabilizing variable \code{A}.
##' @title StabilizedRegression
##' @param X predictor matrix. Numeric matrix of size n times d, where
##'   columns correspond to individual predictors.
##' @param Y response variable. Numeric vector of length n.
##' @param A stabilizing variable. Numeric vector of length n which
##'   can be interpreted as a factor.
##' @param pars list of additional parameters. \code{m} (default
##'   ncol(X)) integer specifying the largest possible subset
##'   size. \code{B} (default 100) integer specifying the number of
##'   random subsets to sample, if NA all subsets will be
##'   used. \code{alpha_stab} (default 0.05) value between 0 and 1
##'   specifiying the stability cutoff. \code{alpha_pred} (default
##'   0.05) value between 0 and 1 specifiying the predictive
##'   cutoff. \code{size_weight} (default "linear") one of the strings
##'   "linear", "constant", "quadratic", "rbf" or numeric weight
##'   vector specifying a probablity for each potential set size from
##'   1 to \code{m}. \code{compute_predictive_model} (default TRUE)
##'   boolean specifying whether to additionally compute SR (pred) and
##'   SR (diff) as well. \code{prescreen_size} (default NA) integer
##'   specifying the number of variables to screen down to before
##'   applying SR, if NA then no screening is
##'   applied. \code{prescreen_type} (default "correlation") one of
##'   the strings "correlation", "ols", "lasso", "deconfounding",
##'   "correlation_env", "deconfounding_env" specifying the type of
##'   screening. \code{stab_test} (default "exact") specifies which
##'   stability test to use. Either "exact" for a Bonferroni-corrected
##'   version of Chow's test, "mean_sres" a mean test based on
##'   resampling of the scaled residuals or "meanvar_sres" a mean and
##'   variance test based on resampling of the scaled
##'   residuals. \code{pred_score} (default "mse") specifies the
##'   prediction score. Either "mse" for the mean squared error,
##'   "mse_env" for the environment-wise best mean squared error,
##'   "aic" for the Akaike information criterion or "bic" for the
##'   Bayesian information criterion. \code{topk} (default 1) is a
##'   tuning parameter that can be used to increase the number of
##'   predictive sets. It should be an integer value, where higher
##'   values lead to more accepted sets based on the predictive
##'   cutoff. \code{variable_importance} (default
##'   "scaled_coefficient") specifies the type of variable
##'   ranking. Either "weighted" for a weighted average of all
##'   selected subsets, "scaled_coefficient" for a ranking based on
##'   the scaled average regression parameter or "permutation" for a
##'   permutation based ranking.
##' @param verbose 0 for no output, 1 for text output and 2 for text
##'   and diagnostic plots.
##' @param seed fix the seed value at the beginning of the function.
##' 
##' @return Object of class 'StabilizedRegression' consisting of the following
##'   elements
##'
##' \item{learner_list}{List of all fitted linear OLS regressions (fitted R6 'linear_regression' objects).}
##' \item{weighting}{Weighting of the individual regressions in SR.}
##' \item{weighting_pred}{Weighting of the individual regressions in SR (pred). Only computed if compute_predictive_model is TRUE.}
##' \item{variable_importance}{Variable importance measure for all predictors based on SR.}
##' \item{variable_importance_pred}{Variable importance measure for all predictors based on SR (pred). Only computed if compute_predictive_model is TRUE.}
##' \item{variable_importance_diff}{Variable importance measure for all predictors based on difference between SR and SR (pred). Only computed if compute_predictive_model is TRUE.}
##' 
##' @export
##'
##' @import stats utils MASS R6 graphics grDevices
##'
##' @references Pfister, N., E. Williams, R. Aebersold, J. Peters and
##'   P. B{\"u}hlmann (2019). Stabilizing Variable Selection and
##'   Regression. arXiv preprint arXiv:1911.01850.
##'
##' @author Niklas Pfister
##'
##' @examples
##' ## Example
##' set.seed(1)
##' X1 <- rnorm(200)
##' Y <- X1 + rnorm(200)
##' X2 <- 0.5 * X1 + Y + 0.2 * c(rnorm(100), rnorm(100)+2)
##'
##' X <- cbind(X1, X2)
##' A <- as.factor(rep(c(0, 1), each=100))
##'
##' fit_sr <- StabilizedRegression(X, Y, A, pars=list(B=NA))
##' fit_lm <- lm(Y ~ X)
##'
##' print(paste("Coefficients of SR:", toString(coefficients(fit_sr))))
##' print(paste("Coefficients of SR (pred):", toString(coefficients(fit_sr, predictive_model=TRUE))))
##' print(paste("Coefficients of OLS:", toString(coefficients(fit_lm))))


StabilizedRegression <- function(X, Y, A,
                                 pars=list(m=ncol(X),
                                           B=100,
                                           alpha_stab=0.05,
                                           alpha_pred=0.05,
                                           size_weight="linear",
                                           compute_predictive_model=TRUE,
                                           use_resampling=FALSE,
                                           prescreen_size=NA,
                                           prescreen_type="correlation",
                                           stab_test="exact",
                                           pred_score="mse",
                                           topk=1,
                                           variable_importance="scaled_coefficient"),
                                 verbose=0,
                                 seed=NA){

  ## Set defaults for missing parameters
  if(!exists("m", pars)){
    pars$m <- ncol(X)
  }
  if(!exists("B", pars)){
    pars$B <- 100
  }
  if(!exists("alpha_stab", pars)){
    pars$alpha_stab <- 0.05
  }
  if(!exists("alpha_pred", pars)){
    pars$alpha_pred <- 0.05
  }
  if(!exists("size_weight", pars)){
    pars$size_weight <- "linear"
  }
  if(!exists("compute_predictive_model", pars)){
    pars$compute_predictive_model <- TRUE
  }
  if(!exists("use_resampling", pars)){
    pars$use_resampling <- FALSE
  }
  if(!exists("prescreen_size", pars)){
    pars$prescreen_size <- NA
  }
  if(!exists("prescreen_type", pars)){
    pars$prescreen_type <- "correlation"
  }
  if(!exists("stab_test", pars)){
    pars$stab_test <- "exact"
  }
  if(!exists("pred_score", pars)){
    pars$pred_score <- "mse"
  }
  if (!exists("topk", pars)) {
    pars$topk <- 1
  }
  if(!exists("variable_importance", pars)){
    pars$variable_importance <- "scaled_coefficient"
  }

  ## Set seed
  if(!is.na(seed)){
    set.seed(seed)
  }

  ## Read out parameters and perform checks
  n <- nrow(X)
  d <- ncol(X)
  if(is.na(A)[1] | length(unique(A)) == 1){
    warning("A either not specified correctly, contains only one environment or is set to NA.\n No stability test will be performed and all sets will be considered stable.")
    A <- as.factor(rep(1, n))
    pars$stab_test  <-  "none"
  }
  else if(!is.factor(A)){
    warning("A has been converted to a factor variable.")
    A <- as.factor(A)
  }
  if(length(pars$pred_score) == 1){
    pars$pred_score <- rep(pars$pred_score, 2)
  }
  else if(length(pars$pred_score) != 2){
    stop("Incorrect length of pred_score.")
  }
  regression_pars <- list(test=pars$stab_test,
                          pred_score=pars$pred_score)
  
  ## Remove row and column names as these can cause matching issues
  X <- unname(X)
  A <- unname(A)
  Y <- unname(Y)
  

  ## Prescreen
  if(is.numeric(pars$prescreen_size)){
    pars$prescreen_size <- min(c(pars$prescreen_size, ncol(X)))
    if(pars$prescreen_type == "lasso"){
      fit <- glmnet(X, Y)
      sel.matrix <- (fit$beta != 0)
      first.entrance <- apply(sel.matrix, 1, which.max)
      first.entrance[which(apply(sel.matrix, 1, sum) == 0)] <- Inf
      tmpvec <- order(first.entrance, decreasing=FALSE)
      num_noinf <- sum(first.entrance != Inf)
      if(length(tmpvec)-num_noinf>1){
        tmpvec[(num_noinf+1):length(tmpvec)] <- sample(tmpvec[(num_noinf+1):length(tmpvec)])
      }
      screened_vars <- tmpvec[1:pars$prescreen_size]
    }
    else if(pars$prescreen_type == "ols"){
      if(n <= d){
        stop(" n <= d: ols screening can only be used in low-dimensional settings")
      }
      fit <- lm(Y ~ X)
      pval <- summary(fit)$coefficients[-1,'Pr(>|t|)']
      screened_vars <- which(pval <= 0.1)
    }
    else if(pars$prescreen_type == "correlation"){
      pval <- apply(X, 2, function(x) cor.test(Y, x)$p.value)
      screened_vars <- order(pval)[1:pars$prescreen_size]
    }
    else if(pars$prescreen_type == "correlation_env"){
      Alist <- lapply(unique(A), function(a) which(A == a))
      pval <- apply(X, 2, function(x) sapply(Alist,
                                             function(Aind)
                                               cor.test(Y[Aind], x[Aind])$p.value))
      pval <- apply(pval, 2, min)
      screened_vars <- order(pval)[1:pars$prescreen_size]
    }
    else if(pars$prescreen_type == "deconfounding"){
      corr_deconf <- abs(deconfounding_correlation(X, Y)[-1])
      screened_vars <- order(corr_deconf, decreasing=TRUE)[1:pars$prescreen_size]
    }
    else if(pars$prescreen_type == "deconfounding_env"){
      Alist <- lapply(unique(A), function(a) which(A == a))
      corr_deconf <- sapply(Alist, function(Aind)
        abs(deconfounding_correlation(X[Aind,,drop=FALSE], Y[Aind])[-1]))
      corr_deconf <- apply(corr_deconf, 1, max)
      screened_vars <- order(corr_deconf, decreasing=TRUE)[1:pars$prescreen_size]
    }
    else{
      stop("Selected prescreen_type does not exist.")
    }
  }
  else{
    screened_vars <- 1:d
  }
  m <- min(c(pars$m, length(screened_vars), nrow(X)-1))



  ## Define function for a single iteration
  Xscreened <- X[,screened_vars]
  single_iteration <- function(S, extra){
    # initialize R6 regressor (linear_regressor)
    regobj <- linear_regressor$new(S=S, pars=regression_pars)
    # fit model
    if(pars$use_resampling){
      ind <- sample(1:n, floor(n/2), replace=TRUE)
    }
    else{
      ind <- 1:n
    }
    regobj$fit(Xscreened[ind,, drop=FALSE], Y[ind], A[ind], extra)
    return(regobj)
  }
                                   
  ## Compute list of weak learners
  if(is.numeric(pars$B)){
    if(!pars$use_resampling & pars$B >= 2^(length(screened_vars))){
      warning("Using all sets, since pars$B is larger than maximum number of permitted subsets and pars$use_resampling is FALSE.")
      pars$B <- NA
    }
  }
  if(!is.numeric(pars$B)){
    sets <- list()
    for(i in 1:m){
      sets <- c(sets, combn(1:length(screened_vars), i, simplify=FALSE))
    }
    pars$B <- length(sets)
  }
  else{
    # Determine how to weight set sizes
    if(m == 1){
      pars$size_weight <- 1
    }
    else if(pars$size_weight == "const"){
      pars$size_weight <- rep(1/m, m)
    }
    else if(pars$size_weight == "linear"){
      pars$set_size_weight <- c(1:floor(m/2), ceiling(m/2):1)
      pars$set_size_weight <- pars$set_size_weight/sum(pars$set_size_weight)
    }
    else if(pars$size_weight == "quadratic"){
      pars$set_size_weight <- c(1:floor(m/2), ceiling(m/2):1)^2
      pars$set_size_weight <- pars$set_size_weight/sum(pars$set_size_weight)
    }
    else if(pars$size_weight == "rbf"){
      pars$set_size_weight <- exp(-c(-(ceiling(m/2):1), 1:floor(m/2))^2/m)
      pars$set_size_weight <- pars$set_size_weight/sum(pars$set_size_weight)
    }
    else if(is.numeric(pars$size_weight)){
      pars$set_size_weight <- rep(0, m)
      pars$set_size_weight[pars$size_weight] <- 1
      pars$set_size_weight <- pars$set_size_weight/sum(pars$set_size_weight)
    }
    else{
      stop("Specified size_weight does not exist.")
    }
    # Sample subsets
    sets <- lapply(1:pars$B, function(i){
      set_size <- sample(1:m, 1, prob=pars$set_size_weight)
      return(sort(sample(1:length(screened_vars), set_size)))
    })
  }
  sets <- c(list(numeric()), sets)

  ## Fit each weak learner
  ptm <- proc.time()
  learner_list <- lapply(sets, single_iteration, extra=NA)
  if(verbose > 0){
    print(proc.time()-ptm)
  }

  ## Adjust sets for screening
  learner_list <- lapply(learner_list,
                         function(x){
                           x$S <- screened_vars[x$S]
                           return(x)
                         })

  ## Compute scores for each weak learner and determine weighting
  ptm <- proc.time()
  numL <- length(learner_list)
  numB <- 500
  
  scores <- sapply(learner_list, function(x) x$scores)
  pval_stab <- scores[1,]
  mse <- scores[2,]
  mse_pred <- scores[3,]

  # Compute stable models
  if(!is.na(pars$alpha_stab)){
    stabmods <- pval_stab >= pars$alpha_stab
  }
  else{
    warning("alpha_stab is NA, all sets were selected")
    stabmods <- rep(TRUE, length(pval_stab))
  }
  if(sum(stabmods) == 0){
    warning(paste("No significantly stable sets at level alpha_stab.",
                  "Top 10 percent of sets were selected"))
    stabmods <- pval_stab >= pval_stab[
      order(pval_stab, decreasing=TRUE)[ceiling(0.1*length(pval_stab))]]
  }

  ## Compute predictive models
  if(is.na(pars$alpha_pred)){
    predmods_all <- rep(TRUE, length(mse))
    predmods_stab <- stabmods
    weighting <- rep(0, length(learner_list))
    weighting[stabmods & predmods_stab] <- 1
    weighting <- weighting/sum(weighting)
    weighting_pred <- rep(0, length(learner_list))
    weighting_pred[predmods_all] <- 1
    weighting_pred <- weighting_pred/sum(weighting_pred)
  }
  else if(pars$pred_score[1] %in% c("mse", "mse_env", "expvar_env")){
    if(!(pars$pred_score[2] %in% c("mse", "mse_env", "expvar_env"))){
      warning("Combination of prediction score does not seem to be correct!")
    }
    if(pars$use_resampling){
      topk <- 0.1*sum(stabmods)
    }
    else{
      topk <- min(c(pars$topk, sum(stabmods)))
    }
    bestmod_all <- order(mse_pred, decreasing=FALSE)[1:topk]
    bestmod_stab <- order(mse[stabmods], decreasing=FALSE)[1:topk]
    cutoff_stab <- -Inf
    cutoff_all <- -Inf
    for(k in 1:topk){
      # All models
      Xtmp <- X[,learner_list[[bestmod_all[k]]]$S,drop=FALSE]
      bootstrap_all <- bootstrap_mse(Y, cbind(rep(1, nrow(Xtmp)), Xtmp), A, numB,
                                     pars$pred_score[2])
      cutoff_all <- max(c(quantile(bootstrap_all, 1-pars$alpha_pred), cutoff_all))
      # Stable models
      Xtmp <- X[,learner_list[stabmods][[bestmod_stab[k]]]$S,drop=FALSE]
      bootstrap_stab <- bootstrap_mse(Y, cbind(rep(1, nrow(Xtmp)), Xtmp), A, numB,
                                      pars$pred_score[1])
      cutoff_stab <- max(c(quantile(bootstrap_stab, 1-pars$alpha_pred), cutoff_stab))
    }
    predmods_all <- mse_pred <= cutoff_all
    predmods_stab <- mse <= cutoff_stab
    
    # Ensure that at least one predictive model was selected
    if(sum(predmods_all) == 0){
      predmods_all <- mse_pred == min(mse_pred)
    }
    if(sum(predmods_stab & stabmods) == 0){
      predmods_stab <- rep(FALSE, length(mse))
      predmods_stab[stabmods] <- mse[stabmods] == min(mse[stabmods])
    }
    # Set weighting
    weighting <- rep(0, length(learner_list))
    weighting[stabmods & predmods_stab] <- 1
    weighting <- weighting/sum(weighting)
    weighting_pred <- rep(0, length(learner_list))
    weighting_pred[predmods_all] <- 1
    weighting_pred <- weighting_pred/sum(weighting_pred)
  }
  else{
    weighting <- rep(0, length(learner_list))
    deltaAIC <- mse[stabmods]-min(mse[stabmods])
    weighting[stabmods] <- exp(-0.5*deltaAIC)/sum(exp(-0.5*deltaAIC))
    deltaAIC <- mse_pred-min(mse_pred)
    weighting_pred <- exp(-0.5*deltaAIC)/sum(exp(-0.5*deltaAIC))
  }

  
  if(verbose > 0){
    print(proc.time()-ptm)
  }

  ## Plot results
  if(verbose > 1){
    par(mfrow=c(1, 3))
    plot(mse, pval_stab)
    plot(mse, pval_stab,
         pch=19, col=rgb(0, 0, 0, weighting/max(weighting)))
    plot(mse_pred, pval_stab,
         pch=19, col=rgb(0, 0, 0, weighting_pred/max(weighting_pred)))
  }
 
  ## Collect output
  res <- list(learner_list=learner_list,
              weighting=weighting)
  if(pars$compute_predictive_model){
    res$weighting_pred <- weighting_pred
  }
  class(res) <- "StabilizedRegression"


  ## Compute variable importance
  if(pars$variable_importance=="weighted"){
    res$variable_importance <- rep(0, d)
    for(i in 1:length(learner_list)){
      set <- learner_list[[i]]$S
      res$variable_importance[set] <- res$variable_importance[set] + weighting[i]
    }
    if(pars$compute_predictive_model){
      res$variable_importance_pred <- rep(0, d)
      for(i in 1:length(learner_list)){
        set <- learner_list[[i]]$S
        res$variable_importance_pred[set] <- res$variable_importance_pred[set] + weighting_pred[i]
      }
    }
  }
  else if(pars$variable_importance=="permutation"){
    Btmp <- 100
    res$variable_importance <- VariableImportance(X, Y, res,
                                                  B=Btmp)
    if(pars$compute_predictive_model){
      res$weighting <- weighting_pred
      res$variable_importance_pred <- VariableImportance(X, Y, res,
                                                         B=Btmp)
      res$weighting <- weighting
    }
  }
  else if(pars$variable_importance=="scaled_coefficient"){
    coef <- rep(0, ncol(X))
    for(i in 1:length(learner_list)){
      S <- learner_list[[i]]$S
      coef[S] <- coef[S] + weighting[i]*learner_list[[i]]$estimator[-1]
    }
    res$variable_importance <- abs(coef*apply(X, 2, sd))
    if(pars$compute_predictive_model){
      coef <- rep(0, ncol(X))
      for(i in 1:length(learner_list)){
        S <- learner_list[[i]]$S
        coef[S] <- coef[S] + weighting_pred[i]*learner_list[[i]]$estimator[-1]
      }
      res$variable_importance_pred <- abs(coef*apply(X, 2, sd))
    }
  }
  else{
    warning("Specified variable_importance does not exist! No variable importance was computed.")
    res$variable_importance <- NA
    if(pars$compute_predictive_model){
      res$variable_importance_pred <- NA
    }
  }


  ## Compute variable importance for predictive but non-stable variables
  if(pars$compute_predictive_model){
    res$variable_importance_diff <- res$variable_importance_pred - res$variable_importance
  }
  else{
    res$variable_importance_diff <- NA
  }

  
  return(res)
}
