##' Stability analysis based on stabilized regression used to analyze
##' the trade-off between stability and predictivness of individual
##' predictors.
##'
##' This function performs two version of StabilizedRegression: SR
##' which selects a stable and predictive model and SRpred which fits
##' a plain predictive model. Stability selection is then performed
##' using the variable importance measures from both these methods and
##' from their difference SRdiff as variable selection criterion. This
##' allows to distinguish between which predictive variables are
##' stable and which are unstable with respect to the stabilizing
##' variable A. The results can be visualized by plotting the
##' resulting object using the plot() function.
##'
##' Due to the resampling this function can be quite computationally
##' involved, we therefore recommend making use of the \code{cores}
##' parameter for parallel computations.
##' 
##' @title Stability analysis
##' @param X predictor matrix. Numeric matrix of size n times d, where
##'   columns correspond to individual predictors.
##' @param Y response variable. Numeric vector of length n.
##' @param A stabilizing variable. Numeric vector of length n which
##'   can be interpreted as a factor.
##' @param num_reps number of resamples to use in stability selection.
##' @param pred_scores characeter vector of length 2, specifying the
##'   \code{pred_score} for SR and SRpred.
##' @param prescreen_types characeter vector of length 2, specifying
##'   the \code{prescreen_type} for SR and SRpred.
##' @param pars_SR list of all remaining parameters going into
##'   StabilizedRegression. \code{compute_predictive},
##'   \code{pred_score} and \code{prescreen_type} are ignored.
##' @param threshold numeric value between 0 and 1, specifying in
##'   stability selection at which value to select variables.
##' @param cores number of cores used in mclapply.
##' @param verbose 0 for no output, 1 for text output and 2 for text
##'   and diagnostic plots.
##' @param seed fix the seed value at the beginning of the function.
##' 
##' @return Object of class 'SRanalysis' consisting of the following
##'   elements
##'
##' \item{results}{List of stability selection results for for SR, SRpred and SRdiff.}
##' \item{varnames}{Vector of variable names taken from the column names of X.}
##' \item{avgcoefsign_SR}{Vector of average coefficient signs for SR}
##' \item{avgcoefsign_SRpred}{Vector of average coefficient signs for SRpred}
##' 
##' @export
##'
##' @import stats utils MASS R6
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
##' X2 <- 0.5 * X1 + Y + 0.2 * c(rnorm(100), rnorm(100)+3)
##'
##' X <- cbind(X1, X2)
##' A <- as.factor(rep(c(0, 1), each=100))
##'
##' obj <- SRanalysis(X, Y, A, 10,
##'                   pars_SR=list(B=NA))
##' plot(obj, varnames = c("X1", "X2"), labels=TRUE)
##' print(obj$results)


SRanalysis <- function(X, Y, A,
                       num_reps=100,
                       pred_scores=c("mse", "mse_env"),
                       prescreen_types=c("correlation", "correlation_env"),
                       pars_SR=list(m=ncol(X),
                                    B=100,
                                    alpha_stab=0.05,
                                    alpha_pred=0.05,
                                    size_weight="linear",
                                    use_resampling=FALSE,
                                    prescreen_size=NA,
                                    stab_test="exact",
                                    variable_importance="scaled_coefficient"),
                       threshold=0,
                       cores=1,
                       verbose=0,
                       seed=NA){


  ## Set seed
  if(!is.na(seed)){
    set.seed(seed)
  }

  ## Set parameters
  n <- nrow(X)

  ## Define function to compute single iteration of variable selection
  if(length(unique(pred_scores)) == 1 & length(unique(prescreen_types)) == 1){
    single_iteration <- function(Xsub, Ysub, Asub){
      # SR
      pars_SR$compute_predictive_model <- TRUE
      pars_SR$pred_score <- pred_scores[1]
      pars_SR$prescreen_type <- prescreen_types[1]
      fit <- StabilizedRegression(Xsub, Ysub, Asub,
                                  pars_SR,
                                  verbose=0)
      vs_stab <- fit$variable_importance
      beta_stab <- coef(fit)
      # SRpred
      vs_pred <- fit$variable_importance_pred
      beta_pred <- coef(fit, predictive_model=TRUE)
      # SRdiff
      vs_diff <- vs_pred - vs_stab
      
      return(list(SR=vs_stab,
                  SRpred=vs_pred,
                  SRdiff=vs_diff,
                  beta_stab=beta_stab,
                  beta_pred=beta_pred))
    }
  }
  else{
    single_iteration <- function(Xsub, Ysub, Asub){
      # SR
      pars_SR$compute_predictive_model <- FALSE
      pars_SR$pred_score <- pred_scores[1]
      pars_SR$prescreen_type <- prescreen_types[1]
      fit <- StabilizedRegression(Xsub, Ysub, Asub,
                                  pars_SR,
                                  verbose=0)
      vs_stab <- fit$variable_importance
      beta_stab <- coef(fit)
      # SRpred
      pars_SR$compute_predictive_model <- TRUE
      pars_SR$pred_score <- pred_scores[2]
      pars_SR$prescreen_type <- prescreen_types[2]
      fit <- StabilizedRegression(Xsub, Ysub, Asub,
                                  pars_SR,
                                  verbose=0)
      vs_pred <- fit$variable_importance_pred
      beta_pred <- coef(fit, predictive_model=TRUE)
      # SRdiff
      vs_diff <- vs_pred - vs_stab
      
      return(list(SR=vs_stab,
                  SRpred=vs_pred,
                  SRdiff=vs_diff,
                  beta_stab=beta_stab,
                  beta_pred=beta_pred))
    }
  }

  ## Apply stability selection resampling
  indlist <- lapply(1:num_reps, function(i) sample(1:n, floor(n/2), replace=TRUE))
  resample_res <- mclapply(indlist,
                           function(ind) single_iteration(X[ind,,drop=FALSE], Y[ind], A[ind]),
                           mc.cores=cores, mc.preschedule=FALSE)
  ## resample_res <- lapply(indlist,
  ##                        function(ind) single_iteration(X[ind,,drop=FALSE], Y[ind], A[ind]))
  
  ## Compute selection probabilities
  vs_ind <- startsWith(names(resample_res[[1]]), "SR")
  num_methods <- sum(vs_ind)
  method_name <- vector("numeric", num_methods)
  results <- vector("list", num_methods)
  for(i in 1:num_methods){
    method_name[i] <- as.character(names(resample_res[[1]][vs_ind][i]))
    vs_score <- t(sapply(resample_res, function(x) x[vs_ind][[i]]))
    ## Compute selection probabilities
    selection_ind <- matrix(FALSE, nrow(vs_score), ncol(vs_score))
    q_est <- mean(apply(cbind(rowSums(vs_score > threshold),
                              rep(floor(sqrt(ncol(vs_score)*0.5*1)),
                                  nrow(vs_score))),
                        1, min))
    for(k in 1:nrow(vs_score)){
      tmpvec <- order(vs_score[k,], decreasing=TRUE)
      num_select <- min(c(sum(vs_score[k,] > threshold),
                          floor(sqrt(ncol(vs_score)*0.5*1))), na.rm=TRUE)
      selection_ind[k, tmpvec[1:num_select]] <- TRUE
    }
    selection_probs <- colMeans(selection_ind)
    false_discovery <- ceiling(1/(2*selection_probs-1)*(q_est^2)/(ncol(vs_score)))
    false_discovery[selection_probs <= 0.5] <- NA  
    siglevel <- (q_est^2/ncol(vs_score)+1)/2
    ## Collect results
    results[[i]]$selection_probs <- selection_probs
    results[[i]]$false_discovery <- false_discovery
    results[[i]]$siglevel <- siglevel
    results[[i]]$q_est <- q_est
  }
  names(results) <- method_name

  ## Compute average coefficient sign
  betamat_pred <- sapply(resample_res, function(x) sign(x$beta_pred))
  avgcoefsign_SRpred <- apply(betamat_pred, 1,
                              function(x) ifelse(sum(x != 0) == 0, 0.5, mean(x[x != 0]==1)))
  betamat_stab <- sapply(resample_res, function(x) sign(x$beta_stab))
  avgcoefsign_SR <- apply(betamat_stab, 1,
                          function(x) ifelse(sum(x != 0) == 0, 0.5, mean(x[x != 0]==1)))

  ## Collect results
  res <- list(results=results,
              varnames=colnames(X),
              avgcoefsign_SR=avgcoefsign_SR,
              avgcoefsign_SRpred=avgcoefsign_SRpred)
  class(res) <- "SRanalysis"
  
  
  return(res)
}
