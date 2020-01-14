##' Coefficients functions for 'StabilizedRegression' objects.
##'
##' @title coefficients function
##' @param object object of class 'StabilizedRegression'.
##' @param predictive_model boolean specifying whether to use the
##' @param ... additional arguments affecting the summary produced.
##'
##' @author Niklas Pfister
##'
##' @export


coef.StabilizedRegression <- function(object, predictive_model=FALSE, ...){
  stopifnot(inherits(object, "StabilizedRegression"))

  if(!predictive_model){
    weighting <- object$weighting
  }
  else{
    if(length(object$weighting_pred) == 0){
      stop("Predictive model has not been computed.")
    }
    weighting <- object$weighting_pred
  }
  
  ## Aggregate coefficients by weighting
  d <- length(object$variable_importance)
  beta <- rep(0, d + 1)
  non_zeros <- which(weighting>0) 
  for(i in 1:length(non_zeros)){
    w <- weighting[non_zeros[i]]
    estimator <- object$learner_list[[non_zeros[i]]]$estimator
    S <- object$learner_list[[non_zeros[i]]]$S
    beta[c(1, S+1)] <- beta[c(1, S+1)] + estimator*w
  }
  
  return(beta) 
}
