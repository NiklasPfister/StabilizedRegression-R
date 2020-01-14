##' Predict functions for 'StabilizedRegression' objects.
##'
##' @title predict function
##' @param object object of class 'StabilizedRegression'.
##' @param newdata matrix or data.frame for which the response should
##'   be predicted.
##' @param predictive_model boolean specifying whether to use the
##' @param ... additional arguments affecting the prediction produced.
##'
##' @author Niklas Pfister
##'
##' @export


predict.StabilizedRegression <- function(object, newdata, predictive_model=FALSE, ...){
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
  
  ## Aggregate predictions either using averages or majority vote
  Ypred <- rep(0, nrow(newdata))
  non_zeros <- which(weighting>0)
  for(i in 1:length(non_zeros)){
    w <- weighting[non_zeros[i]]
    regobj <- object$regressor$new()
    regobj$estimator <- object$learner_list[[non_zeros[i]]]$estimator
    regobj$S <- object$learner_list[[non_zeros[i]]]$S
    Ypred <- Ypred + w * regobj$predict(newdata)
  }
  
  return(Ypred) 
}
