##' Plot functions for 'SRanalysis' objects. Allows to visualize the
##' stability and predictiveness trade-off of individual predictors.
##'
##' @title plot function
##' @param x object of class 'SRanalysis'.
##' @param x_axis either "SRdiff" or "SRpred".
##' @param varnames vector of variables names given in same ordering
##'   as columns of X. If NA the variable names saved in the
##'   SRanalysis object are used.
##' @param labels boolean specifying whether to print names for all
##'   variables with selection probability greater than 0.5. Only
##'   works if varnames has been specified.
##' @param ... arguments to be passed to or from other methods.
##'
##' @import ggplot2 ggrepel
##'
##' @author Niklas Pfister
##'
##' @export

plot.SRanalysis <- function(x, x_axis="SRdiff", varnames=NA, labels=FALSE, ...){
  stopifnot(inherits(x, "SRanalysis"))

  ## Read out results
  sp_stab <- x$results$SR$selection_probs
  siglevel_stab <- x$results$SR$siglevel
  if(x_axis == "SRdiff"){
    sp_pred <- x$results$SRdiff$selection_probs
    siglevel_pred <- x$results$SRdiff$siglevel
    xlabel <- "selection probability (SRdiff)"
  }
  else if(x_axis == "SRpred"){
    sp_pred <- x$results$SRpred$selection_probs
    siglevel_pred <- x$results$SRpred$siglevel
    xlabel <- "selection probability (SRpred)"
  }
  else{
    stop("x_axis needs to be either SRpred or SRdiff")
  }
  beta_sign <- x$avgcoefsign_SR[-1]
  if(length(varnames) != length(sp_stab)){
    if(is.na(varnames[1])){
      varnames <- x$varnames
    }
    else{
      warning("Length of supplied varnames does not correspond to results. Names from SRanalysis object are used.")
      varnames <- x$varnames
    }
  }

  df <- data.frame(sp_pred=sp_pred, sp_stab=sp_stab,
                   names=varnames,
                   beta_sign=beta_sign)
  splot <- ggplot()
  splot <- splot + geom_rect(aes(xmin=-Inf, xmax=0.5, ymin=-Inf, ymax=0.5), 
                             fill=rgb(1,0,0), alpha=0.2)
  splot <- splot + geom_rect(aes(xmin=siglevel_pred, xmax=Inf, ymin=-Inf, ymax=Inf), 
                             fill=rgb(0,1,0), alpha=0.2)
  splot <- splot + geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=siglevel_stab, ymax=Inf), 
                             fill=rgb(0,1,0), alpha=0.2)
  splot <- splot + geom_point(data=df, aes(x=sp_pred, y=sp_stab, color=beta_sign),
                              cex=1.5)
  splot <- splot + scale_color_gradient(low = "red", high = "black", limits=c(0,1))
  if(labels){
    splot <- splot + geom_text_repel(data=df, aes(x=sp_pred, y=sp_stab,
                                                  label=ifelse((sp_pred > 0.5 |
                                                                  sp_stab > 0.5),
                                                               as.character(names),'')),
                                     segment.size=0.2,
                                     force=0.1, hjust=0,
                                     segment.color="grey50",
                                     nudge_x=0, nudge_y=0, size=3.5)
  }
  splot <- splot + xlim(c(0,1)) + ylim(c(0,1))
  splot <- splot + xlab(xlabel)
  splot <- splot + ylab("selection probability (SR)")
  splot <- splot + theme(legend.position="none")
  return(splot)
}
