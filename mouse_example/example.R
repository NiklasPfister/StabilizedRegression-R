## Load functions and packages
library(StabilizedRegression)


## Load data
load("data/pathway_info.Rda")
data <- read.csv("data/preprocessed_protein_data.csv")
pathway_list <- pathway_info$pathway_list_prot
Xmat <- as.matrix(data[,-(1:2)])
env <- data[,2]
varnames <- colnames(Xmat)

## Set variables
X <- Xmat[,varnames != "Hmgcs1"]
Y <- Xmat[,varnames == "Hmgcs1"]
A <- env

## Fix parameters
cores <- 4
pars <- list(m=6,
             B=100,
             alpha_stab=0.1,
             alpha_pred=0.01,
             size_weight="linear",
             prescreen_size=20,
             use_resampling=FALSE,
             stab_test="exact",
             variable_importance="scaled_coefficient")

fit_obj <- SRanalysis(X, Y, A, 12,
                      pars_SR=pars, cores=cores)

## Plot with SRdiff on the x_axis
plot(fit_obj, x_axis="SRdiff", varnames=colnames(X), labels=TRUE)

## Plot with SRpred on the x_axis
plot(fit_obj, x_axis="SRpred", varnames=colnames(X), labels=TRUE)
