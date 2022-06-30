# StabilizedRegression-R

Implementation of StabilizedRegression, a framework for
multi-environment regression. Analyzes the functional dependence of
the response on a set of predictors and assess whether this functional
relationship remains stable or unstable across different environments.

N. Pfister, E. Williams, R. Aebersold, P. Bühlmann: *Stabilizing Variable Selection and Regression*. Annals of Applied Statistics, 15(3):1220-1246. [https://doi.org/10.1214/21-AOAS1487](https://doi.org/10.1214/21-AOAS1487)


## Getting started

The package is available on CRAN. To install and load the package use the following command:
```R
install.packages(StabilizedRegression)
library(StabilizedRegression)
```

The main functions are StabilizedRegression(), which performs a multi-environment regression.

```R
help(StabilizedRegression)
example(StabilizedRegression)
```

The diagnostic plots can be created by using the function SRanalysis(). This will generate an object of class 'SRanalysis' which can then be plotted with using the plot function.

```R
help(SRanalysis)
example(SRanalysis)
```
