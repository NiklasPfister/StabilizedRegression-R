# StabilizedRegression-R

Implementation of StabilizedRegression, a framework for
multi-environment regression. Analyzes the functional dependence of
the response on a set of predictors and assess whether this functional
relationship remains stable or unstable across different environments.

N. Pfister, E. G. Williams, J. Peters, R. Aebersold, P. Bühlmann: Stabilizing Variable Selection and Regression. ArXiv preprint: 1911.01850. [https://arxiv.org/abs/1911.01850](https://arxiv.org/abs/1911.01850)

## Getting started

Clone this repository. Then, the package can be installed and loaded with the following commands:

```R
install.packages("PATHTOPACKAGE/StabilizedRegression_1.0.tar.gz")
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

To try out the package on the mouse data, run the example in mouse_example/example.R.
