# culex

A discretized stochastic and deterministic model of temperature and photoperiod driven culex dynamics with movement of adults based on the mathematical model(s) presented in:

1.  Paper: Uncovering mechanisms behind mosquito seasonality by integrating mathematical models and daily empirical population data: Culex pipiens in the UK [here](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-019-3321-2), Code repo: [here](https://github.com/davewi13/Mosquito-seasonality-paper)
2.  Paper: Modelling the effect of temperature on the seasonal population dynamics of temperate mosquitoes [here](https://www.sciencedirect.com/science/article/pii/S0022519316300285), Code repo: [here](https://github.com/davewi13/Temperate-Mosquito-DDE)
3.  Paper: A novel approach for predicting risk of vector-borne disease establishment in marginal temperate environments under climate change: West Nile virus in the UK [here](https://doi.org/10.1098/rsif.2021.0049), Code repo: [here](https://github.com/davewi13/WNV_model)

The models in this package are run in discrete time, but the size of the time step (`dt` in code) may be arbitrarily small. The package implements deterministic and stochastic lifecycle models, as well as models with SEI infection dynamics in mosquitoes, where the force of infection upon mosquitoes (decomposed into `f`, `q`, and `kappa`; feeding rate, proportion of bloodmeals on each host species, and net infectiousness of each host species, respectively) is input from an external model. The use case is for this model of mosquito dynamics to be embedded in a larger simulation framework, perhaps including humans, birds, or other host species for the pathogen of interest. These models also implement movement of adult mosquitoes between patches based on a dispersal matrix.

Rates are converted into probabilities by an exponential CDF, which are interpreted as proportions of the population moving between states in the deterministic models and probabilities for Binomial draws of movement between states in the stochastic models.

# Installation

Installing this package requires [RcppArmadillo](https://dirk.eddelbuettel.com/code/rcpp.armadillo.html) to be installed on the users machine, which provides fast linear algebra. If you are having issues on MacOS, please read this guide ["R Compiler Tools for Rcpp on macOS"](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/).

You can install it either by cloning this repo and building from within RStudio, or directly via:

```
library(remotes)
remotes::install_github(repo = "dd-harp/culex-sim", dependencies = TRUE)
```

# Documentation

There are several vignettes in the package demonstrating use of the models.

# Simulation

To transform the delay differential equations into discrete time models, we note that
the backward looking delays can be transformed into forward looking queues of completion
times by first integrating the _o_rdinary differential equations for delay duration.

These ODEs are of the form:

```
$\dot{\tau} = 1 - \frac{g(t)}{g(t - \tau(t))}$
```

Where $\tau$ is the length of the maturation delay; that is, the time that individuals
maturing at time $t$ would have needed to wait before maturing (assuming survival). Because
we want to know the time an individual needs to wait before maturing if entering a stage
at time $t$, we first solve $\tau$ over a time horizon, discretize the solution to `dt`,
and subtract $t - \tau$, to get the duration of the forward looking delay.

![](man/figures/delays.pdf)
