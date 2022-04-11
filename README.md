# culex

A discretized stochastic and deterministic model of temperature and photoperiod driven culex dynamics with movement of adults based on the mathematical model(s) presented in:

1. Paper: Uncovering mechanisms behind mosquito seasonality by integrating mathematical models and daily empirical population data: Culex pipiens in the UK [here](https://parasitesandvectors.biomedcentral.com/articles/10.1186/s13071-019-3321-2), Code repo: [here](https://github.com/davewi13/Mosquito-seasonality-paper)
2. Paper: Modelling the effect of temperature on the seasonal population dynamics of temperate mosquitoes [here](https://www.sciencedirect.com/science/article/pii/S0022519316300285), Code repo: [here](https://github.com/davewi13/Temperate-Mosquito-DDE)
3. Paper: A novel approach for predicting risk of vector-borne disease establishment in marginal temperate environments under climate change: West Nile virus in the UK [here](https://doi.org/10.1098/rsif.2021.0049), Code repo: https://github.com/davewi13/WNV_model

Running this model requires [RcppArmadillo](https://dirk.eddelbuettel.com/code/rcpp.armadillo.html) to be installed on the users machine, which provides fast linear algebra. The model is run in discrete time, but the size of the time step (`dt` in code) may be arbitrarily small. Rates are converted into probabilities by an exponential CDF.