# Inferring the causes of animal social network structure from time-series data
This project models dyadic social interactions using generative simulations and Bayesian inference.
We first encode a generative model in R: a  Continuous-Time Markov Chain on top of which we layer focal-animal sampling.
We then develop a statistical model that takes the simulated data as input and returns posterior probability distributions for the parameters of the generative model.

## Structure of the repo
The core of this repository consists of Quarto notebooks (`.qmd` files) that implement both the generative and inferential models.
Each notebook can be read as an `.html` file in any web browser (see hyperlinks below), or opened and executed in _R_ as a `.qmd` document.
These notebooks are self-contained and can be run independently to reproduce each step of the workflow.

-   [`01.1_core_model_sim`](https://htmlpreview.github.io/?https://github.com/BenKawam/ASN-time-series/blob/main/01.1_core_model_sim.html). Core generative and statistical models. Includes validation of the statistical model using the generative simulation, prior and posterior predictive checks.
-   [`01.2_core_model_empirical`](https://htmlpreview.github.io/?https://github.com/BenKawam/ASN-time-series/blob/main/01.2_core_model_empirical.html). We update the statistical model from `01.1_core_model_sim` with empirical data. Includes computation of behavioural expectations and posterior predictive checks.
-   [`02.1_extended_model_sim`](https://htmlpreview.github.io/?https://github.com/BenKawam/ASN-time-series/blob/main/02.1_extended_model_sim.html). Extended generative and statistical models. Includes validation of the statistical model using the generative simulation and posterior predictive checks.
-   [`02.2_extended_model_empirical`](https://htmlpreview.github.io/?https://github.com/BenKawam/ASN-time-series/blob/main/02.2_extended_model_empirical.html). We update the statistical model from `02.1_extended_model_sim` with empirical data. Includes computation of average causal effects and posterior predictive checks.

The notebooks use and produce data in the following folders:
- `empirical_data`: Focal-animal sampling data collected at Phu Khieo Wildlife Sanctuary between June 2017 and July 2018.
- `fitted_models`: MCMC samples (obtaining by running the notebooks; not tracked).
- `functions`: functions used in the Quarto notebooks.
- `sim_data`: Syhtnetic data generated in the notebooks (obtaining by running the notebooks; tracked).
- `stan_models`: stan models called in the the notebooks.
