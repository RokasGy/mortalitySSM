# mortalitySSM
State-space models for statistical mortality projections
The code was developed to perform stachastic Lee-Carter model mortality projections using state-space model set up. Two separate models were developed:
- basic linear model (DLM -Dynamic Linear Model)
- DLM with regime switching between low and high volatility regimes. 
Parameter fitting is performed using MCMC Gibbs sampler. The model useses R dlm package to perform Kalman filtering.
As model application, the code is provided for calculation of mortality VAR (Value-at-Risk).

As the imput the code uses mortality data obtained from Human Mortality Database https://mortality.org/

The following files are uploaded:
- **Lee-Carter DLM with switching 2022 03.R**  The code was used to derive mortality projections for the Swedish population in the article: https://www.mdpi.com/2227-7390/8/7/1053.
- **Lee-Carter DLM with switching v2 2022 07.R** Updated version of the model, which also includes stochastic modelling of alpha(x) parameters.
- **Mortality risk VAR model.R** Model used the calculate VAR rates. See article https://www.mdpi.com/2227-9091/7/2/58 for the description of the underlying methodology. 
- **Particle filter for likelihood estimation.R** The code used to estimate log-likelihood conditional on estimated parameters. See article: https://www.mdpi.com/2227-7390/8/7/1053 for the description of the methodology.
