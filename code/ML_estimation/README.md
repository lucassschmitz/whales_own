This folder contains the code for the Maximum Likelihood Estimation of a production model. 
Files that start with 'est' are for estimating the model in the data and the ones starting with 'sim' are to estimate the model in simulated data. The ones that have 'ind' in the name estimate the model assuming independent errors and the ones with 'corr' allow for correlation in the errors. For example, `est\_corr.m` estimates the model using the real data and allows for correlation. 

`est\_corr.m` is the main script.

**Folder Contents**

-Estimation Scripts: `est\_ind.m`, `est\_corr.m`, `sim\_corr.m`, `sim\_ind.m`

-Likelihood Functions: `voyageLik.m`, `captainLik.m`, `globalLik.m`, `L\_v\_corr.m`, `L\_c\_corr.m`, `L\_c\_corr\_int.m`, `globalLik\_corr.m`

-Data Simulation Functions: `simulate\_data.m`,  `simulate\_partial\_data.m`

-Gauss-Hermite Quadrature Functions: `HG.m`, `HG2D.m`, `HG3D.m`

-Parameter Transformation Functions: `to\_theta.m`, `to\_chol\_theta.m`, `from\_theta.m`



**Files Content:**

* `est\_corr.m`: The main script for estimating the production model on \*\*real\*\* data, allowing for \*\*correlation\*\* in the errors between different products.
* `sim\_corr.m` and `sim\_ind.m`: Script for estimating the production model on simulated data, with and without allowing for \*\*correlation\*\* in the errors. Used to test the model and estimation procedure.
* `est\_ind.m`: Script for estimating the production model assuming \*\*independent\*\* errors across products. It loads the data, defines the model parameters, and uses `fminunc` to find the maximum likelihood estimates.
* `voyageLik.m` \[function]\*\*: calculates the likelihood for a single voyage, assuming independent errors. 
* `captainLik.m` \[function]\*\*: This function computes the likelihood for a single captain by integrating the voyage likelihoods over the captain's unobserved skill (random effect), assuming independent errors. 
* `globalLik.m` \[function]\*\*: This function calculates the total log-likelihood for the entire dataset by summing the log-likelihoods of each captain, assuming independent errors.
* `L\_v\_corr.m` \[function]\*\*: A function that calculates the voyage-level likelihood, allowing for \*\*correlation\*\* in the shocks across different products.
* `L\_c\_corr.m` \[function]: calculates the captain-level likelihood conditional on a specific skill level, allowing for \*\*correlation\*\*.
* `L\_c\_corr\_int.m` \[function]: integrates the captain-level likelihood over the distribution of captain skills to get the unconditional likelihood for a captain, allowing for \*\*correlation\*\*.
* `globalLik\_corr.m` \[function]: This function calculates the total log-likelihood for the entire dataset, allowing for \*\*correlation\*\* in the errors.  
* `simulate\_data.m`,  `simulate\_partial\_data.m` \[functions]\*\*: generate simulated data for testing the estimation code. 'simulate\_data' simulates the data from scratch whereas 'simulate\_partial\_data' uses duration and weights of real voyages to simulate the remaining data. 
* `HG.m`, `HG2D.m`, `HG3D.m` \[functions]\*\*:  generate nodes and weights for one, two, and three-dimensional Gauss-Hermite quadrature. Used for integration. 

\* \*\*`to\_theta.m`, `from\_theta.m`, `to\_chol\_theta.m` \[functions]\*\*: These are helper functions to manage the parameter vector `theta`. They allow to use the whole covariance matrix in the `theta` vector or using the parameter vector only with the Cholesky decomposition of the covariance matrix. 