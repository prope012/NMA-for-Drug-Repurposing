[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10069879.svg)](https://doi.org/10.5281/zenodo.10069879)

# Title of Paper:
Network Meta-Analysis to Predict the Efficacy of an Approved Treatment in a New Indication

# Purpose of these Files:
These files allow the user to replicate the simulation study in the above paper. They also demonstrate how to apply the proposed methods to a data set in practice. 

# Brief Description of these Files:
# 1. NMA-Simulation-Functions.R
These functions can be used to:
- Randomly generate a set of indication-specific basic parameters using a desired amount of correlation.
- Randomly generate a data set with binary outcomes according to a user-specified network.
- Fit a standard contrast-based NMA model to data from each indication.
- Fit the proposed NMA models for 2 indications to a given data set.
- Compute summary statistics (Bias, variance, RMSE, average width of the 95% credible interval, and coverage probability).
- Compute probability of success for a future clinical trial.  
  
The purpose, inputs, and outputs of each function are included in this R script.

# 2. NMA-Simulation-Example.R:
This R script demonstrates how to use the functions in NMA-Simulation-Functions.R to replicate the simulation study described in Section 5.1 of the paper. In particular, the R object "Mixed.Mod.Trt2.Corr5" will contain the posterior samples for $d^2_12$ arising from the mixed effects model with a half-t prior, coefficient set 1, and a correlation of 0.5. Then, the R object "Mixed.Mod.Trt2.Corr5.Ht.Stats" will contain the following summary statistics for the posterior median of $d^2_12$: bias, variance, RMSE, expected width of the 95% credible interval, and the coverage probability of the 95% credible interval. These summary statistics were used to create Figures 3, 4, and 5 in this manuscript. All simulation results can be replicated by altering the input to the functions used to generate the R object "Mixed.Mod.Trt2.Corr5."

In addition, this file demonstrates how to apply the proposed methodology in practice, including how to compute PoS for our Motivating Example using the file pso-psa-data.csv [Section 6 of the manuscript]. The data contained in pso-psa-data.csv was curated from publicly available data for multi-arm trials in psoriasis and psoriatic arthritis. For more information, please see Section 2 of the manuscript.

# 3. "JAGS CODE" Folder:
This folder contains the Jags scripts required to implement each standard CB-NMA model and each proposed NMA model for 2 indications. 

# 4. example-data.csv:
This CSV file contains an artifical dataset that is suitable for this methodology. This dataset was used in our simulation study.

# 5. pso-psa-data.csv:
This CSV file contains the dataset that was used in our motivating example. This dataset has been de-identified for confidentiality reasons.

# 6. renv.lock:
This is a lock file produced by the Rpackage 'renv' that allows users to run our code using the same package versions.
