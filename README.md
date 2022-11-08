# Title of Paper:
Network Meta-Analysis to Predict the Efficacy of an Approved Treatment in a New Indication

# Purpose of these Files:
These files allow the user to replicate the simulation study in the above paper. They also demonstrate how to apply the proposed methods to a data set in practice. 

# Brief Description of these Files:
# 1. NMA-Simulation-Functions.R
These functions can be used to:
- Randomly generate a set of indication-specific basic parameters using a desired amount of correlation
- Randomly generate a data set with binary outcomes according to a user-specified network
- Fit a standard contrast-based NMA model to data from each indication
- Fit the proposed NMA models for 2 indications to a given data set
- Compute summary statistics (Bias, variance, RMSE, average width of the 95% credible interval, and coverage probability)
- Compute probability of success for a future clinical trial  
  
The purpose, inputs, and outputs of each function are included in this R script.

# 2. NMA-Simulation-Example.R:
This R script demonstrates how to use the functions in NMA-Simualtion-Functions.R to replicate our simulation study. They also demonstrate how to apply the proposed methodology in practice, including how to compute PoS. 

# 3. Jags Scripts Folder:
This folder contains the Jags scripts required to implement each standard CB-NMA model and each proposed NMA model for 2 indications. 
