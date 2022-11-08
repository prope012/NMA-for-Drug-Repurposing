# Title of Paper:
Network Meta-Analysis to Predict the Efficacy of an Approved Treatment in a New Indication

# Purpose of these Files:
These files allow the user to replicate the simulation study in the above paper. They also demonstrate how to apply the proposed methods to a data set in practice. 

# Brief Description of these Files:
# 1. NMA-Simulation-Functions.R
These functions can be used to randomly generate a set of indication-specific basic parameters with a certain amount of correlation, randomly generate a data set with binary outcomes according to a user-specified network, fit a standard contrast-based NMA model to data from each indication, fit each proposed NMA model for 2 indications, compute summary statistics if performing a simulation study (bias, variance, RMSE, average width of the 95% credible interval, and coverage probability), and compute probability of success. 

# 2. NMA-Simulation-Example.R:
This R script demonstrates how to use the functions in NMA-Simualtion-Functions.R to replicate our simulation study. They also demonstrate how to apply the proposed methodology in practice, including how to compute PoS. 

# 3. Jags Scripts Folder:
This folder contains the Jags scripts required to implement each standard CB-NMA model and each proposed NMA model for 2 indications. 
