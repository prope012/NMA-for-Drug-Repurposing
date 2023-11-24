###################################################################
#### Functions to replicate the simulation study in the paper: ####
###### "NMA to predict the efficacy of an approved treatment ######
##################### in a new indication" ########################
###################################################################

#renv::init()
#------------------------------------------------------#
#-------------- Load Required Packages ----------------#
#------------------------------------------------------#
library(R2jags)
library(dplyr)
library(MASS)
library(locfit)
library(mcmcplots)
library(ggthemes)
library(ggplot2)

# specify directory for jags scripts
work.dir <- "./"
jags.script.dir <- paste(work.dir, "JAGS-Scripts/",sep="")

#------------------------------------------------------#
#---------- Generate Correlated Basic Parameters ------#
#------------------------------------------------------#

### Purpose: This function generates correlated basic parameters (LORs relative to placebo) across indications using a linear mixed model of the following form:
#            LORs = beta0 + beta1*I{class = 2} + beta2*I{indication = 2} + gamma_k + Sigma_{ki}.
#            Using these randomly generated values, the values of the LORs comparing the active treatment to the baseline treatment in each study in the network diagram in Figure 2 are determined.
#            Baseline effect parameters are also generated according to the network diagram in Figure 2.

### Inputs:
#   1) beta0: The value of the intercept in the linear mixed model.
#   2) beta1: The value of the fixed class effect in the linear mixed model.
#   3) beta2: The value of the fixed indication effect in the linear mixed model.
#   4) sigma.g: The standard deviation of gamma_k (the random treatment effect where gamma_k are iid N(0,sigma.g^2)).
#   5) sigma.e: The standard deviation of Sigma_{ki} (the random error term where Sigma_{ki} are idd N(0,sigma.e^2)).
#   6) class.ind: A 1x8 vector containing the drug class of each active treatment. 
#             e.g., class = c(1,1,2,2) indicates that active treatments 1 and 2 are in class 1 and 3 and 4 are in class 2.
#   7) heterogeneous.CE = an indicator that equals 1 if a heterogeneous class effect is desired and 0 otherwise.

### Outputs: A list with 5 components:
#   1) A 8x2 data.frame containing the randomly generated values for the basic parameters. 
#      Column corresponds to indication and Row corresponds to treatment (e.g., column 1 = indication 1 and row 2 = the LOR comparing treatment 2 vs. placebo).
#   2) A 40x2 matrix containing the true values of the LORs comparing the treatment in each arm to the baseline treatment for indication 1. 
#      Rows represent studies and columns represent study arms.
#   3) A 35x2 matrix containing the true values of the LORs comparing the treatment in each arm to the baseline treatment for indication 2. 
#      Rows represent studies and columns represent study arms.
#   4) The baseline effects for the studies in indication 1.
#   5) The baseline effects for the studies in indication 2.

sim.pars <- function(beta0, beta1, beta2, sigma.g, sigma.e, class.ind, heterogeneous.CE = F){
  
  ### draw random error term for each treatment and indication
  err1 <- rnorm(8, 0, sd = sigma.e); err2 <- rnorm(8, 0, sd = sigma.e)
  
  ### draw random treatment effect for each treatment
  gam <- rnorm(8, 0, sd = sigma.g)
  
  ### find LORs for indication 1 and indication 2; add 0 for the 1 v 1 comparison
  d1 <- round(c(0, beta0 + beta1*(class.ind-1) + gam + err1),1)
  if(heterogeneous.CE == T){
      d2 <- round(c(0, beta0 + beta2 + gam + err2),1)} else{
        d2 <- round(c(0, beta0 + beta1*(class.ind-1) + beta2 + gam + err2),1)}
  
  ### Establish matrices containing the true values of the population-averaged LORs comparing the treatment in each arm to the baseline treatment
  #   Note that column 1 should contain all 0s.
  d1_vals <- matrix(c(rep(d1[1],40), rep(d1[2],8), rep(d1[3],4), rep(d1[4], 4),
                      rep(d1[5],2), rep(d1[6],5), rep(d1[7],3), rep(d1[9],2),
                      rep(d1[8]-d1[2],3), rep(d1[4]-d1[3],3), rep(d1[6]-d1[5],3), 
                      rep(d1[5]-d1[2],2), rep(d1[8]-d1[7],1)), ncol=2)
  
  d2_vals <- matrix(c(rep(d2[1],35), rep(d2[2],3), rep(d2[3],2), rep(d2[4],2),
                      rep(d2[5],6), rep(d2[7],3), rep(d2[8],5),  rep(d2[9],6),
                      rep(d2[6]-d2[3],2), rep(d2[8]-d2[4],3), rep(d2[7]-d2[5],3)), ncol = 2)
  
  ### Find reasonable baseline values assuming the mean log(odds) = -3 for placebo in indication 1 and -2.3 for placebo in indication 2
  #   Note that log(OR k vs. b) = log(odds k) - log(odds b)
  log.odds.ind1 <- d1-3; log.odds.ind2 <- d2-2.3
  mu1 <- round(c(rnorm(28,log.odds.ind1[1], 0.25), rnorm(3, log.odds.ind1[2], 0.25), rnorm(3, log.odds.ind1[3], 0.25),
                 rnorm(3, log.odds.ind1[5], 0.25), rnorm(2, log.odds.ind1[2], 0.25), rnorm(1, log.odds.ind1[7], 0.25)),1)
  mu2 <- round(c(rnorm(27,log.odds.ind2[1], 0.25), rnorm(2, log.odds.ind2[3], 0.25), rnorm(3, log.odds.ind2[4], 0.25),
                 rnorm(3, log.odds.ind2[5], 0.25)),1)
  
  ### return mu1, mu2, d_bk1, d_bk2
  return(list(d = data.frame(d1=d1,d2=d2), d1_vals = d1_vals, d2_vals = d2_vals, mu1 = mu1, mu2 = mu2))
}


#------------------------------------------------------#
#------------------- Generate Data --------------------#
#------------------------------------------------------#

### Purpose: This function generates binary outcome data for 2 indications that will be analyzed by an NMA model.

### Inputs: 
#   1) J: The total number of studies in the data set.
#   2) num.arms: A 1xJ vector containing the number of arms in each study.
#   3) indication: A 1xJ vector containing the indication associated with each study (2 for indication 2, 1 for indication 1).
#   4) trt: A J by max(num.arms) dataframe containing the treatments studied in each arm of each study.
#           Rows represent studies and columns represent study arms. Treatments should be assigned a numeric number with 1 denoting the network reference treatment.
#           The baseline treatment for each study should be in arm 1. NA should be used for studies having less than max(num.arms) arms.
#   5) n: A J by max(num.arms) dataframe containing the number of participants in each arm of each study. 
#         Rows represent studies and columns represent study arms. NA should be used for studies having less than max(num.arms) arms.
#   6) d: A J by max(num.arms) dataframe containing the true values of the population-averaged LORs comparing the treatment in each arm to the baseline treatment.
#         Rows represent studies and columns represent study arms. NA should be used for studies having less than max(num.arms) arms.
#   7) mu: A 1xJ vector containing the true values of the log-odds of response for the baseline treatment in each study. 
#   8) sigma: A 1x2 vector containing the true values of the between-study standard deviations for indications 1 and 2. Default is c(0.25,0.25). 

### Outputs: 
#   1) A dataframe with the following variables: study number, indication, treatment (trt), number of successes (r), and number of participants (n) in each arm of each study.

get.dat <- function(J, num.arms, indication, trt, n, d, mu, sigma = c(0.25,0.25)){
  
  ### Format data set as desired in long format
  dat <- data.frame(study = rep(1:J,num.arms),
                    indication = rep(indication, num.arms),
                    trt = as.numeric(na.omit(as.vector(t(trt)))),
                    n = as.numeric(na.omit(as.vector(t(n)))),
                    mu = rep(mu,num.arms),
                    d = as.numeric(na.omit(as.vector(t(d))))
  )
  
  ### Simulate study-specific effects using delta_jbk ~ N(d_bki, sigma_i^2)
  for(i in 1:nrow(dat)){
    if(dat$d[i] != 0 & dat$indication[i] == 1){dat$d[i] <- rnorm(1,dat$d[i],sigma[1])}
    if(dat$d[i] != 0 & dat$indication[i] == 2){dat$d[i] <- rnorm(1,dat$d[i],sigma[2])}
  }
  
  ### Solve for the probability of response in each arm in each trial (i.e. p_jk)
  #   Then, simulate the number of successes in each arm assuming r_jk ~ Bin(n_jk,p_jk)
  dat <- dat %>% mutate(p = expit(mu + d),
                        r = rbinom(n = nrow(dat),size = n, prob = p)) 
  
  ### Grab variables of interest and return data set
  return(dat %>% dplyr::select(study, indication, trt, n, r))
}


#------------------------------------------------------#
#------------- Fit Standard NMA Model  ----------------#
#------------------------------------------------------#

### Purpose: This function fits a standard random effects NMA model to a data set with ONE indication. 

### Inputs:
#   1) dat: A dataframe with the format outputted by the "get.dat" function but for only one indication. 
#   2) sd.prior: The prior distribution for the between-study standard deviation parameter in the random effects model.
#                Users can choose from the following 4 options: "Unif" (Uniform(0,5)), "HN" (standard half-normal), "Ht" (half-t with 7 degrees of freedom and a scale of 2.5), and "LN" (log-normal(2.70,1.52)).

### Output: A dataframe containing the estimated mean, sd, and 2.5th, 25th, 50th, 75th, and 97.5th quantiles for all model parameters.
fit.std.mod <- function(dat, sd.prior = c("unif","HN","Ht","LN")){
  
  ### Prepare data for NMA
  ns <- length(unique(dat$study)) # number of studies
  na <- (dat %>% group_by(study) %>% dplyr::summarize(na = length(study)))$na # number of arms per study
  r <- NULL; for(i in 1:max(na)){r <- cbind(r, (dat %>% group_by(study) %>% summarise(r=r[i]))$r)} # wide format for number of successes
  n <- NULL; for(i in 1:max(na)){n <- cbind(n, (dat %>% group_by(study) %>% summarise(n=n[i]))$n)} # wide format for number participants
  t <- NULL; for(i in 1:max(na)){t <- cbind(t, (dat %>% group_by(study) %>% summarise(t=trt[i]))$t)} # wide format for treatments
  nt <- max(dat$trt) # number of treatments in the network
  
  ### Identify appropriate jags script based on type of model and desired prior
  jags.script <- case_when(sd.prior == "unif" ~ "Standard Models/std-random-effects-unif.txt",
                           sd.prior == "HN" ~ "Standard Models/std-random-effects-HN.txt",
                           sd.prior == "Ht" ~ "Standard Models/std-random-effects-Ht.txt",
                           sd.prior == "LN" ~ "Standard Models/std-random-effects-LN.txt")
  
  ### Fit NMA model
  params = c("d","sd")
  data <- list(ns = ns, r = r, t = t, n = n, nt = nt, na = na)
  jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                 n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""))
  
  ### Grab desired output information
  foo <- jags.fit$BUGSoutput$summary
  foo_return <- data.frame(sd = foo[,"sd"], "2.5%" = foo[,"2.5%"], Median = foo[,"50%"],
                           "97.5%" = foo[,"97.5%"], Rhat = foo[,"Rhat"])
  rm(jags.fit)
  return(foo_return)
}


#------------------------------------------------------#
#------ Fit Proposed NMA Models for 2 Indications -----#
#------------------------------------------------------#

### Purpose: This function fits our proposed NMA models to a data set with TWO indications. 

### Inputs:
#   1) dat: A dataframe with the format outputted by the "get.dat" function.
#   2) sd.prior: The prior distribution for the between-study standard deviation parameter in the random effects model.
#                Users can choose from the following 4 options: "Unif" (Uniform(0,5)), "HN" (standard half-normal), "Ht" (half-t with 7 degrees of freedom and a scale of 2.5), and "LN" (log-normal(2.70,1.52)).
#   3) class: A 1xK vector containing the drug class of each treatment. It is assumed that the network reference treatment is placebo or standard of care with class[1] = NA.
#             e.g., class = c(NA,1,1,2,2) indicates that treatments 2 and 3 are in class 1 and 4 and 5 are in class 2.
#   4) mod.type: Indicates what two-indication model to fit. There are 9 options: 
#               - "bvn" for the bivariate normal model, "bvn-class" for the bivariate normal model with a homogeneous class effect, and "bvn-int" for the bivariate normal model with a heterogeneous class effect
#               - "mvn" for the multivariate normal model, "mvn-class" for the multivariate normal model with a homogeneous class effect, and "mvn-int" for the multivariate normal model with a heterogeneous class effect
#               - "mixed" for the mixed effects model, "mixed-class" for the mixed effects model with a homogeneous class effect, and "mixed-int" for the mixed effects model with a heterogeneous class effect

### Outputs: A dataframe containing:
#     1) Posterior SD & quantiles for all relevant model parameters.
#     2) The Gelman-Rubin convergence diagnostic (Rhat) for each parameter. 

fit.2ind.mod <- function(dat, random = F, sd.prior = c("unif","HN","Ht","LN"), class = NULL,
                         mod.type = c("bvn", "bvn-class", "bvn-int", "mvn", "mvn-class", "mvn-int", "mixed", "mixed-class", "mixed-int")){
  
  ### Identify appropriate jags script based on model type and prior
  jags.script <- paste(mod.type,"-random-effects-",sd.prior,".txt",sep="")

  ### Prepare data for NMA
  ns <- length(unique(dat$study)) # number of studies
  na <- (dat %>% group_by(study) %>% dplyr::summarize(na = length(study)))$na # number of arms per study
  r <- NULL; for(i in 1:max(na)){r <- cbind(r, (dat %>% group_by(study) %>% summarise(r=r[i]))$r)} # wide format for number of successes
  n <- NULL; for(i in 1:max(na)){n <- cbind(n, (dat %>% group_by(study) %>% summarise(n=n[i]))$n)} # wide format for number participants
  t <- NULL; for(i in 1:max(na)){t <- cbind(t, (dat %>% group_by(study) %>% summarise(t=trt[i]))$t)} # wide format for treatments
  nt <- max(dat$trt) # number of treatments in the network 
  indication <- (dat %>% group_by(study) %>% dplyr::summarize(ind = indication[1]))$ind # indication indicator
  
  ### Fit desired model:
  if(mod.type == "bvn"){
    params = c("d","beta0","beta2","sd","rho","sdd","d.new")
    init.vals <- list(d = matrix(rep(c(NA, rep(0, nt-1)),2), ncol = 2))
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt, indication=indication, na = na)
    jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  if(mod.type == "bvn-class"){
    params = c("d","beta0","beta1","beta2","sd","rho","sdd","d.new")
    init.vals <- list(d = matrix(rep(c(NA, rep(0, nt-1)),2), ncol = 2))
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt, indication=indication, na = na, class = class)
    jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  if(mod.type == "bvn-int"){
    params = c("d","beta0","beta1","beta2","beta3","sd","rho","sdd","d.new")
    init.vals <- list(d = matrix(rep(c(NA, rep(0, nt-1)),2), ncol = 2))
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt, indication=indication, na = na, class = class)
    jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  if(mod.type == "mvn"){
    params = c("d","sd","rho","sdd","d.new","beta0","beta2")
    init.vals <- list(d = c(NA,NA,rep(0,(nt-1)*2)))
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt,indication=indication, na = na)
    jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  if(mod.type == "mvn-class"){
    params = c("d","sd","rho","sdd","d.new","beta0","beta2","beta1")
    init.vals <- list(d = c(NA,NA,rep(0,(nt-1)*2)))
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt,indication=indication, na = na, class = class)
    jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  if(mod.type == "mvn-int"){
    params = c("d","sd","rho","sdd","d.new","beta0","beta2","beta1")
    init.vals <- list(d = c(NA,NA,rep(0,(nt-1)*2)))
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt,indication=indication, na = na, class = class)
    jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  if(mod.type == "mixed"){
    params = c("d","beta0","beta2","gamma","sdg","sde","sdd","d.new")
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt, indication=indication, na = na)
    jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  if(mod.type == "mixed-class"){
    params = c("d","beta0","beta1","beta2","gamma","sdg","sde","sdd","d.new")
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt, indication=indication, na = na, class = class)
    jags.fit =jags(data = data, n.burnin =5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  if(mod.type == "mixed-int"){
    params = c("d","beta0","beta1","beta2","beta3","gamma","sdg","sde","sdd","d.new")
    data <- list(ns = ns, r = r, t = t, n = n, nt=nt, indication=indication, na = na, class=class)
    jags.fit =jags(data = data, n.burnin =5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
                   n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
                   quiet = TRUE, progress.bar = "none")
  }
  ### Grab posterior percentiles of interest & Rhat
  a <- jags.fit$BUGSoutput$sims.matrix
  b <- data.frame(t(apply(a,2,function(x) quantile(x, probs = c(0.025,0.05,0.075,0.25,0.50,0.75,0.925,0.95,0.975)))))
  b$sd <- apply(a, 2, sd)
  b$Rhat <- jags.fit$BUGSoutput$summary[,"Rhat"]
  
  ### Grab desired output information
  foo_return <- data.frame(sd = b[,"sd"], "2.5%" = b[,"X2.5."], "5%" = b[,"X5."], "7.5%" = b[,"X7.5."], "25%" = b[,"X25."], Median = b[,"X50."], "75%" = b[,"X75."], 
                           "92.5%" = b[,"X92.5."], "95%" = b[,"X95."],"97.5%" = b[,"X97.5."], Rhat = b[,"Rhat"])
  rownames(foo_return) <- colnames(a)
  rm(jags.fit)
  return(foo_return)
}


#------------------------------------------------------#
#------------- Compute Summary Statistics -------------#
#------------------------------------------------------#

### Purpose: The purpose of this function is to compute the bias, variance, and RMSE of the posterior median for each basic parameter,
#            as well as the expected width and coverage probability of the 95% credible interval. 
#            Posterior convergence is also assessed by finding the proportion of simulations where the gelman-rubin statistic <1.1 for each model parameter.

### Inputs:
#   1) sim.list: A list of simulation results arising from the "fit.std.mod" or "fit.2ind.mod" functions
#   2) true values: a vector containing the true values of the basic parameters for treatments 1 to K. 
#                   For two indications, begin with the true values for indication 1 (e.g., c(d111,d121,...,d1K1,d112,d122,...,d1k2))
#   3) nt: the number of treatments in the list.
#   4) inds: denotes what indication(s) this list pertains to; "Ind1" = indication 1, "Ind2" = indication 2, and "Both" = both indications 1 and 2.
#   5) mod.type: The type of model used in the simulation. This can be equal to one of the two-indication models or "standard" for a standard NMA model for one indication. 
#   6) LOOCV: TRUE or FALSE; TRUE if sim.list corresponds to LOOCV and FALSE otherwise.
#   7) lower.CrI: the lower quantile for the desired credible interval. Options: X2.5., X5., X10., X15., X25., X50., X75., X85., X90., X95., X97.5.
#   8) upper.CrI: the upper quantile for the desired credible interval. Options: X2.5., X5., X10., X15., X25., X50., X75., X85., X90., X95., X97.5.

### Outputs: A dataframe containing the bias, variance, and RMSE of the posterior median for each basic parameter.
#            expected width and coverage probability of the 95% credible interval is also provided.

get.summary.stats <- function(sim.list, true.vals, nt, inds = c("Ind1","Ind2","Both"),
                              mod.type = c("standard","bvn", "bvn-class", "bvn-int", "mvn", "mvn-class", "mvn-int", "mixed", "mixed-class", "mixed-int"),
                              lower.CrI = "X2.5.", upper.CrI = "X97.5."){
  
  ### grab basic parameters from sim.list
  d.pars <- lapply(sim.list, function(x) x[grep("^d\\[", rownames(x)),])
  
  ### reformat d.pars for mvn models
  if((mod.type %in% c("mvn","mvn-class"))){
    ind <- c(seq(1,nt*2-1,2),seq(2,2*nt,2))
    d.pars <- lapply(d.pars, function(x) x[ind,])
  }
  
  ### Grab posterior medians (which we approximated using MCMC methods)
  med.est <- do.call(rbind,lapply(d.pars, function(x) x[,"Median"]))
  
  ### Compute bias, variance, RMSE using the posterior median
  bias <- colMeans(med.est) - true.vals
  var <- apply(med.est, 2, var)
  RMSE <- sapply(1:length(true.vals), function(x) sqrt(mean((med.est[,x]-true.vals[x])^2)))
  
  ### Compute coverage probability and average credible interval width using desired percentiles
  CP <- colMeans(do.call(rbind,lapply(d.pars, function(x) true.vals >= x[,lower.CrI] & true.vals <= x[,upper.CrI])))
  CrI.width <- colMeans(do.call(rbind,lapply(d.pars, function(x) x[,upper.CrI] - x[,lower.CrI])))
  
  ### Compute proportion of simulations where Rhat < 1.1.
  prop.Rhat <- colMeans(do.call(rbind,lapply(d.pars, function(x) x[,"Rhat"] < 1.1)))
  
  ### Format Results
  if(inds == "Ind1"){names = paste("d",1:nt,sep = "",1)}
  if(inds == "Ind2"){names = paste("d",1:nt,sep = "",2)}
  if(inds == "Both"){names = c(paste("d",1:nt,sep = "",1),paste("d",1:nt,sep = "",2))}
  
  res <- data.frame(Bias = bias, Variance = var, RMSE = RMSE, CovProb = CP, CrIWidth = CrI.width, Rhat = prop.Rhat)
  rownames(res) <- names

  return(res)
}


#------------------------------------------------------#
#---- Get Data for Leave-One-Out Cross-Validation -----#
#------------------------------------------------------#

### Purpose: The purpose of this function is to generate data that can be used for LOOCV.
#            Data sets are generated after sequentially removing treatments from the network for indication 2 (by removing all indication 2 studies that contain a given treatment).

### Inputs: 
#   1) J: The total number of studies in the COMPLETE data set.
#   2) num.arms: A 1xJ vector containing the number of arms in each study in the COMPLETE data set.
#   3) indication: A 1xJ vector containing the indication associated with each study (2 for indication 2, 1 for indication 1) in the COMPLETE data set.
#   4) trt: A J by max(num.arms) dataframe containing the treatments studied in each arm of each study in the COMPLETE data set.
#           Rows represent studies and columns represent study arms. Treatments should be assigned a numeric number with 1 denoting the network reference treatment.
#           The baseline treatment for each study should be in arm 1. NA should be used for studies having less than max(num.arms) arms.
#   5) n: A J by max(num.arms) dataframe containing the number of participants in each arm of each study in the COMPLETE data set. 
#         Rows represent studies and columns represent study arms. NA should be used for studies having less than max(num.arms) arms.
#   6) d: A J by max(num.arms) dataframe containing the true values of the population-averaged log-odds ratios comparing the treatment in each arm to the baseline treatment.
#         Rows represent studies and columns represent study arms. NA should be used for studies having less than max(num.arms) arms.
#   7) mu: A 1xJ vector containing the true values of the log-odds of response for the baseline treatment in each study in the COMPLETE data set.
#   8) sigma: A 1x2 vector containing the true values of the between-study standard deviations for indications 1 and 2. Default is c(0.25,0.25). 
#   9) niters: The number of data sets to generate for each LOOCV scenario. 

### Outputs: A list containing a list of niters data sets for each LOOCV scenario.

get.LOOCV.dat <- function(J, num.arms, indication, trt, n, d, mu, sigma = c(NA,NA), niters){
  
  ### Get vector of treatments (excluding reference treatment)
  nt <- max(as.numeric(na.omit(as.vector(t(trt)))))
  
  ### Get data sets after excluding 1 treatment in indication 2; repeat for all treatments
  foo <- lapply(2:nt, function(k){
    
    ### identify rows with treatment k in indication 2; simulate data without treatment k in indication 2
    ind <- apply(trt, 1, function(x) k %in% x) & indication == 2
    
    ### generate niters data sets
    dat.sub <- lapply(1:niters, function(x){
      
      set.seed(x); get.dat(J = J-sum(ind), num.arms = num.arms[!ind], indication = indication[!ind], trt = trt[!ind,], n = n[!ind,], 
                           d = d[!ind,], mu = mu[!ind], sigma = sigma)})
    return(dat.sub)})
  
  ### return LOOCV data sets
  return(setNames(foo, paste("Minus Trt",2:nt)))
}


#------------------------------------------------------#
#-------------------- Compute PoS ---------------------#
#------------------------------------------------------#

### Purpose: The purpose of this function is to estimate the probability of success (PoS) for a candidate treatment.
#            A trial is deemed successful if it achieves statistical significance (two-sided p-value <= alpha) and clinical significance (observed treatment effect > some threshold).
#            The primary analysis uses fisher's exact test. 

### Inputs:
#   1) fitted.jags.mod: A fitted jags model from the "fit.2ind.mod" function. 
#   2) trt.num: the number corresponding to the candidate treatment.
#   3) mvn: An indicator for whether a mvn model was fit (TRUE for mvn and FALSE otherwise).
#   4) nt: The total number of treatments analyzed by fitted.jags.mod.
#   5) piC: The probability of response for the reference treatment obtained from external cohort studies or expert opinion for indication 2.
#   6) n: The sample size per arm.
#   7) pow: The desired power of the design. Default is 0.90.
#   8) two.sided.TIE: the two-sided type I error rate of the design. Default is 0.05.
#   9) niters: Number of posterior predictive samples to use when approximating PoS via MCMC integration. Default is 10000.
#   10) delta: a clinically meaningful log-odds ratio; that is, the efficacy threshold that you wish to surpass in the future study.

### Outputs: The PoS for the candidate treatment.

get.pos <- function(fitted.jags.mod, trt.num, mvn = FALSE, nt, piC, n, two.sided.TIE = 0.05, delta, niters = 10000){
  
  ### Grab samples from the posterior predictive distribution for the candidate treatment.
  if(mvn == F){
    post.samps.deltaRR <- sapply(fitted.jags.mod$BUGSoutput$sims.list$d.new[,,2][,trt.num], function(x) (exp(x)*(piC/(1-piC)))/(1+(exp(x)*(piC/(1-piC)))))[1:niters]
  } else{
    post.samps.deltaRR <- sapply(fitted.jags.mod$BUGSoutput$sims.list$d.new[,seq(2,nt*2,2)][,trt.num], function(x) (exp(x)*(piC/(1-piC)))/(1+(exp(x)*(piC/(1-piC)))))[1:niters]
  }
  
  ### Simulate trial using n, piC, and each piE_hat
  RR <- sapply(post.samps.deltaRR, function(x){
    trial.dat <- data.frame(y = c(rbinom(n, 1, piC), rbinom(n, 1, x)), grp = c(rep(0,n),rep(1,n)))
    trial.dat$y <- factor(trial.dat$y, levels = c(0,1)); trial.dat$grp <- factor(trial.dat$grp, levels = c(0,1))
    rr <- tapply(trial.dat$y==1, trial.dat$grp, mean)
    pval <- fisher.test(matrix(table(trial.dat),ncol = 2))$p.value
    return(c(rr,pval))})
  
  # transform to LOR
  LOR <- log(RR[2,]/(1-RR[2,])) - log(RR[1,]/(1-RR[1,]))
  
  # determine if pval < 0.05
  StatSig <- RR[3,] < 0.05
  
  # determine if clinically significant
  ClinSig <- LOR > delta
  
  ### Pr(Stat Sig + Clin Sig)
  return(mean(StatSig * ClinSig))
  
}
#renv::snapshot()
