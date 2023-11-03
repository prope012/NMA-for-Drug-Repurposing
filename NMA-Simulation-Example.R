###################################################################
#### An example showing how to use the functions in  ##############
################### NMA-Simulation-Functions.R ####################
###################################################################

#renv::init()

#-----------------------------------------------#
#------ Source Functions & Load Packages -------#
#-----------------------------------------------#

### source in the NMA-Simulation-Functions.R script
work.dir <- "placeholder"
source(paste(work.dir,"NMA-Simulation-Functions.R",sep=""))
library(tidyr)
library(dplyr)
library(locfit)

### Establish matrices containing the treatments studied in each arm for each indication.
ind1.trt <- matrix(c(rep(1,28),rep(2,3),rep(3,3),rep(5,3),rep(2,2),7,
                     rep(2,8),rep(3,4),rep(4,4),rep(5,2),rep(6,5),rep(7,3),
                     rep(9,2),rep(8,3),rep(4,3),rep(6,3),rep(5,2),8), ncol=2)

ind2.trt <- matrix(c(rep(1,27),rep(3,2),rep(4,3),rep(5,3),rep(2,3),rep(3,2),
                     rep(4,2),rep(5,6),rep(7,3),rep(8,5),rep(9,6),rep(6,2),
                     rep(8,3),rep(7,3)), ncol = 2)

### Find random variances such that sg^2+se^2 = 0.50 and sg^2/(sg^2+se^2) = 0.1, 0.3, 0.5, 0.7, and 0.9
#	Solving this system of equations, you get: sg2 = 0.5*rho and se2 = 0.5*(1-rho)
get.var <- function(rho){
  varg <- 0.5*rho; vare <- 0.5*(1-rho)
  return(c(varg,vare))}	
corr <- seq(0.1,0.9,0.2)
var.g <- sapply(corr, function(x) get.var(x)[1])
var.e <- sapply(corr, function(x) get.var(x)[2])

### Number of randomly generated sets of d's
nsims <- 2

### Define number of iterations (number of data sets within a set of d's)
niters <- 10


#---------------------------------------------------------------------#
#--------------- Generate data using Cor(d^1_1k,d^2_1k) = 0.5 --------#
#------ and excluding treatment 2 from the indication 2 network.------#
#------ Then, fit the mixed effects model to each data set. ----------#
#---------------------------------------------------------------------#

Mixed.Mod.Trt2.Corr5 <- lapply(1:nsims, function(y){
  
  # set seed
  set.seed(y)
  
  # Randomly generate basic parameters and baseline effect parameters; get population-averaged LORs for both indications
  foo <- sim.pars(beta0 = 3.4, beta1 = 1.3, beta2 = -1.6, sigma.g = sqrt(var.g[3]), sigma.e = sqrt(var.e[3]), class.ind = c(1,2,2,2,1,2,1,2))
  
  # get data sets that exclude each treatment from indication 2 one at a time
  LOOCV.dat <- get.LOOCV.dat(J = 75, 
                             num.arms = rep(2,75),
                             indication = c(rep(1,40), rep(2,35)),
                             trt = rbind(ind1.trt, ind2.trt),
                             n = matrix(c(rep(60,80),rep(38,70)), nrow = 75, ncol = 2, byrow = T),
                             d = rbind(foo$d1_vals, foo$d2_vals),
                             mu <- c(foo$mu1,foo$mu2),
                             sigma = c(0.25,0.25),
                             niters = niters)
  
  # Fit model to the data sets that exclude treatment 2 from indication 2
  res <- lapply(LOOCV.dat[["Minus Trt 2"]], function(x){
    
    # Fit model
    mod <- fit.2ind.mod(dat = x, mod.type = "mixed", class = c(NA,1,2,2,2,1,2,1,2), sd.prior = "Ht")
    
    # Return fitted model results
    return(list("Ht" = mod))
  })
  
  # count
  print(y)
  
  return(res)
  
})


#-----------------------------------------------#
#---------- Compute summary statistics ---------#
#-----------------------------------------------#

### Get true values for the basic parameters
true.d.vals.cor5 <- lapply(1:nsims, function(x){
  set.seed(x)
  sim.pars(beta0 = 3.4, beta1 = 1.3, beta2 = -1.6, sigma.g = sqrt(var.g[3]), sigma.e = sqrt(var.e[3]), class.ind = c(1,2,2,2,1,2,1,2))
})

### Get summary stats
Mixed.Mod.Trt2.Corr5.Ht <- lapply(Mixed.Mod.Trt2.Corr5, function(x){foo <- lapply(x, function(y) y[["Ht"]]); return(foo)})
Mixed.Mod.Trt2.Corr5.Ht.Stats <- list()
for(i in 1:nsims){
  Mixed.Mod.Trt2.Corr5.Ht.Stats[[i]] <- get.summary.stats(sim.list = Mixed.Mod.Trt2.Corr5.Ht[[i]], true.vals = c(true.d.vals.cor5[[i]]$d[,1],true.d.vals.cor5[[i]]$d[,2]),
                                                nt = 9, inds = "Both", mod.type = "mixed")
}

#-------------------------------------------------------#
#----- Compute PoS for Treatment 2 in Indication 2 -----#
#-------------------------------------------------------#

### PoS should be computed using a real data set. For this example, a data set will be randomly generated.
set.seed(28)
foo <- sim.pars(beta0 = 3.4, beta1 = 1.3, beta2 = -1.6, sigma.g = sqrt(var.g[3]), sigma.e = sqrt(var.e[3]), class.ind = c(1,2,2,2,1,2,1,2))
pos.dat <- get.LOOCV.dat(J = 75, 
                       num.arms = rep(2,75),
                       indication = c(rep(1,40), rep(2,35)),
                       trt = rbind(ind1.trt, ind2.trt),
                       n = matrix(c(rep(60,80),rep(38,70)), nrow = 75, ncol = 2, byrow = T),
                       d = rbind(foo$d1_vals, foo$d2_vals),
                       mu <- c(foo$mu1,foo$mu2),
                       sigma = c(0.25,0.25),
                       niters = 1)$'Minus Trt 2'

### Prepare data for NMA
pos.dat <- pos.dat[[1]]
ns <- length(unique(pos.dat$study)) # number of studies
na <- (pos.dat %>% group_by(study) %>% dplyr::summarize(na = length(study)))$na # number of arms per study
r <- NULL; for(i in 1:max(na)){r <- cbind(r, (pos.dat %>% group_by(study) %>% summarise(r=r[i]))$r)} # wide format for number of successes
n <- NULL; for(i in 1:max(na)){n <- cbind(n, (pos.dat %>% group_by(study) %>% summarise(n=n[i]))$n)} # wide format for number participants
t <- NULL; for(i in 1:max(na)){t <- cbind(t, (pos.dat %>% group_by(study) %>% summarise(t=trt[i]))$t)} # wide format for treatments
nt <- max(pos.dat$trt) # number of treatments in the network 
indication <- (pos.dat %>% group_by(study) %>% dplyr::summarize(ind = indication[1]))$ind # indication indicator

### Fit desired model
jags.script <- paste("Two Indication Models/","mixed","-random-effects-","Ht",".txt",sep="")
params = c("d","beta0","beta2","gamma","sdg","sde","sdd","d.new")
data <- list(ns = ns, r = r, t = t, n = n, nt=nt, indication=indication, na = na)
jags.fit =jags(data = data, n.burnin = 5000, n.iter = 15000, n.thin = 1, jags.module = "glm",
               n.chains = 2, parameters.to.save = params, model.file = paste(jags.script.dir,jags.script,sep=""),
               quiet = TRUE, progress.bar = "none")

### Compute PoS, where success is defined as achieving statistical and clinical significance.
get.pos(fitted.jags.mod = jags.fit, trt.num = 2, mvn = FALSE, nt = 9, piC = 0.09, n = 60, delta = 2.5, niters=1000)

#renv::snapshot()
