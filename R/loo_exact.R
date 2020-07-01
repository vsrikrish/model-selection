#####################################################################
# loo-exact.R                                                       #
#   This file generates data from a known GEV distribution. It then #
#   computes the exact leave-one-out cross-validation score.        #
#####################################################################

## load packages
library(extRemes) # GEV random number generation and density
library(DEoptim) # identify the MLE using differential evolution
library(adaptMCMC) # run the MCMC chains
library(doParallel) # set up cluster
library(foreach) # parallelized loop over the data
library(plyr) # useful for converting data to different structures

## load scripts
source('R/gev-functions.R')
source('R/gev-priors.R')
source('R/model-selection.R')

## set random number generating seed for reproducability
set.seed(1000)

## set up experiment
dat_len <- 100

# set true model parameters
# get job array ID and look up if the true model is stationary or nonstationary
aid <- as.numeric(Sys.getenv('PBS_ARRAYID'))
true_mod <- c('st', 'nsloc')

loc <- 800 + 0.6*(1:dat_len)*(true_mod == 'nsloc') # location
scale <- 25 # scale
shape <- 0.6 # shape

# generate pseudodata
dat <- extRemes::revd(n=dat_len, loc=loc, scale=scale, shape=shape, type='GEV')

## set vector of parameter names
parnames_all <- c('loc', 'loc1', 'loc2', 'scale', 'scale1', 'scale2', 'shape')

models <- c('st', 'nsloc', 'nslocscale')

# set parameter name subset for each model
parnames <- list()
parnames[['st']] <- c('loc', 'scale', 'shape')
parnames[['nsloc']] <- c('loc1', 'loc2', 'scale', 'shape')
parnames[['nslocscale']] <- c('loc1', 'loc2', 'scale1', 'scale2', 'shape')

# set upper and lower bounds for all parameters
lower_bd <- c('loc'=0, 'loc1'=0, 'loc2'=-50, 'scale'=0, 'scale1'=0, 'scale2'=-20, 'shape'=-2)
upper_bd <- c('loc'=2000, 'loc1'=2000, 'loc2'=50, 'scale'=75, 'scale1'=75, 'scale2'=20, 'shape'=2)

# run DEoptim to identify the MLE for each model
mle <- list()
for (m in models) {
  mle[[m]] <- DEoptim::DEoptim(neg_log_lik_gev, parnames=parnames[[m]], dat=dat, type=m,
                               lower=lower_bd[parnames[[m]]], upper=upper_bd[parnames[[m]]],
                               control=DEoptim.control(itermax=2000, NP=25*length(parnames[[m]]), trace=FALSE))
}

# get priors for each model
priors <- list()
for (m in models) {
  priors[[m]] <- set_priors(parnames[[m]])
}

# set up cluster for fitting each LOO model
cl <- makeCluster(4)
registerDoParallel(cl)
loo_out <- foreach(i=1:4, .packages=c('truncnorm', 'extRemes')) %dopar% {
  # run MCMC for the models without the ith data point, which is held out
  mcmc_out <- list()
  samp <- list()
  dens <- list()
  for (m in c('st', 'nsloc', 'nslocscale')) {
    acc_rate <- 0.234 + ((0.44 - 0.234) / length(parnames[[m]]))
    mcmc_out[[m]] <- adaptMCMC::MCMC(log_post_gev, n=5e5, init=mle[[m]]$optim$bestmem, adapt=TRUE, acc.rate=acc_rate, 
                  gamma=0.67, list=TRUE, parnames=parnames[[m]], dat=dat[-i], type=m, priors=priors[[m]])
    samp[[m]] <- mcmc_out[[m]]$samples[sample(1:nrow(mcmc_out[[m]]$samples[(1e5+1):5e5]), 1e5, replace=TRUE),] # sub-sample MCMC output
    # compute the log-likelihood for the held-out data point and for each MCMC sample
    dens[[m]] <- apply(samp[[m]], 1, log_lik_gev, type=m, parnames=colnames(samp[[m]]), dat=dat[i])
  }
  
  unlist(lapply(dens, function(l) log(mean(exp(l)))))
}
stopCluster(cl)

saveRDS(plyr::ldply(loo_out), paste0('output/loo_exact-', true_mod, '.rds'))
