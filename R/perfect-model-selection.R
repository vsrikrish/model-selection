#####################################################################
# perfect-model-selection.R                                         #
#   This file generates data from a known GEV distribution. It then #
#   finds the maximum likelihood estimates and does Bayesian        #
#   inversion using Markov Chain Monte Carlo (MCMC) for a few       #
#   different GEV structureas. Various model selection criteria are # 
#   then calculated for each model.                                 #
#####################################################################

## load packages
library(extRemes) # GEV random number generation and density
library(DEoptim) # identify the MLE using differential evolution
library(adaptMCMC) # run the MCMC chains

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
loc <- 800 # location
scale <- 25 # scale
shape <- 0.6 # shape

# generate pseudodata
dat <- extRemes::revd(n=dat_len, loc=loc, scale=scale, shape=shape, type='GEV')

## set vector of parameter names
parnames_all <- c('loc', 'loc1', 'loc2', 'scale', 'scale1', 'scale2', 'shape')

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
mle_all <- list()
for (m in models) {
  mle[[m]] <- DEoptim::DEoptim(neg_log_lik_gev, parnames=parnames[[m]], dat=dat, type=m,
                               lower=lower_bd[parnames[[m]]], upper=upper_bd[parnames[[m]]],
                               control=DEoptim.control(itermax=2000, NP=25*length(parnames[[m]]), trace=FALSE))
  mle_all[[m]] <- mle[[m]]$optim$bestval
}

saveRDS(mle, paste0('output/mle-', true_mod, '.rds'))

# get priors for each model
priors <- list()
for (m in models) {
  priors[[m]] <- set_priors(parnames[[m]])
}

# run MCMC for each model
mcmc_out <- list()
mcmc_all <- list()
for (m in c('st', 'nsloc', 'nslocscale')) {
  acc_rate <- 0.234 + ((0.44 - 0.234) / length(parnames[[m]]))
  mcmc_out[[m]] <- adaptMCMC::MCMC(log_post_gev, n=1e5, init=mle[[m]]$optim$bestmem, adapt=TRUE, acc.rate=acc_rate, 
                                   gamma=0.67, list=TRUE, parnames=parnames[[m]], dat=dat, type=m, priors=priors[[m]])
  mcmc_all[[m]] <- mcmc_out[[m]]$samples[(5e4+1):1e5,] # just assume burnin is half the run for the time being
  
}  

saveRDS(mcmc, paste0('output/mcmc-', true_mod, '.rds'))

np <- list('st'=3, 'nsloc'=4, 'nslocscale'=5) # we need the number of parameters for AIC and BIC

## compute model selection metrics
# Divide AIC and BIC by 2 to make them comparable to LOO-CV, as opposed to on the deviance scale
aic_all <- mapply(aic, mle_ll = mle_all, np = np) / 2
bic_all <- mapply(bic, mle_ll = mle_all, np=np, dat_len=length(dat)) / 2
waic_all <- unlist(lapply(mcmc_all, waic, dat=dat, lik_fun='log_lik_gev'))
isloo_all <- unlist(lapply(mcmc_all, imp_loo, dat=dat, lik_fun='log_lik_gev'))
psloo_all <- unlist(lapply(mcmc_all, ps_loo, dat=dat, lik_fun='log_lik_gev'))
# load exact LOO for comparison
loo_out <- readRDS(paste0('output/loo_exact-', true_mod, '.rds'))
loo_all <- colSums(loo_out)

# combine into a table
metrics <- data.frame('LOO'=loo_all, 'AIC'=aic_all, 'BIC'=bic_all, 'WAIC'=waic_all, 'IS-LOO'=isloo_all, 'PS-LOO'=psloo_all)
print(metrics)