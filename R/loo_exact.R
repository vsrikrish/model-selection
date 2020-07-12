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
dat_len <- 100 # pseudodata length

# MCMC and posterior analysis parameters
burnin <- 1e5 # length of burnin
niter <- 5e5 # number of MCMC iterations
nsamp <- 1e5 # number of posterior samples for LOO estimation

models <- c('st', 'nsloc', 'nslocscale')

# generate pseudodata
trend <- (1:(dat_len*365))*1.5/365
per_semi <- sin(2*pi*(1:(dat_len*365))/180 + 90)*80
per_ann <- sin(2*pi*(1:(dat_len*365))/360 - 45)*300
per_all <- per_semi + per_ann + trend
# generate white noise from AR process
wn <- numeric(dat_len*365)
set.seed(1000)
wn[1] <- 25
for (t in 2:(100*365)) {
  wn[t] = wn[t-1]*0.6 + rnorm(1, mean=0, sd=150)
}
dat_all <- per_all + wn
# find annual maxima from the pseudodata
dat_m <- matrix(dat_all, nrow=365)
dat <- apply(dat_m, 2, max)

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

# get priors for each model
priors <- list()
for (m in models) {
  priors[[m]] <- set_priors(parnames[[m]])
}

# set up cluster for fitting each LOO model
cl <- makeCluster(detectCores())
registerDoParallel(cl)
# loop over each time point and fit models for each time point within the loop
loo_out <- foreach(i=1:dat_len, .packages=c('truncnorm', 'extRemes', 'DEoptim')) %dopar% {
  # initialize lists for storage across models for this time
  dens <- list()
  for (m in c('st', 'nsloc', 'nslocscale')) {
    mle <- DEoptim::DEoptim(neg_log_lik_gev, parnames=parnames[[m]], dat=dat[-i], type=m,
                                 lower=lower_bd[parnames[[m]]], upper=upper_bd[parnames[[m]]],
                                 control=DEoptim.control(itermax=2000, NP=25*length(parnames[[m]]), trace=FALSE))

    acc_rate <- 0.234 + ((0.44 - 0.234) / length(parnames[[m]]))
    mcmc_out <- adaptMCMC::MCMC(log_post_gev, n=niter, init=mle$optim$bestmem, adapt=TRUE, acc.rate=acc_rate, 
                  gamma=0.67, list=TRUE, parnames=parnames[[m]], dat=dat[-i], type=m, priors=priors[[m]])
    idx <- sample((burnin+1):niter, nsamp, replace=TRUE) # sample post-burnin indices
    samp <- mcmc_out$samples[idx,] # sub-sample MCMC output
    # compute the log-likelihood for the held-out data point and for each MCMC sample
    dens[[m]] <- apply(samp, 1, log_lik_gev, type=m, parnames=colnames(samp), dat=dat[i])
  }
  
  unlist(lapply(dens, function(l) log(mean(exp(l)))))
}
stopCluster(cl)

saveRDS(plyr::ldply(loo_out), paste0('output/loo_exact.rds'))
