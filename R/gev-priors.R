#####################################################################
# gev-functions.R                                                   #
#   This file contains functions specifying the likelihood          #
#   functions for several GEV models as well as a function to       #
#   evaluate priors and posteriors.                                 #
#####################################################################

## specify prior distributions
# we do this for each parameter that might be included in any given model
set_priors <- function(parnames) {
  priors <- vector('list', length(parnames))
  names(priors) <- parnames
  for (p in parnames) {
    if (p %in% c('loc', 'loc1')) {
      priors[[p]] <- vector('list', 3)
      names(priors[[p]]) <- c('dens_fun', 'mean', 'sd')
      # stationary location parameter or nonstationary location intercept
      priors[[p]][['dens_fun']] <- 'dnorm'
      priors[[p]][['mean']] <- 1000
      priors[[p]][['sd']] <- 300
    } else if (p == 'loc2') {
    # nonstationary location slope
      priors[[p]] <- vector('list', 3)
      names(priors[[p]]) <- c('dens_fun', 'mean', 'sd')
      priors[[p]][['dens_fun']] <- 'dnorm'
      priors[[p]][['mean']] <- 0
      priors[[p]][['sd']] <- 50
    } else if (p %in% c('scale', 'scale1')) {
      # stationary scale parameter or nonstationary scale intercept
      priors[[p]] <- vector('list', 4)
      names(priors[[p]]) <- c('dens_fun', 'a', 'mean', 'sd')
      priors[[p]][['dens_fun']] <- 'dtruncnorm'
      priors[[p]][['a']] <- 0
      priors[[p]][['mean']] <- 10
      priors[[p]][['sd']] <- 20
    } else if (p == 'scale2') {
      priors[[p]] <- vector('list', 3)
      names(priors[[p]]) <- c('dens_fun', 'mean', 'sd')
      # nonstationary scale slope
      priors[[p]][['dens_fun']] <- 'dnorm'
      priors[[p]][['mean']] <- 0
      priors[[p]][['sd']] <- 10
    } else if (p == 'shape') {
      priors[[p]] <- vector('list', 3)
      names(priors[[p]]) <- c('dens_fun', 'mean', 'sd')
      # shape parameter
      priors[[p]][['dens_fun']] <- 'dnorm'
      priors[[p]][['mean']] <- 0
      priors[[p]][['sd']] <- 0.5
    }
  }
  priors
}