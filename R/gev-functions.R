#####################################################################
# gev-functions.R                                                   #
#   This file contains functions specifying the likelihood          #
#   functions for several GEV models as well as a function to       #
#   evaluate priors and posteriors.                                 #
#####################################################################

#####################################################################
# log_lik_gev: function to evaluate the log-likelihood for          #
#               GEV models. This function just calls                #
#               extRemes::qevd after unpacking the parameters       #
#               based on the model type, then sums across the       # 
#               data.                                               #
#                                                                   #
#     Inputs: 1) params: vector of parameters, which is then        #
#                        unpacked based on the model type.          #
#             2) parnames: names of parameters.                     #
#             3) dat: data stored in a structure which can be       #
#                     coerced to a numeric vector                   #
#             4) type: string corresponding to the type of GEV.     #
#                     Supported types are 'st' (stationary),        #
#                     'nsloc' (non-stationary location),            #
#                     and 'nslocscale' (non-stationary location     #
#                     and scale).                                   #
#     Outputs: 1) log-likelihood value for the given parameters.    #
#####################################################################

library(extRemes) # load package for GEV density function
library(truncnorm) # load package for truncated normal density function

log_lik_gev <- function(params, parnames, dat, type) {
  # compute location, scale, and shape based on passed-in parameters and type
  if (type == 'st') {
    loc <- params[match('loc', parnames)]
    scale <- params[match('scale', parnames)]
    shape <- params[match('shape', parnames)]
  } else if (type == 'nsloc') {
    # non-stationary location is the has the form loc1 + t * loc2
    loc <- params[match('loc1', parnames)] + (1:length(dat)) * params[match('loc2', parnames)]
    scale <- params[match('scale', parnames)]
    shape <- params[match('shape', parnames)]
  } else if (type == 'nslocscale') {
    # non-stationary location is the has the form loc1 + t * loc2
    loc <- params[match('loc1', parnames)] + (1:length(dat)) * params[match('loc2', parnames)]
    # non-stationary scale is the has the form scale1 + t * scale2
    scale <- params[match('scale1', parnames)] + (1:length(dat)) * params[match('scale2', parnames)]
    shape <- params[match('shape', parnames)]
  } else {
    stop('GEV type must be st, nsloc, or nslocscale!')
  }
  # check if the scale parameter is positive and return -Inf if not
  if (sum(scale > 0) < length(scale)) {
    return(-Inf)
  }
  # get the log-likelihood and return
  extRemes::devd(dat, loc=loc, scale=scale, shape=shape, type='GEV', log=TRUE)
}

#####################################################################
# log_lik_gev: function to evaluate the log-likelihood for          #
#               GEV models. This function just calls                #
#               extRemes::qevd after unpacking the parameters       #
#               based on the model type, then sums across the       # 
#               data.                                               #
#                                                                   #
#     Inputs: 1) params: vector of parameters, which is then        #
#                        unpacked based on the model type.          #
#             2) parnames: names of parameters.                     #
#             3) priors: list containing information for each prior #
#                       distribution. Each element should be named  #
#                       based on the appropriate parname entry,     #
#                       and should be a list containing the         #
#                       following information:                      #
#                         'dens_fun': string name of density        #
#                                     function to evaluate;         #
#                          necessary arguments for the density      #
#                             function, named for the argument.     #
#     Outputs: 1) log-prior value for the given parameters.         #
#####################################################################

log_prior <- function(params, parnames, priors) {
  # this function evaluates the log-prior density for a given parameter
  log_dens <- function(name) {
    val <- params[match(name, parnames)]
    if (priors[[name]][['dens_fun']] == 'dtruncnorm') {
      return(log(do.call(match.fun(priors[[name]][['dens_fun']]),
              c(list(x=val),
                priors[[name]][-which(names(priors[[name]]) %in% c('dens_fun'))])
              )))
    } else {
      return(do.call(match.fun(priors[[name]][['dens_fun']]),
              c(list(x=val, log=TRUE),
                priors[[name]][-which(names(priors[[name]]) %in% c('dens_fun'))])
      ))
    }
  }
  # evaluate log-prior densities for each parameter
  lp <- vapply(parnames, log_dens, numeric(1))
  
  # return sum of log-priors
  lp
}

#####################################################################
# log_post_gev: function to evaluate the log-posterior for          #
#               GEV models.                                         #
#                                                                   #
#     Inputs: 1) params: vector of parameters, which is then        #
#                        unpacked based on the model type.          #
#             2) parnames: names of parameters.                     #
#             3) dat: data stored in a structure which can be       #
#                     coerced to a numeric vector                   #
#             4) type: string corresponding to the type of GEV.     #
#                     Supported types are 'st' (stationary),        #
#                     'nsloc' (non-stationary location),            #
#                     and 'nslocscale' (non-stationary location     #
#                     and scale).                                   #
#             5) priors: list containing information for each prior #
#                       distribution. Each element should be named  #
#                       based on the appropriate parname entry,     #
#                       and should be a list containing the         #
#                       following information:                      #
#                         'dens_fun': string name of density        #
#                                     function to evaluate;         #
#                          necessary arguments for the density      #
#                             function, named for the argument.     #
#     Outputs: 1) log-posterior value for the given parameters.     #
#####################################################################

log_post_gev <- function(params, parnames, dat, type, priors) {
  # get log-prior density value
  lp <- sum(log_prior(params, parnames, priors))
  # if log-prior density is -Inf, no need to evaluate the likelihood
  if (lp == -Inf) {
    return(-Inf)
  }
  # evaluate log-likelihood
  ll <- log_lik_gev(params, parnames, dat, type)
  # return log-posterior
  sum(ll + lp)
}

#####################################################################
# neg_log_lik_gev: function to evaluate the negative of the         #
#               log-likelihood for GEV models. This function is     #
#               intended to be used with an optimization routine    #
#               such as DEoptim.                                    #
#                                                                   #
#     Inputs: 1) params: vector of parameters, which is then        #
#                        unpacked based on the model type.          #
#             2) parnames: names of parameters.                     #
#             3) dat: data stored in a structure which can be       #
#                     coerced to a numeric vector                   #
#             4) type: string corresponding to the type of GEV.     #
#                     Supported types are 'st' (stationary),        #
#                     'nsloc' (non-stationary location),            #
#                     and 'nslocscale' (non-stationary location     #
#                     and scale).                                   #
#     Outputs: 1) negative log-likelihood value for the given       #
#                 parameters.                                       #
#####################################################################

neg_log_lik_gev <- function(params, parnames, dat, type) {
  -1 * sum(log_lik_gev(params, parnames, dat, type))
}
