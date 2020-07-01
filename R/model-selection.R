#####################################################################
# model-selection.R                                                 #
#   This file contains functions computing several model            #
#   selection metrics. Included are AIC, BIC, WAIC, Importance      #
#   Sampled Bayesian Leave-One-Out Cross Validation, and            #
#   Pareto-Smoothed Importance Sampled Bayesian Leave-One-Out       #
#   Cross Validation.                                               #
#####################################################################

#####################################################################
# aic: function to evaluate the Akaike Information Criterion (AIC)  #
#      from Akaike (1974). This uses the maximum-likelihood         #
#      estimate only.                                               #
#                                                                   #
#     Inputs: 1) mle_ll: log-likelihood for the MLE.                #
#             2) np: number of model parameters.                    #
#     Outputs: 1) AIC value.                                        #
#####################################################################
aic <- function(mle_ll, np) { 2 * (np - mle_ll) }

#####################################################################
# bic: function to evaluate the Bayesian Information Criterion      #
#      (BIC) from Schwarz (1978). This uses the maximum-likelihood  #
#      estimate only.                                               #
#                                                                   #
#     Inputs: 1) mle_ll: log-likelihood for the MLE.                #
#             2) np: number of model parameters.                    #
#             3) dat_len: length of the data used to fit the model. #
#     Outputs: 1) BIC value.                                        #
#####################################################################
bic <- function(mle_ll, np, dat_len) { np * log(dat_len) - 2 * mle_ll }

#####################################################################
# waic: function to evaluate the Watanabe-Akaike Information        #
#       Criterion (WAIC) from Watanabe (2010), also discussed in    #
#       Vehtari et al (2017). This is a Bayesian metric which has   #
#       other theoretical advantages over AIC.                      #
#                                                                   #
#     Inputs: 1) mcmc_out: MCMC samples from the posterior. This    #
#                should be non-thinned, but burned-in if desired.   #
#             2) dat: data used to fit the model. Required to       #
#                compute the likelihood for each MCMC sample.       #
#             3) lik_fun: string with name of likelihood function.  #
#     Outputs: 1) WAIC value.                                       #
#####################################################################
waic <- function(mcmc_out, dat, lik_fun) {
  samp <- mcmc_out[sample(1:nrow(mcmc_out), 1e5, replace=TRUE),] # sub-sample MCMC output
  # determine the model type based on the number of parameters
  if (ncol(samp) == 3) {
    type <- 'st'
  } else if (ncol(samp) == 4) {
    type <- 'nsloc'
  } else if (ncol(samp) == 5) {
    type <- 'nslocscale'
  }
  # compute the log-likelihood for each data point and for each MCMC sample
  dens <- apply(samp, 1, match.fun(lik_fun), type=type, parnames=colnames(samp), dat=dat)
  lpd <- sum(log(rowMeans(exp(dens)))) # compute the log-predictive density for each point
  # bias-correct the log-predictive density with an estimate of effective parameter size and return
  lpd - sum(apply(dens, 1, var)) 
}

#####################################################################
# imp_loo: function to evaluate the Bayesian leave-one-out          #
#          cross-validation (LOO-CV) score using raw importance     #
#          sampling to avoid re-fitting the model for each held-out #
#          point (Gelfand 1992, Vehtari et al (2017)).              #
#                                                                   #
#     Inputs: 1) mcmc_out: MCMC samples from the posterior. This    #
#                should be non-thinned, but burned-in if desired.   #
#             2) dat: data used to fit the model. Required to       #
#                compute the likelihood for each MCMC sample.       #
#             3) lik_fun: string with name of likelihood function.  #
#     Outputs: 1) LOO-CV estimate.                                  #
#####################################################################
imp_loo <- function(mcmc_out, dat, lik_fun) {
  samp <- mcmc_out[sample(1:nrow(mcmc_out), 1e5, replace=TRUE),] #sub-sample MCMC chain
  # determine model type using number of parameters
  if (ncol(samp) == 3) {
    type <- 'st'
  } else if (ncol(samp) == 4) {
    type <- 'nsloc'
  } else if (ncol(samp) == 5) {
    type <- 'nslocscale'
  }
  # compute the log-likelihood for each data point and each MCMC sample
  dens <- apply(samp, 1, match.fun(lik_fun), type=type, parnames=colnames(samp), dat=dat)
  imp_wght <- 1 / exp(dens)   # compute the importance weights: the reciprocal of the likelihood
  sum(log(1 / rowMeans(imp_wght))) # return the LOO estimate across all of the data
}

#####################################################################
# ps_loo: function to evaluate the Bayesian leave-one-out           #
#          cross-validation (LOO-CV) score using Pareto-smoothed    #
#          importance sampling. This fits a GPD distribution to the #
#          tails of the importance weights from IS-LOO to smooth    #
#          the extreme weights and stabilize the variance of their  # 
#          distribution.                                            #
#     Inputs: 1) mcmc_out: MCMC samples from the posterior. This    #
#                should be non-thinned, but burned-in if desired.   #
#             2) dat: data used to fit the model. Required to       #
#                compute the likelihood for each MCMC sample.       #
#             3) lik_fun: string with name of likelihood function.  #
#     Outputs: 1) LOO-CV estimate.                                  #
#####################################################################
ps_loo <- function(mcmc_out, dat, lik_fun) {
  samp <- mcmc_out[sample(1:nrow(mcmc_out), 1e5, replace=TRUE),]  # sub-sample MCMC chain
  # determine model type based on number of samples
  if (ncol(samp) == 3) {
    type <- 'st'
  } else if (ncol(samp) == 4) {
    type <- 'nsloc'
  } else if (ncol(samp) == 5) {
    type <- 'nslocscale'
  }
  # compute the log-likelihood for each data point and each MCMC sample
  dens <- apply(samp, 1, log_lik_gev, type=type, parnames=colnames(samp), dat=dat)
  imp_wght <- 1 / exp(dens) # compute the importance weights: the recriprocal of the likelihood
  # we use the top 20% of the importance weights as the extremes to be smoothed
  nwght <- ncol(imp_wght) / 5 # find the number of weights to be smoothed
  gpd_shape <- numeric(nrow(imp_wght)) # store the shape distributions of the GPDs for diagnostics
  smooth_wght <- matrix(0, ncol=ncol(imp_wght), nrow=nrow(imp_wght)) # set storage for smoothed weights
  # smooth the weights using a fitted generalized Pareto distribution
  for (i in 1:nrow(imp_wght)) {
    # sort the weights to identify the top 20%
    wght_srt <- sort(imp_wght[i, ], decreasing=TRUE)
    thresh <- wght_srt[nwght] # set the threshold for the GPD
    gpd_fit <- extRemes::fevd(imp_wght[i, ], threshold=thresh, type='GP')$results$par # fit the GPD
    gpd_shape[i] <- gpd_fit['shape'] 
    # find the expected order statistics of the fitted GPDs which will replace the extreme weights
    repl_wght <- extRemes::qevd(((1:nwght) - 0.5) / nwght, scale=gpd_fit['scale'], 
                                shape=gpd_fit['shape'], type='GP', threshold=thresh)
    wght_srt[1:nwght] <- repl_wght
    # truncate the weights to ensure a finite variance
    truncval <- length(wght_srt)^(3/4) * mean(wght_srt)
    wght_srt[wght_srt > truncval] <- truncval
    smooth_wght[i, ] <- wght_srt # store the smoothed weights
  }
  # return the PS-LOO estimate
  sum(log(rowSums(smooth_wght * exp(dens)) / rowSums(smooth_wght)))
}
