#' Fit a random effect accelerated failure time BART model
#'
#' This function implements the random effect accelerated failure time BART (riAFT-BART) algorithm.
#'
#' @param M.burnin A numeric value indicating the number of MCMC iterations to be treated as burn in.
#' @param M.keep A numeric value indicating the number of MCMC posterior draws after burn in.
#' @param M.thin A numeric value indicating the thinning parameter.
#' @param status A vector of event indicators: status = 1 indicates that the event was observed while status = 0 indicates the observation was right-censored.
#' @param y.train A vector of follow-up times.
#' @param x.train A dataframe or matrix, including all the covariates but not treatments for training data, with rows corresponding to observations and columns to variables.
#' @param trt.train A numeric vector representing the treatment groups for the training data.
#' @param x.test A dataframe or matrix, including all the covariates but not treatments for testing data, with  rows corresponding to observations and columns to variables.
#' @param SA A logical indicating whether to conduct sensitivity analysis. The default is FALSE.
#' @param prior_c_function_used Prior confounding functions used for SA, which is inherited from the sa function. The default is NULL.
#' @param gps Generalized propensity score, which is inherited from the sa function. The default is NULL.
#' @param trt.test A numeric vector representing the treatment groups for the testing data.
#' @param cluster.id A vector of integers representing the clustering id. The cluster id should be an integer and start from 1.
#' @param verbose A logical indicating whether to show the progress bar. The default is FALSE
#' @return A list with the following elements:
#' \item{b:}{A matrix including samples from the posterior of the random effects.}
#' \item{tree:}{A matrix with M.keep rows and nrow(x.train) columns represnting the predicted log survival time for x.train.}
#' \item{tree.pred:}{A matrix with M.keep rows and nrow(x.test) columns represnting the predicted log survival time for x.test.}
#' \item{tau:}{A vector representing the posterior samples of tau, the standard deviation of the random effects.}
#' \item{sigma:}{A vector representing the posterior samples of sigma, the residual/error standard deviation.}
#' \item{vip:}{A matrix with M.keep rows and ncol(x.train) columns represnting the variable inclusion proportions for each variable.}
#' @export
#'
#' @examples
#' \donttest{
#' library(riAFTBART)
#' set.seed(20181223)
#' n = 5       # number of clusters
#' k = 50      # cluster size
#' N = n*k     # total sample size
#' cluster.id = rep(1:n, each=k)

#' tau.error = 0.8
#' b = stats::rnorm(n, 0, tau.error)

#' alpha = 2
#' beta1 = 1
#' beta2 = -1
#' sig.error = 0.5
#' censoring.rate = 0.02

#' x1 = stats::rnorm(N,0.5,1)
#' x2 = stats::rnorm(N,1.5,0.5)
#' trt.train = sample(c(1,2,3), N, prob = c(0.4,0.3,0.2), replace = TRUE)
#' trt.test = sample(c(1,2,3), N, prob = c(0.3,0.4,0.2), replace = TRUE)
#' error = stats::rnorm(N,0,sig.error)

#' logtime = alpha + beta1*x1 + beta2*x2 + b[cluster.id] + error
#' y = exp(logtime)
#' C = rexp(N, rate=censoring.rate) # censoring times
#' Y = pmin(y,C)
#' status = as.numeric(y<=C)
#' res <- riAFTBART_fit(M.burnin = 10, M.keep = 10, M.thin = 1, status = status,
#'                       y.train = Y, trt.train = trt.train, trt.test = trt.test,
#'                       x.train = cbind(x1,x2),
#'                       x.test = cbind(x1,x2),
#'                       cluster.id = cluster.id)
#'}

riAFTBART_fit <- function(M.burnin, M.keep, M.thin = 1, status, y.train, x.train, trt.train, x.test, trt.test, cluster.id, verbose = FALSE, SA = FALSE,prior_c_function_used = NULL, gps = NULL)
{
  # initial values
  N = length(y.train)  # total sample size
  n = length(unique(cluster.id)) # number of clusters
  null_aft_model <- survival::survreg(survival::Surv(y.train, status) ~ 1, dist="lognormal") # Follow Henderson 2020, we first fit a parametric AFT model that only has an intercept in the model and that assumes log-normal residual to get the priors for the intercept and scale
  null_intercept <- null_aft_model$coefficients # Estimate of intercept from parametric AFT model
  aft_model_scale <- null_aft_model$scale # Estimate of the scale from parametric AFT model
  aft_model <- survival::survreg(survival::Surv(y.train, status) ~ ., dist="lognormal", data = as.data.frame(cbind(x.train, y.train, status))) # Fit a parametric AFT model that uses all of the predictors and that assumes log-normal residual to get the priors for random intercepts and sigmas
  b.init = tapply(stats::resid(aft_model), cluster.id, mean) # Initial random intercepts: mean of the lognormal residuals from the fitted AFT model
  sigma.init <- stats::sd(tapply(stats::resid(aft_model), cluster.id, mean)) # Initial sigma: Standard deviation of the model residuals over K cluster
  x.train <- cbind(x.train, trt.train) # add the treatment indicator to the x.train matrix
  ncol.x.train <- ncol(x.train) # Number of columns for the new x.train matrix
  colnames(x.train) <- NULL # Remove the column names for x.train
  x.test <- cbind(x.test, trt.test) # add the treatment indicator to the x.test matrix
  colnames(x.test) <- NULL # Remove the column names for x.test
  tau.init = 1 # Set the initial value of tau to 1
  alpha.init = 1 # Set the initial value of alpha to 1
  tree.init = stats::rnorm(N) # Set the initial values from the BART trees from rnorm(0,1)
  logT.init = log(y.train)-null_intercept # Follow Henderson 2020, use the transformed "centered" outcome for the MCMC iterations
  # Define and initialize the list of things to keep track of in the "current state" of the chain
  cur <- list(tree=tree.init, b=b.init, sigma=sigma.init, tau=tau.init, alpha = alpha.init, logT=logT.init)

  # Define matrix to store MCMC results
  P <- ncol(x.train)     # number of covariates
  N.pred <- nrow(x.test) # sample size of test data

  chain.keep <- list() # List to store all the outcomes
  chain.keep$vip <- array(NA, dim = c(ncol(dbarts::makeModelMatrixFromDataFrame(as.data.frame(x.train))), M.keep)) # Matrix to store VIP
  chain.keep$b <- array(NA, dim = c(n, M.keep))       # Matrix to store random effects
  chain.keep$tree <- array(NA, dim = c(N, M.keep))    # Matrix to store predicted tree values from training data
  chain.keep$tree.pred <- array(NA, dim = c(N.pred, M.keep))    # Matrix to store predicted tree values from test data
  chain.keep$tau <- array(NA, dim = c(1,M.keep))      # Matrix to store standard deviation of random intercept
  chain.keep$sigma <- array(NA, dim = c(1,M.keep))    # Matrix to store the standard deviation of residual error

  rownames(chain.keep$tau) <- 'tau' # Add the row names for tau
  rownames(chain.keep$sigma) <- 'sigma' # Add the row names for sigma
  rownames(chain.keep$b) <- paste('b', 1:n, sep="") # Add the row names for random intercept
  rownames(chain.keep$tree) <- paste('tree', 1:N, sep="") # Add the row names for predicted tree values from training data
  rownames(chain.keep$tree.pred) <- paste('tree.pred', 1:N.pred, sep="") # Add the row names for predicted tree values from test data
  if (verbose == TRUE){
  cat("Running the MCMC algorithm. It may be long, keep cool :)", "\n") # Words to be printed when verbose == TRUE
  cat("\n")
}
  # Burn-in phase: do not keep these results
  if (M.burnin > 0) {
    if (verbose == TRUE){
    cat("Burn-in phase", "\n") # Words to be printed when verbose == TRUE
    }
    for (i in 1:M.burnin) {
      cur <- blocked.mcmc.update(status, y.train, x.train, x.test, cluster.id, cur,SA = SA,trt.train,prior_c_function_used,gps) # Use the blocked.mcmc.update function to update the results for the burn in stage
      if (verbose == TRUE){
        progressBar(i, M.burnin) # progressBar() function to show the progress
      }

    } # end of i for-loop
  } # end of if-loop

  # Converged phase: keep these results
  if (verbose == TRUE){
  cat("\n")
  cat("Converged phase", "\n") # Words to be printed when verbose == TRUE
}
  for (m in 1:M.keep) {
    if (M.thin > 1) {	# Skip the "thinned" pieces of the chain, if M.thin > 1
      for (j in 1:(M.thin-1)) {
        cur <- blocked.mcmc.update(status, y.train, x.train, x.test, cluster.id, cur,SA = SA,trt.train,prior_c_function_used,gps) # Use the blocked.mcmc.update function to update the results
      }
    }
    cur <- blocked.mcmc.update(status, y.train, x.train, x.test, cluster.id, cur,SA = SA,trt.train,prior_c_function_used,gps) # Use the blocked.mcmc.update function to update the results for the posterior draws
    if (verbose == TRUE){
      progressBar(m, M.keep) # progressBar function to show the progress
    }
    chain.keep$vip[,m] <- cur$vip # Store the VIP
    chain.keep$b[,m] <- cur$b # Store the random intercepts
    chain.keep$tree[,m] <- cur$tree # Store the predicted tree values from training data
    chain.keep$tau[m] <- cur$tau # Store the standard deviation of random intercept
    chain.keep$sigma[m] <- cur$sigma # Store the standard deviation of residual error
    chain.keep$tree.pred[,m] <- cur$tree.pred # Store the predicted tree values from test data

  } # end of m for-loop
  chain.keep$vip <- t(chain.keep$vip) # Final VIP output
  chain.keep$tau <- c(chain.keep$tau) # Final standard deviation of random intercept
  chain.keep$sigma <- c(chain.keep$sigma) # Final standard deviation of residual error
  chain.keep$tree <- chain.keep$tree + null_intercept # Final predicted tree values from training data (in logT scale and we need center back)
  chain.keep$tree.pred <- chain.keep$tree.pred + null_intercept # Final predicted tree values from test data (in logT scale and we need center back)
  class(chain.keep) <- "riAFTBART_estimate" # Add class value for the plot generic function
  return(chain.keep)

}
