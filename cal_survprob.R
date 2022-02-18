#' Calculate the survival probability from a fitted riAFT-BART model
#'
#' This function calculates the individual survival probability from a fitted riAFT-BART model at desired values of times
#'
#' @param object A fitted object from riAFTBART_estimate() function.
#' @param time.points A numeric vector representing the points at which the survival probability is computed.
#' @param test.only A logical indicating whether or not only data from the test set should be computed. The default is FALSE.
#' @param train.only A logical indicating whether or not only data from the training set should be computed. The default is FALSE.
#' @param cluster.id A vector of integers representing the cluster id. The cluster id should be an integer and start from 1.
#'
#' @return
#' A list with the following two components
#' \item{Surv:}{A matrix of survival probabilities for each individual.}
#' \item{time.points:}{The time point entered.}
#' @export
#'
#' @examples
#' \donttest{
#' library(riAFTBART)
#' set.seed(20181223)
#' n = 50      # number of clusters
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
#' res <- riAFTBART_fit(M.burnin = 50, M.keep = 50, M.thin = 1, status = status,
#'                       y.train = Y, trt.train = trt.train, trt.test = trt.test,
#'                       x.train = cbind(x1,x2),
#'                       x.test = cbind(x1,x2),
#'                       cluster.id = cluster.id)
#'
#' surv_prob_res <- cal_surv_prob(object = res, time.points = sort(exp(logtime)),
#' test.only = TRUE, cluster.id = cluster.id)
#' }
cal_surv_prob <- function(object, time.points, test.only = FALSE, train.only = FALSE, cluster.id){
  ### First check that all elements of time.point > 0
  if(sum(time.points <= 0)) {
    stop("All time points must be positive")
  }
  ### Computation of mean survival curves.
  nsubjects <- nrow(object$tree) # Number of individuals
  nsamples <- ncol(object$tree) # Number of posterior samples
  ngrid <- length(time.points) # Number of time points used to calculate the survival probability
  log.time.points <- log(time.points) # Log time point value
  if(test.only) { # If we only need to calculate the survival probability for the test data
    ntest <- nrow(object$tree.pred) # Number of individuals in the test data
    xind.test <- 1:ntest # Set the id for the ntest individuals in the test data
    SS.test <- matrix(0.0, nrow=ntest, ncol=ngrid) # Set the output matrix with nrow = number of individuals and number of columns = Number of time points used to calculate the survival probability
      for(i in xind.test) {
        for(k in 1:ngrid) {
          Amat <- (log.time.points[k] - object$tree.pred[i,] - object$b[cluster.id[i],])/object$sigma
          SS.test[i,k] <- sum(stats::pnorm(Amat, lower.tail=FALSE))/nsamples # Calculate the survival time using normal approximation
        }
        # print(i)
      }
      SS.test.mean <- colMeans(SS.test) # Calculate the mean survival time across all the individuals
    SS.train <- SS.train.mean <- NULL # Since test.only = T, we set NULL to results with training data
  } else if(train.only) {
      ntrain <- nrow(object$tree) # Number of individuals in the training data
      xind.train <- 1:ntrain # Set the id for the ntrain individuals in the training data
      SS.train <- matrix(0.0, nrow=nsubjects, ncol=ngrid) # Set the output matrix with nrow = number of individuals in the training data and number of columns = Number of time points used to calculate the survival probability
      for(i in xind.train) {
        for(k in 1:ngrid) {
          Amat <- (log.time.points[k] - object$tree[i,] - object$b[cluster.id[i],])/object$sigma
          SS.train[i,k] <- sum(stats::pnorm(Amat, lower.tail=FALSE))/nsamples # Calculate the survival time using normal approximation
        }
        # print(i)
      }
      SS.train.mean <- colMeans(SS.train)  # Calculate the mean survival time across all the individuals
      SS.test <- SS.test.mean <- NULL # Since train.only = T, we set NULL to results with test data
  } else if(!train.only & !test.only) { # Need to calculate the survival probability for both training and test data
      ntest <- ncol(object$m.test)# Number of individuals in the test data
      ntrain <- nrow(object$tree) # Number of individuals in the training data
      xind.train <- 1:ntrain # Set the id for the ntrain individuals in the training data
      xind.test <- 1:ntest # Set the id for the ntest individuals in the test data
      SS.test <- matrix(0.0, nrow=ntest, ncol=ngrid)  # Set the output matrix with nrow = number of individuals in the test data and number of columns = Number of time points used to calculate the survival probability
      SS.train <- matrix(0.0, nrow=ntrain, ncol=ngrid) # Set the output matrix with nrow = number of individuals in the training data and number of columns = Number of time points used to calculate the survival probability
      for(i in xind.test) {
        for(k in 1:ngrid) {
          Amat <- (log.time.points[k] - object$tree.pred[i,] - object$b[cluster.id[i],])/object$sigma
          SS.test[i,k] <- sum(stats::pnorm(Amat, lower.tail=FALSE))/nsamples# Calculate the survival time using normal approximation in the test data
        }
      }
      for(i in xind.train) {
        for(k in 1:ngrid) {
          Amat <- (log.time.points[k] - object$tree[i,] - object$b[cluster.id[i],])/object$sigma
          SS.train[i,k] <- sum(stats::pnorm(Amat, lower.tail=FALSE))/nsamples# Calculate the survival time using normal approximation in the training data
        }
      }
      SS.test.mean <- colMeans(SS.test) # Calculate the mean survival time across all the individuals in the test data
      SS.train.mean <- colMeans(SS.train) # Calculate the mean survival time across all the individuals in the training data
  }
  ans <- list() # List to store the outputs
  class(ans) <- "riAFTBART_survProb"
  ans$Surv.train <- SS.train # Store the individual survival time for the training data
  ans$Surv.test <- SS.test # Store the individual survival time for the test data
  ans$Surv.train.mean <- SS.train.mean # Store the mean individual survival time for the training data
  ans$Surv.test.mean <- SS.test.mean# Store the mean individual survival time for the test data
  ans$time.points <- time.points # Store the time points used for the analysis
  return(ans)
}
