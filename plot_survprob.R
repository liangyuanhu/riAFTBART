#' Plot the fitted survival curves from riAFT-BART model
#'
#' This function plot the mean/individual survival curves from a fitted riAFT-BART model
#'
#' @param x An object from cal_surv_prob() function.
#' @param test.only A logical indicating whether or not only data from the test set should be computed. The default is FALSE.
#' @param train.only A logical indicating whether or not only data from the training set should be computed. The default is FALSE.
#' @param id A vector representing the IDs for the individual survival curves to plot. The default is NULL and the mean survival curves will be plotted.
#' @param ... further arguments passed to or from other methods.
#'
#' @return A plot
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
#' surv_prob_res <- cal_surv_prob(object = res, time.points = sort(exp(logtime)),
#' test.only = TRUE, cluster.id = cluster.id)
#' plot(x = surv_prob_res, test.only = TRUE, train.only = FALSE)
#' }
plot.riAFTBART_survProb <- function(x, test.only = FALSE, train.only = TRUE, id = NULL,...){
  if (test.only) { # If we only need to plot the survival probability for the test data
    if (is.null(id)) { # If we only need to plot the mean survival probability for the test data
      plot(x$time.points, x$Surv.test.mean, type = "l", xlab  = "Time points", ylab = "Predicted survival probability")
    } else { # Plot the individual survival curves for the specific individual
      plot(x$time.points, x$Surv.test[id,], type = "l", xlab  = "Time points", ylab = "Predicted survival probability")
    }
  } else { # If we need to plot the survival probability for the training data
    if (is.null(id)) { # If we only need to plot the mean survival probability
      plot(x$time.points, x$Surv.train.mean, type = "l", xlab  = "Time points", ylab = "Predicted survival probability")
    } else { # Plot the individual survival curves for the specific individual
      plot(x$time.points, x$Surv.train[id,], type = "l", xlab  = "Time points", ylab = "Predicted survival probability")
    }
  }
}
