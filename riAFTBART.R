#' A flexible approach for causal inference with multiple treatments and clustered survival outcomes
#'
#' This function implements the random effect accelerated failure time BART (riAFT-BART) for causal inference with multiple treatments and clustered survival outcomes.
#' @param M.burnin A numeric value indicating the number of MCMC iterations to be treated as burn in.
#' @param M.keep A numeric value indicating the number of MCMC posterior draws after burn in.
#' @param M.thin A numeric value indicating the thinning parameter.
#' @param status A vector of event indicators: status = 1 indicates that the event was observed while status = 0 indicates the observation was right-censored.
#' @param y A vector of follow-up times.
#' @param x A dataframe or matrix, including all the covariates but not treatments with rows corresponding to observations and columns to variables.
#' @param trt A numeric vector representing the treatment groups.
#' @param cluster.id A vector of integers representing the clustering id. The cluster id should be an integer and start from 1.
#' @param verbose A logical indicating whether to show the progress bar for riAFT-BART. The default is FALSE
#' @param estimand A character string representing the type of causal estimand. Only \code{"ATT"} or \code{"ATE"} is allowed. When the \code{estimand = "ATT"}, users also need to specify the reference treatment group by setting the \code{reference_trt} argument.
#' @param reference_trt A numeric value indicating reference treatment group for ATT effect.
#' @return A list of causal estimands in terms of log T between different treatment groups.
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
#' res_ate <- riAFTBART(M.burnin = 10, M.keep = 10, M.thin = 1, status = status,
#'                       y = Y, trt = trt.train,
#'                       x = cbind(x1,x2),
#'                       cluster.id = cluster.id, estimand = "ATE")
#' }
riAFTBART <- function(M.burnin, M.keep, M.thin = 1, status, y, x, trt, cluster.id, verbose = FALSE, estimand = "ATE", reference_trt = NULL){
  n_trt <- length(unique(trt))
  if (estimand == "ATE") { # causal estimand is ATE
    for (i in 1:length(unique(trt))) {
      assign(paste0("riAFTBART_ATE_",i), riAFTBART_fit(M.burnin = M.burnin,
                                                            M.keep = M.keep,
                                                            M.thin = M.thin,
                                                            status = status,
                                                            y.train = y,
                                                            x.train = x,
                                                            trt.train = trt,
                                                            x.test = x,
                                                            trt.test = i, # Counterfactual for treatment i
                                                            cluster.id = cluster.id,
                                                            verbose = verbose))
    }
    for (i in 1:(n_trt - 1)) { # Calculate the pairwise ATE effect based on the predicted counterfactual survival time in log scale
      for (j in (i + 1):n_trt) {
        assign(paste0("logT", i, j, "_est"), mean(eval(parse(text = (
          paste0("riAFTBART_ATE_", i))))[["tree.pred"]]) - mean(eval(parse(text = (
            paste0("riAFTBART_ATE_", j))))[["tree.pred"]])
        )
      }
    }
    result <- NULL # Save the results in a list
    counter <- 1 # Set the counter for naming the list
    for (i in 1:(n_trt-1)){
      for (j in (i + 1):n_trt){
        result <- c(result, (eval(parse(text = (paste0("logT", i, j, "_est")))))) #  Store the pairwise ATEs in the result list
        names(result)[[counter]] <- paste0("logT", i, j) # Naming the elements
        counter <- counter +1
      }
    }
    return(result)
  } else if (estimand == "ATT"){ # causal estimand is ATT
    trt_indicator = 1:n_trt
    trt_indicator_no_reference <- trt_indicator[trt_indicator!=reference_trt]
    for (i in 1:length(unique(trt))) {
      assign(paste0("riAFTBART_ATT_",i), riAFTBART_fit(M.burnin = M.burnin,
                                                            M.keep = M.keep,
                                                            M.thin = M.thin,
                                                            status = status,
                                                            y.train = y,
                                                            x.train = x,
                                                            trt.train = trt,
                                                            x.test = x[trt == reference_trt,],# Counterfactual for treatment i only among those in the reference treatment group
                                                            trt.test = i,
                                                            cluster.id = cluster.id,
                                                            verbose = verbose))
    }

      for (j in 1:length(trt_indicator_no_reference)){# Calculate the pairwise ATT effect based on the predicted counterfactual survival time in log scale
        assign(paste0("logT", reference_trt, trt_indicator_no_reference[j], "_est"), mean(eval(parse(text = (
          paste0("riAFTBART_ATT_", reference_trt))))[["tree.pred"]]) - mean(eval(parse(text = (
            paste0("riAFTBART_ATT_", trt_indicator_no_reference[j]))))[["tree.pred"]])
        )
      }

    result <- NULL # Save the results in a list
    counter <- 1 # Set the counter for naming the list
      for (j in 1:length(trt_indicator_no_reference)){
        result <- c(result, (eval(parse(text = (paste0("logT",reference_trt, trt_indicator_no_reference[j], "_est")))))) #  Store the pairwise ATTs in the result list
        names(result)[[counter]] <- paste0("logT", reference_trt, trt_indicator_no_reference[j]) # Naming the elements
        counter <- counter +1
      }
    return(result)
  }

}
