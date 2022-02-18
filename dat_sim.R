#' Simulate data with multiple treatments and clustered survival outcomes
#'
#' This function simulate data with multiple treatments and clustered survival outcomes. Users can adjust the following 11 design factors: (1) The number of clusters, (2) the sample size in each cluster, (3) ratio of units across treatment groups, (4) whether the treatment assignment model and the outcome generating model are linear or nonlinear, (5) whether the covariates that best predict the treatment also predict the outcome well, (6) whether the response surfaces are parallel across treatment groups, (7) degree of covariate overlap, (8) Whether the proportional hazards assumption is satisfied, (9) mean follow up time for each treatment group, (10) censoring proportion and (11) Standard deviation for the cluster effect in the treatment assignment and outcome generating model.
#'
#' @param nK A numeric value indicating the number of clusters.
#' @param K A numeric value indicating the sample size in each cluster.
#' @param n_trt A numeric value indicating the number of treatments.
#' @param X A vector of characters representing covariates, with each covariate being generated from the standard probability \code{\link[stats:Distributions]{distributions}} in the \code{\link[stats:stats-package]{stats}} package.
#' @param lp_y A vector of characters of length \code{n_trt}, representing the linear effects in the outcome generating model.
#' @param nlp_y A vector of characters of length \code{n_trt}, representing the nonlinear effects in the outcome generating model.
#' @param align A logical indicating whether the predictors in the treatment assignment model are the same as the predictors for the outcome generating model. The default is \code{TRUE}. If the argument is set to \code{FALSE}, users need to specify additional two arguments \code{lp_w} and \code{nlp_w}.
#' @param eta A numeric value to induce proportional hazards assumption or a character including linear combination of Xs to induce nonproportional hazards assumption.
#' @param lambda A numeric vector of length \code{n_trt} inducing different follow up time across treatment groups.
#' @param delta A numeric vector of length \code{n_trt}-1 inducing different ratio of units across treatment groups.
#' @param psi A numeric value for the parameter governing the sparsity of covariate overlap.
#' @param lp_w A vector of characters of length \code{n_trt} - 1, representing the treatment assignment model.
#' @param nlp_w A vector of characters of length \code{n_trt} - 1, representing the treatment assignment model.
#' @param censor_rate A numeric value for the rate parameter governing the proportion of censoring.
#' @param sigma_w A numeric value representing the standard deviation for the cluster effect in the treatment assignment model.
#' @param sigma_y A numeric value representing the standard deviation for the cluster effect in the outcome generating model.
#'
#' @return A list with 7 elements for simulated data. It contains
#' \item{covariates:}{X matrix}
#' \item{w:}{treatment indicators}
#' \item{Tobs:}{observed follow up time for the simulated right censored data}
#' \item{status:}{the censoring indicator}
#' \item{cluster:}{the clustering indicator}
#' \item{censor_prop:}{the censoring proportion}
#' \item{T_mean:}{mean observed follow up time}
#' \item{ratio_of_units:}{the proportions of units in each treatment group}
#' @export
#' @import dplyr
#'
#' @examples
#' library(riAFTBART)
#' lp_w_all <-
#'   c(".4*x1 + .1*x2  - .1*x4 + .1*x5",    # w = 1
#'     ".2 * x1 + .2 * x2  - .2 * x4 - .3 * x5")  # w = 2
#' nlp_w_all <-
#'   c("-.5*x1*x4  - .1*x2*x5", # w = 1
#'     "-.3*x1*x4 + .2*x2*x5")# w = 2
#' lp_y_all <- rep(".2*x1 + .3*x2 - .1*x3 - .1*x4 - .2*x5", 3)
#' nlp_y_all <- rep(".7*x1*x1  - .1*x2*x3", 3)
#' X_all <- c(
#'   "rnorm(1000, 0, 0.5)",# x1
#'   "rbeta(1000, 2, .4)",   # x2
#'   "runif(1000, 0, 0.5)",# x3
#'   "rweibull(1000,1,2)",  # x4
#'   "rbinom(1000, 1, .4)"# x5
#' )
#' set.seed(111111)
#' data <- dat_sim(
#'   nK = 20,
#'   K = 50,
#'   n_trt = 3,
#'   X = X_all,
#'   eta = 2,
#'   lp_y = lp_y_all,
#'   nlp_y  = nlp_y_all,
#'   align = FALSE,
#'   lp_w = lp_w_all,
#'   nlp_w = nlp_w_all,
#'   lambda = c(1000,2000,3000),
#'   delta = c(0.5,0.5),
#'   psi = 1,
#'   sigma_w = 1,
#'   sigma_y = 2,
#'   censor_rate = 0.1
#' )
dat_sim <- function(nK,
                    K,
                    n_trt,
                    X,
                    lp_y,
                    nlp_y,
                    align = TRUE,
                    eta,
                    lambda,
                    delta,
                    psi,
                    lp_w,
                    nlp_w,
                    sigma_w,
                    sigma_y,
                    censor_rate)
{
  if (align == TRUE) { # If align = TRUE, the treatment assignment model will be the same as the predictors for the outcome generating model
    lp_w <- lp_y
    nlp_w <- nlp_y
  }
  if (is.null(lp_y)) { # If lp_y is null, we set it as 0s
    lp_y <- rep(0, n_trt)
  }
  if (is.null(nlp_y)) { # If nlp_y is null, we set it as 0s
    nlp_y <- rep(0, n_trt)
  }
  if (is.null(lp_w)) { # If lp_w is null, we set it as 0s
    lp_w <- rep(0, n_trt)
  }
  if (is.null(nlp_w)) { # If nlp_w is null, we set it as 0s
    nlp_w <- rep(0, n_trt)
  }
  sample_size <- nK * K # Total sample size is number of clusters times the sample size in each cluster
  for (i in 1:length(X)){
    assign(paste0("x",i), eval(parse(text = X[i]))) # assign xs to the distribution users specified
  }

  X_matrix <- matrix(NA, nrow = sample_size, ncol = length(X))
  for (i in 1:dim(X_matrix)[2]){
    X_matrix[,i] <- eval(parse(text = paste0("x",i))) # Create an X matrix to store all the generated Xs
  }

  xi = stats::rnorm(K, 0,sigma_w^2) # Generate random intercepts for the treatment assignment model

  cl = rep(1:K, each=nK) # Generate cluster labels
  xiX  = rep(xi, each=nK) # Generate random intercepts of the treatment assignment model for all the individuals

  treatment_exp_matrix_all <- NULL
  for (i in 1:(n_trt-1)){ # Generate treatment exp values based on psi, delta, lp_w, nlp_w and xiX
    treatment_exp_matrix <- exp(eval(parse(text = paste0( psi, "*(",delta[i], "+", lp_w[i], "+", nlp_w[i], "+xiX", ")"))))
    treatment_exp_matrix_all <- cbind(treatment_exp_matrix_all, treatment_exp_matrix)
  }
  treatment_exp_matrix_all[is.infinite(treatment_exp_matrix_all)] <- 10^8 # Cap the treatment exp values to 10^8
  probs <- sweep(treatment_exp_matrix_all, 1, (rowSums(treatment_exp_matrix_all)+1), '/') # Calculate the assignment probability for each treatment group except the last treatment group

  prob_last_all <- rep(NA, dim(probs)[1])
  for (i in 1:dim(probs)[1]){
    prob_last_all[i] <- 1-sum(probs[i,]) # Calculate the assignment probability for the last treatment group
  }
  w <- rep(NA, dim(probs)[1])
  for (i in 1:dim(probs)[1]){ # Assign treatment to each individual based on the treatment assignment probability for each individual
    w[i] <- sample(x = 1:n_trt, size = 1, prob = c(probs[i,], 1-sum(probs[i,])), replace = TRUE)
  }
  w
  bk = stats::rnorm (K,0,sigma_y^2) # Generate random intercepts for the outcome generating model
  bkT = rep(bk, each = nK) # Generate random intercepts of the outcome generating model for all the individuals

  U = stats::runif(sample_size,0,1) # Generate U ~ unif(0,1) used for the outcome generating model
  for (i in 1:n_trt){
    assign(paste0("LP",i,"_final"),NULL)
    assign(paste0("T",i,"_final"),NULL)
  }

  for (j in 1:n_trt){ # Generate the outcome survival time based on lp_y, nlp_y, bkT, lambda and eta
    assign(paste0("LP",j),  eval(parse(text = paste0(lp_y[j], "+", nlp_y[j], "+bkT"))))
    assign(paste0("T",j),  eval(parse(text = paste0("(",lambda[j], "*(-log(U))/exp(",paste0("LP",j),"))","^(1/eta)" ))))
  }
  T_true <- matrix(NA, nrow = sample_size, ncol = n_trt)
  for (i in 1:n_trt){ # Put the true survival time for each individual into a matrix
    T_true[,i] <- eval(parse(text = paste0("T",i)))
  }
  T_true_with_treatment <- cbind(T_true, w) # cbind to get the matrix for both the true survival time and treatment indicator
  Tobs = apply(T_true_with_treatment, 1, function(x) x[1:n_trt][x[n_trt+1]])  # Observed survival time for each individual
  C <- stats::rexp(sample_size, rate = censor_rate) # Generate the censoring time for each individual
  Tobs_C <- pmin(Tobs,C) # Generate the survival/censoring time for each individual
  censor_rate <- sum(Tobs>C)/sample_size # Calculate the censoring proportion
  delta = as.numeric(Tobs>C) # Calculate the censoring indicator
  T_mean <- tibble(w = as.character(w), Tobs_C) %>% # Generate the mean follow up time for each treatment group
    group_by(w) %>%
    summarise(T_mean = mean(Tobs_C)) %>%
    bind_rows(tibble(w = "Overall",
                     T_mean = mean(Tobs_C))) %>%
    dplyr::mutate(T_mean = round(T_mean, 2)) %>%
    as.data.frame()
  return( # Output a list of values
    list(
      covariates = as.data.frame(X_matrix),
      w = w,
      Tobs = Tobs,
      T_mean = T_mean,
      censor_prop = censor_rate,
      delta = delta,
      cluster = cl,
      ratio_of_units = round(table(w) / length(w),2)
    )
  )

}
