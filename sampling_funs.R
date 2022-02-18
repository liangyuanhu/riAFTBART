# update tree structure, BART(X) and residual standard error
sample.tree <- function(status, y.train, x.train, x.test, cluster.id, old, trt.train,SA = SA,prior_c_function_used,gps) {
  logT <- old$logT
  b <- old$b
  # Fit BART using logT-b as the outcome
  Z <- logT - b[cluster.id]
  if (SA){ # Sensitivity analysis to correct for Z if sensitivity analysis is called for
    n_trt <- length(unique(trt.train))
    Z = ifelse(trt.train == sort(unique(trt.train))[1], Z - (unlist(prior_c_function_used[1]) * gps[2, ] + unlist(prior_c_function_used[4]) * gps[3, ]), # Scenario 1: Treatment == 1
                         ifelse(trt.train == sort(unique(trt.train))[2], Z - (unlist(prior_c_function_used[2]) * gps[1, ] + unlist(prior_c_function_used[3]) * gps[3, ]),# Scenario 2: Treatment == 2
                                Z - (unlist(prior_c_function_used[5]) * gps[1, ] + unlist(prior_c_function_used[6]) * gps[2, ]))) # Scenario 3: Treatment == 3
  }
  # Fit BART using dbarts::bart function
  mod <- dbarts::bart(x.train=x.train, y.train=Z, x.test=x.test, verbose=FALSE)
  tree <- mod$yhat.train.mean # Update using posterior mean of the tree value from BART using training data
  sigma <- mean(mod$sigma) # Update using posterior mean of sigma value from BART
  tree.pred <- mod$yhat.test.mean # Update using posterior mean of the tree value from BART using test data
  cur <- old
  cur$tree <- tree # Assign updated tree value to the current list
  cur$sigma <- sigma # Assign updated sigma value to the current list
  cur$tree.pred <- tree.pred # Assign updated tree value using test data to the current list
  cur$SA <- SA # Assign SA value, which is the same throughout
  cur$vip <- apply(mod$varcount,2,mean)/200 # Mean VIP from the fitted BART (the total count of the number of times that variable is used in a tree decision rule divided by ntree, which is 200 by default)
  return(cur)
}



# Update random intercept, b based on derived formula; email jj869@sph.rutgers.edu to get the complete derivation.
sample.b <- function(status, y.train, x.train, x.test, cluster.id, old) {
  n <- length(unique(cluster.id))
  cluster.size <- as.numeric(table(cluster.id))
  logT <- old$logT
  tau2 <- old$tau^2
  alpha <- old$alpha
  sigma2 <- old$sigma^2
  tree <- old$tree
  b.new <- NULL
  for (k in 1:n) {
    b.var <- (cluster.size[k]/sigma2 + 1/(tau2 * alpha))^(-1)
    yy <- logT[cluster.id==k]
    tt <- tree[cluster.id==k]
    resid <- yy - tt
    b.mu <- b.var * sum(resid) / sigma2
    b.new[k] <- stats::rnorm(1, b.mu, sqrt(b.var))
  }
  cur <- old
  cur$b <- b.new
  return(cur)
}

# Update parameter expansion term based on derived formula; email jj869@sph.rutgers.edu to get the complete derivation.
sample.alpha <- function(status, y.train, x.train, x.test, cluster.id, old) {
  b <- old$b
  tau2 <- old$tau^2
  n <- length(unique(cluster.id))
  rss <- t(b) %*% b
  alpha <- MCMCpack::rinvgamma(1, shape = 1, scale = rss/(2 * tau2) + 1)
  cur <- old
  cur$alpha <- sqrt(alpha)
  return(cur)
}



# Update standard deviation of random intercept, tau based on derived formula; email jj869@sph.rutgers.edu to get the complete derivation.
sample.tau <- function(status, y.train, x.train, x.test, cluster.id, old) {
  hyperpars = list(d1=1, d2=1)
  b <- old$b
  n <- length(unique(cluster.id))
  alpha <- old$alpha
  d1 <- hyperpars$d1
  d2 <- hyperpars$d2
  rss <- t(b) %*% b
  tau2 <- MCMCpack::rinvgamma(1, shape = n/2 + d1, scale = rss/(2*alpha) + d2)
  cur <- old
  cur$tau <- sqrt(tau2)
  return(cur)
}



# Impute censored event time from a truncated normal distribution based on derived formula
sample.time.censored <- function(status, y.train, x.train, cluster.id, old) {
  lower.bound <- log(y.train)
  c.id <- which(status==0) # ID for those who have censored
  n.censored <- length(c.id) # Number of individuals who have censored
  N <- nrow(x.train) # Sample size
  b <- old$b # Random intercept
  tree <- old$tree
  sigma <- old$sigma
  logT <- old$logT
  logT.mean <- tree + b[cluster.id]
  logT.censored <- NULL
  # Censored logT imputed from the a truncated normal distribution
  for(c in 1:n.censored) logT.censored[c] <- msm::rtnorm(1, mean=logT.mean[c.id[c]], sd=sigma, lower=lower.bound[c.id[c]], upper=Inf)
  logT[status==0] <- logT.censored # Only update those who have censored
  cur <- old
  cur$logT <- logT
  return(cur)
}
