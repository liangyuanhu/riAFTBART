blocked.mcmc.update <- function(status, y.train, x.train, x.test, cluster.id, cur,SA = SA,trt.train,prior_c_function_used,gps) {
  # First update the tree structures
  cur <- sample.tree(status, y.train, x.train, x.test, cluster.id, cur,SA = SA,trt.train,prior_c_function_used,gps)
  # Next update the random intercepts
  cur <- sample.b(status, y.train, x.train, x.test, cluster.id, cur)
  # Next update the standard deviation for the random intercept
  cur <- sample.tau(status, y.train, x.train, x.test, cluster.id, cur)
  # Next update the parameter expansion term
  cur <- sample.alpha(status, y.train, x.train, x.test, cluster.id, cur)
  # Next update the log survival time
  cur <- sample.time.censored(status, y.train, x.train, cluster.id, cur)
  return(cur)
}
