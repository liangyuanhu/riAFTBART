# progress bar
progressBar = function(iter, chain.len) {
  Perc = 100 * iter / chain.len
  if ( iter %% (chain.len/100) == 0 ) {
    cat("*", fill=FALSE)
    utils::flush.console()
    if ( iter %% (chain.len/10) == 0 ) {
      cat(":", Perc, "%", sep="", "\n")
      utils::flush.console()
    }
  }
}
