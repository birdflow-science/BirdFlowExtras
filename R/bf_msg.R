# This is a copy of BirdFlowR::bf_msg()
# I didn't want to export it from BirdFlowR but want to use it here
# I deleted the documentation as it's already documented in BirdFlowR

bf_msg <- function(..., sep = "") {
  m <- paste(..., sep = sep)
  if (birdflow_options("verbose")) {
    cat(m)
  }
}

# Conditionally supress messages in code called from BirdFlowExtras based
# on BirdFlowR::birdflow_options("verbose")
bf_suppress_msg <- function(exp) {
  verbose <- birdflow_options("verbose")
  withCallingHandlers(
    message = function(m) {
      if (!verbose)
        tryInvokeRestart("muffleMessage")
    },
    exp
  )
}
