#' Robust Median-Based Bayesian Model Wrapper
#'
#' @param Missing_Type Type of missing mechanism ("MAR", "MNAR", "MCAR", or "no missing")
#' @param data Input data matrix
#' @param seed Random seed
#' @param chain Number of MCMC chains
#' @param Niter Total number of MCMC iterations
#' @param burnIn Number of burn-in iterations
#' @return An mcmc.list object invisibly
#' @export
RMB <- function(Missing_Type, data, seed, chain = 1, Niter = 6000, burnIn = 3000) {
  set.seed(seed)
  tau <- 0.5
  N <- nrow(data)
  Time <- ncol(data)
  dat <- list("N" = N, "y" = data, "tau" = tau, "Time" = Time)
  initial <- list(".RNG.name" = "base::Wichmann-Hill", ".RNG.seed" = seed)

  model_list <- list(
    "MNAR" = model_MNAR,
    "MAR" = model,
    "MCAR" = model,
    "no missing" = model
  )

  if (!(Missing_Type %in% names(model_list))) {
    stop("Invalid Missing_Type. Use 'MNAR', 'MAR', 'MCAR', or 'no missing'.")
  }

  model_string <- model_list[[Missing_Type]]
  message("\nRunning missing data model for Missing_Type = ", Missing_Type)
  
  result <- missing_model(dat, initial, Niter, burnIn, chain, model_string)
  return(result)
  
}
