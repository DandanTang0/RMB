#' Robust Median-Based Bayesian Growth Curve Modeling
#'
#' Fits a median-based Bayesian growth curve model under MCAR, MAR,
#' MNAR or complete-data assumptions.  If \code{K > 0} the first
#' \code{K} columns in \code{data} are treated as auxiliary variables.
#'
#' @param Missing_Type Character; one of \code{MNAR}, \code{MAR},
#'   \code{MCAR}, \code{no missing}.
#' @param data Matrix or data frame containing outcome columns (and
#'   optionally auxiliary variables).
#' @param Time Integer; number of measurement occasions.
#' @param seed Integer seed for reproducibility.
#' @param K Integer; number of auxiliary variables (default 0).
#' @param chain Integer; number of MCMC chains (default 1).
#' @param Niter Integer; iterations per chain (default 6000).
#' @param burnIn Integer; burn-in iterations (default 3000).
#'
#' @return An object of class \code{RomebResult} containing
#'   \describe{
#'     \item{quantiles}{posterior means, SDs and quantiles}
#'     \item{geweke}{Geweke \emph{z}-scores}
#'     \item{credible_intervals}{95% equal-tail credible intervals}
#'     \item{hpd_intervals}{95% highest posterior density intervals}
#'     \item{samps_full}{full \code{coda::mcmc.list} (including burn-in)}
#'   }
#' @export
#'
#' @examples
#'   set.seed(123)
#'   Y <- matrix(rnorm(300), 100, 3)
#'   fit <- Romeb("no missing", data = Y, Time = 3, seed = 123, K = 0,
#'                Niter = 6000, burnIn = 3000)
#'   print(fit)



Romeb <- function(Missing_Type,
                  data,
                  Time,
                  seed,
                  K      = 0,
                  chain  = 1,
                  Niter  = 6000,
                  burnIn = 3000) {
  
  
  if (!Missing_Type %in% c("MNAR","MAR","MCAR","no missing"))
    stop("Missing_Type must be 'MNAR','MAR','MCAR', or 'no missing'.")
  
  if (!is.matrix(data) && !is.data.frame(data))
    stop("data must be a matrix or data.frame.")
  
  if (Time <= 0 || Time != floor(Time))
    stop("Time must be a positive integer.")
  
  if (K < 0 || K != floor(K))
    stop("K must be a non-negative integer.")
  
  if (K + Time > ncol(data))
    stop("data does not have enough columns for K + Time variables.")
  
  
  set.seed(seed)
  tau <- 0.5
  N   <- nrow(data)
  initial <- list(".RNG.name" = "base::Wichmann-Hill",
                  ".RNG.seed" = seed)
  
  if (K == 0) {
    dat <- list(N = N, y = data, tau = tau, Time = Time)
    
    model_string <- switch(Missing_Type,
                           "MNAR"       = model_MNAR,
                           "MAR"        = model,
                           "MCAR"       = model,
                           "no missing" = model)
    message("\nRunning ", Missing_Type, " model  (K = 0)")
  } else {
    X_part <- data[, 1:K, drop = FALSE]
    Y_part <- data[, (K+1):(K+Time), drop = FALSE]
    dat <- list(N = N, y = Y_part, tau = tau,
                K = K, X = X_part, Time = Time)
    
   
    model_string <- model_MNAR_k
    message("\nRunning ", Missing_Type,
            " model  (K = ", K, " auxiliary vars)")
  }
  
  
  jm <- rjags::jags.model(textConnection(model_string),
                          data   = dat,
                          inits  = initial,
                          n.chains = chain,
                          n.adapt  = 1000)
  
  samps_full <- rjags::coda.samples(jm, c("par"), n.iter = Niter)
  total_iter <- nrow(as.matrix(samps_full[[1]]))
  burnIn     <- min(burnIn, total_iter - 1)
  smp_post   <- window(samps_full[[1]], start = burnIn)
  
  
  quantile_summary <- summary(smp_post)
  geweke_z         <- apply(smp_post, 2,
                            function(x) coda::geweke.diag(x)$z)
  credible_int     <- apply(as.matrix(samps_full),
                            2, quantile, probs = c(.025,.975))
  hpd_int          <- coda::HPDinterval(coda::mcmc(as.matrix(samps_full)),
                                        prob = .95)
  
  
  output <- list(
    quantiles          = quantile_summary,
    geweke             = geweke_z,
    credible_intervals = credible_int,
    hpd_intervals      = hpd_int,
    samps_full         = samps_full      
  )
  class(output) <- "RomebResult"         
  return(output)                         
}

#' @export
print.RomebResult <- function(x, ...) {
  cat("Romeb GCM summary\n",
      "==================\n", sep = "")
  cat("\nPosterior medians (50% quantiles):\n")
  print(x$quantiles$quantiles[,"50%"])
  cat("\nGeweke test:\n")
  print(x$geweke)
  cat("\n95% credible intervals:\n")
  print(t(x$credible_intervals))
  cat("\n95% hpd intervals:\n")
  print(t(x$hpd_intervals[,1:2]))
  cat("\nUse x$samps_full to access full MCMC samples,",
      "and coda::traceplot(x$samps_full[,'par[i]']) for the trace plot of par[i].\n")
  invisible(x)
}


#' @export
model <- "model {
  for (i in 1:N)  {
    for(t in 1:Time) {
      V[i,t] ~ dexp(pre_sigma)
      y[i,t] ~ dnorm(muy[i,t], pre_sig2[i,t])
      muy[i,t] <- LS[i,1]+(t-1)*LS[i,2] + zeta*V[i,t]
      pre_sig2[i,t]<- 1/sig2_y[i,t]
      sig2_y[i,t] <- eta^2*V[i,t]/pre_sigma
      loglik[i,t] <- logdensity.norm(y[i,t], muy[i,t], pre_sig2[i,t])
    }
    LS[i,1:2]  ~ dmnorm(muLS[1:2], Inv_cov[1:2,1:2])
  }
  zeta <- (1-2*tau)/(tau*(1-tau))
  eta <- sqrt(2/(tau*(1-tau)))
  
  pre_sigma ~ dgamma(.001, .001)
  sigma <- 1/pre_sigma

  muLS[1] ~ dnorm(0, 0.001)
  muLS[2] ~ dnorm(0, 0.001)

  Inv_cov[1:2,1:2] ~ dwish(R[1:2,1:2], 3)
  Cov_b <- inverse(Inv_cov[1:2,1:2])
  
  R[1,1] <- 1
  R[2,2] <- 1
  R[2,1] <- R[1,2]
  R[1,2] <- 0

  par[1] <- muLS[1]
  par[2] <- muLS[2]
  par[3] <- Cov_b[1,1]
  par[4] <- Cov_b[1,2]
  par[5] <- Cov_b[2,2]
}"

#' @export
model_MNAR <- "model {
  for (i in 1:N)  {
    for(t in 1:Time) {
      V[i,t] ~ dexp(pre_sigma)
      y[i,t] ~ dnorm(muy[i,t], pre_sig2[i,t])
      muy[i,t] <- LS[i,1]+(t-1)*LS[i,2] + zeta*V[i,t]
      pre_sig2[i,t]<- 1/sig2_y[i,t]
      sig2_y[i,t] <- eta^2*V[i,t]/pre_sigma
      loglik[i,t] <- logdensity.norm(y[i,t], muy[i,t], pre_sig2[i,t])
    }
    LS[i,1:2]  ~ dmnorm(muLS[1:2], Inv_cov[1:2,1:2])
  }

  zeta <- (1-2*tau)/(tau*(1-tau))
  eta <- sqrt(2/(tau*(1-tau)))

  for(i in 1:N){
    for(t in 2:Time){
      m[i,t] ~ dbern(q[i,t])
      logit(q[i,t]) <-  r0 + r1*y[i,(t-1)] + r2*y[i,t]
    }
  }

  r0  ~ dnorm(0, 0.001)
  r1  ~ dnorm(0, 0.001)
  r2  ~ dnorm(0, 0.001)

  pre_sigma ~ dgamma(.001, .001)
  sigma <- 1/pre_sigma

  muLS[1] ~ dnorm(0, 0.001)
  muLS[2] ~ dnorm(0, 0.001)

  Inv_cov[1:2,1:2] ~ dwish(R[1:2,1:2], 3)
  Cov_b <- inverse(Inv_cov[1:2,1:2])
  
  R[1,1] <- 1
  R[2,2] <- 1
  R[2,1] <- R[1,2]
  R[1,2] <- 0

  par[1] <- muLS[1]
  par[2] <- muLS[2]
  par[3] <- Cov_b[1,1]
  par[4] <- Cov_b[1,2]
  par[5] <- Cov_b[2,2]
  par[6] <- r0
  par[7] <- r1
  par[8] <- r2
}"

#' @export
model_MNAR_k <- "
model {
  for (i in 1:N) {
    for(t in 1:Time) {
      V[i,t] ~ dexp(pre_sigma)
      
      y[i,t] ~ dnorm(muy[i,t], pre_sig2[i,t])
      
      muy[i,t]      <- LS[i,1] + (t-1)*LS[i,2] + zeta*V[i,t]
      pre_sig2[i,t] <- 1/sig2_y[i,t]
      sig2_y[i,t]   <- eta^2 * V[i,t] / pre_sigma
      
      loglik[i,t]   <- logdensity.norm(y[i,t], muy[i,t], pre_sig2[i,t])
    }
    
    LS[i,1:2]  ~ dmnorm(muLS[1:2], Inv_cov[1:2,1:2])
  }

  for(i in 1:N){
    for(t in 2:Time){
      m[i,t] ~ dbern(q[i,t])
      # logit(q[i,t]) = r0 + r1*y[i,(t-1)] + r2*y[i,t] + inprod( beta_aux[], X[i,] )
      logit(q[i,t]) <- r0 + r1*y[i,(t-1)] + r2*y[i,t] + inprod(beta_aux[], X[i,])
    }
  }

  zeta <- (1 - 2*tau) / (tau*(1 - tau))
  eta  <- sqrt(2 / (tau*(1 - tau)))

  r0 ~ dnorm(0, 0.001)
  r1 ~ dnorm(0, 0.001)
  r2 ~ dnorm(0, 0.001)
  
  for(k in 1:K){
    beta_aux[k] ~ dnorm(0, 0.001)
  }

  pre_sigma ~ dgamma(0.001, 0.001)
  sigma <- 1 / pre_sigma

  muLS[1] ~ dnorm(0, 0.001)
  muLS[2] ~ dnorm(0, 0.001)
  
  Inv_cov[1:2,1:2] ~ dwish(R[1:2,1:2], 3)
  Cov_b <- inverse(Inv_cov[1:2,1:2])
  
  R[1,1] <- 1
  R[2,2] <- 1
  R[1,2] <- 0
  R[2,1] <- R[1,2]

  par[1] <- muLS[1]
  par[2] <- muLS[2]
  par[3] <- Cov_b[1,1]
  par[4] <- Cov_b[1,2]
  par[5] <- Cov_b[2,2]
  par[6] <- r0
  par[7] <- r1
  par[8] <- r2

  for(k in 1:K){
    par[8 + k] <- beta_aux[k]
  }
}
"

