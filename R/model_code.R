#' @export
missing_model <- function(dat, initial, Niter, burnIn, chain = 1, model_string) {
  jags.gmm.md3 <- rjags::jags.model(textConnection(model_string), data = dat, inits = initial, n.chains = chain, n.adapt = 1000)
  params <- c("par")
  samps.gmm.md3 <- rjags::coda.samples(jags.gmm.md3, params, n.iter = Niter)

  total_iter <- nrow(as.matrix(samps.gmm.md3[[1]]))
  burnIn <- min(burnIn, total_iter - 1)
  smp.md3 <- window(samps.gmm.md3[[1]], start = burnIn)
  
  quantile_summary <- summary(smp.md3)
  geweke.md3 <- apply(smp.md3, 2, function(x) coda::geweke.diag(x)$z)
  samples_matrix <- as.matrix(samps.gmm.md3)[, 1:5]
  credible_intervals <- apply(samples_matrix, 2, quantile, probs = c(0.025, 0.975))
  hpd_intervals <- coda::HPDinterval(coda::mcmc(samples_matrix), prob = 0.95)
  output <- list(
   quantiles = quantile_summary,
    geweke = geweke.md3[1:5],
    credible_intervals = credible_intervals,
    hpd_intervals = hpd_intervals
  )
  
  return(output)
  
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
  par[6] <- sigma
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
  par[6] <- sigma
  par[7] <- r0
  par[8] <- r1
  par[9] <- r2
}"


