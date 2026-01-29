


benchmark_horseshoe <- function(
    vy,
    mX,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!requireNamespace("horseshoe", quietly = TRUE)) {
    stop(
      "The 'horseshoe' package is required for this benchmark.\n",
      "Install it with install.packages('horseshoe').",
      call. = FALSE
    )
  }

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }
  if (!is.numeric(trials) || length(trials) != 1L || trials < 1) {
    stop("'trials' must be a positive scalar.", call. = FALSE)
  }

  vy <- as.numeric(vy)
  mX <- as.matrix(mX)
  p <- ncol(mX)


  nmc =  nsamples - nburn

  # Run the Gibbs sampler trials times
  mStat = NULL
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = horseshoe::horseshoe(y = vy, X = mX, method.tau = "truncatedCauchy", method.sigma = "Jeffreys", burn = nburn, nmc = nmc)    })[3]
    # print(time_val)

    vESS <- c()
    for(j in 1:p){
      vESS[j] <- FastGibbsSamplers:::effective_sample_size(res_mcmc$BetaSamples[j,])
    }
    Ef = median(vESS)/time_val
    ESS_sigma2  = FastGibbsSamplers:::effective_sample_size(res_mcmc$Sigma2Samples)
    Ef_sigma2 = ESS_sigma2/time_val
    ESS_lambda2  = FastGibbsSamplers:::effective_sample_size(res_mcmc$TauSamples)
    Ef_lambda2 = ESS_lambda2/time_val

    stats = c(
      100*median(vESS)/nmc,
      Ef,
      100*ESS_sigma2/nmc,
      Ef_sigma2,
      100*ESS_lambda2/nmc,
      Ef_lambda2,
      time_val)

    # print(stats)

    mStat = rbind(mStat,stats)
  }

  colname_vals = c( "mix_beta", "eff_beta",  "mix_sigma2", "eff_sigma2",  "mix_lambda2", "eff_lambda2", "time")
  colnames(mStat) <- colname_vals

  return(FastGibbsSamplers:::mcmc_diagnostics(t(res_mcmc$BetaSamples), res_mcmc$Sigma2Samples, res_mcmc$TauSamples, beta_inds, mStat))
}

