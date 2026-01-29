

benchmark_AHS <- function(
    vy,
    mX,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!requireNamespace("Mhorseshoe", quietly = TRUE)) {
    stop(
      "The 'Mhorseshoe' package is required for this benchmark.\n",
      "Install it with install.packages('Mhorseshoe').",
      call. = FALSE
    )
  }

  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "The 'posterior' package is required for ESS computation in this benchmark.\n",
      "Install it with install.packages('posterior').",
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

  mX <- as.matrix(mX)
  p <- ncol(mX)

  nmc =  nsamples - nburn

  # Run the Gibbs sampler trials times
  mStat = c()
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = Mhorseshoe::approx_horseshoe(
        y = vy,
        X = mX,
        burn = 1000,
        iter = 5000,
        auto.threshold = TRUE,
        threshold = 0,
        tau = 1,
        s = 0.8,
        sigma2 = 1,
        w = 1,
        alpha = 0.05,
        a = 0.2,
        b = 10,
        t = 10,
        adapt_p0 = 0,
        adapt_p1 = -4.6 * 10^(-4)
      )
    })[3]
    # print(time_val)

    vESS <- c()
    for(j in 1:p){
      vESS[j] <- ess_basic(res_mcmc$BetaSamples[,j])
    }
    Ef = median(vESS)/time_val
    ESS_sigma2  = posterior::ess_bulk(res_mcmc$Sigma2Samples)
    Ef_sigma2 = ESS_sigma2/time_val
    ESS_lambda2  = posterior::ess_bulk(res_mcmc$TauSamples)
    Ef_lambda2 = ESS_lambda2/time_val

    N = length(res_mcmc$Sigma2Samples)


    stats = c(
      100*median(vESS)/N,
      Ef,
      100*ESS_sigma2/N,
      Ef_sigma2,
      100*ESS_lambda2/N,
      Ef_lambda2,
      time_val)

    # print(stats)

    mStat = rbind(mStat,stats)
  }

  colname_vals = c( "mix_beta", "eff_beta",  "mix_sigma2", "eff_sigma2",  "mix_lambda2", "eff_lambda2", "time")
  colnames(mStat) <- colname_vals

  return(FastGibbsSamplers:::mcmc_diagnostics(res_mcmc$BetaSamples, res_mcmc$Sigma2Samples, res_mcmc$TauSamples, beta_inds, mStat, doplots = FALSE))
}



benchmark_EHS <- function(
    vy,
    mX,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!requireNamespace("Mhorseshoe", quietly = TRUE)) {
    stop(
      "The 'Mhorseshoe' package is required for this benchmark.\n",
      "Install it with install.packages('Mhorseshoe').",
      call. = FALSE
    )
  }

  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop(
      "The 'posterior' package is required for ESS computation in this benchmark.\n",
      "Install it with install.packages('posterior').",
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

  mX <- as.matrix(mX)
  p <- ncol(mX)


  nmc =  nsamples - nburn

  # Run the Gibbs sampler trials times
  mStat = c()
  for (i in seq_len(trials)) {
    time_val = system.time({
        res_mcmc = Mhorseshoe::exact_horseshoe(
          y = vy,
          X = mX,
          burn = 1000,
          iter = 5000,
          a = 1/5,
          b = 10,
          s = 0.8,
          tau = 1,
          sigma2 = 1,
          w = 1,
          alpha = 0.05
        )
    })[3]
    # print(time_val)

    vESS <- c()
    for(j in 1:p){
      vESS[j] <- ess_basic(res_mcmc$BetaSamples[,j])
    }
    Ef = median(vESS)/time_val
    ESS_sigma2  = posterior::ess_bulk(res_mcmc$Sigma2Samples)
    Ef_sigma2 = ESS_sigma2/time_val
    ESS_lambda2  = posterior::ess_bulk(res_mcmc$TauSamples)
    Ef_lambda2 = ESS_lambda2/time_val

    N = length(res_mcmc$Sigma2Samples)


    stats = c(
      100*median(vESS)/N,
      Ef,
      100*ESS_sigma2/N,
      Ef_sigma2,
      100*ESS_lambda2/N,
      Ef_lambda2,
      time_val)

    # print(stats)

    mStat = rbind(mStat,stats)
  }

  colname_vals = c( "mix_beta", "eff_beta",  "mix_sigma2", "eff_sigma2",  "mix_lambda2", "eff_lambda2", "time")
  colnames(mStat) <- colname_vals

  return(FastGibbsSamplers:::mcmc_diagnostics(res_mcmc$BetaSamples, res_mcmc$Sigma2Samples, res_mcmc$TauSamples, beta_inds, mStat, doplots = FALSE))
}


