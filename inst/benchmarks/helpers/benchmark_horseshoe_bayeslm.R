

benchmark_horseshoe_bayeslm <- function(
    vy,
    mX,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!requireNamespace("bayeslm", quietly = TRUE)) {
    stop(
      "The 'bayeslm' package is required for this benchmark.\n",
      "Install it with install.packages('bayeslm').",
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

  block_vec = rep(1, p)

  # Run the Gibbs sampler trials times
  mStat = NULL
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = bayeslm::bayeslm(Y=vy, X=mX, prior = 'horseshoe', icept = FALSE, singular = TRUE,
                         block_vec = block_vec, N =(nsamples - nburn), burnin=nburn)
    })[3]
    # print(time_val)


    inds_use = seq_len(nsamples - nburn)

    # Calculate summary statistics of efficiencies and mixing rates
    stats = FastGibbsSamplers:::mcmc_stats(res_mcmc$beta, res_mcmc$sigma, res_mcmc$vglobal, time_val, inds_use)
    # print(stats)

    mStat = rbind(mStat,stats)
  }

  colname_vals = c("eff_beta", "mix_beta", "eff_sigma2", "mix_sigma2", "eff_lambda2", "mix_lambda2", "time")
  colnames(mStat) <- colname_vals

  return(FastGibbsSamplers:::mcmc_diagnostics(res_mcmc$beta, res_mcmc$sigma, res_mcmc$vglobal, beta_inds, mStat, doplots=FALSE))

}
