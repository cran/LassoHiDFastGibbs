
benchmark_pcg_lambda2_col_sigma2 <- function(
    vy,
    mX,
    penalty_type,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }
  if (!is.numeric(trials) || length(trials) != 1L || trials < 1) {
    stop("'trials' must be a positive scalar.", call. = FALSE)
  }


  # Run the Gibbs sampler trials times
  mStat = NULL
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = penalized_pcg_lambda2_sigma2(vy, mX,
                                                     penalty_type,
                                                     a,
                                                     b,
                                                     u,
                                                     v,
                                                     nsamples,
                                                     lambda_init,
                                                     sigma2_init,
                                                     verbose=nsamples/5)
      })[3]
    # print(time_val)

    inds_use = (nburn + 1):nsamples

    # Calculate summary statistics of efficiencies and mixing rates
    stats = mcmc_stats(res_mcmc$mBeta, res_mcmc$vsigma2, res_mcmc$vlambda2, time_val, inds_use)
    # print(stats)

    mStat = rbind(mStat,stats)
  }

  return(mcmc_diagnostics(res_mcmc$mBeta[inds_use,], res_mcmc$vsigma2[inds_use], res_mcmc$vlambda2[inds_use], beta_inds, mStat))
}

################################################################################

benchmark_pcg_sigma2_col_lambda2 <- function(
    vy,
    mX,
    penalty_type,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }
  if (!is.numeric(trials) || length(trials) != 1L || trials < 1) {
    stop("'trials' must be a positive scalar.", call. = FALSE)
  }

  # Run the Gibbs sampler trials times
  mStat = NULL
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = penalized_pcg_sigma2_lambda2(vy,
                                                     mX,
                                                     penalty_type,
                                                     a,
                                                     b,
                                                     u,
                                                     v,
                                                     nsamples,
                                                     lambda_init,
                                                     va_init = rep(1,p),
                                                     verbose = nsamples/5,
                                                     lower = 1.0E-12,
                                                     upper = 5000)
    })[3]
    # print(time_val)

    inds_use = (nburn + 1):nsamples

    # Calculate summary statistics of efficiencies and mixing rates
    stats = mcmc_stats(res_mcmc$mBeta, res_mcmc$vsigma2, res_mcmc$vlambda2, time_val, inds_use)
    # print(stats)

    mStat = rbind(mStat,stats)
  }

  return(mcmc_diagnostics(res_mcmc$mBeta[inds_use,], res_mcmc$vsigma2[inds_use], res_mcmc$vlambda2[inds_use], beta_inds, mStat))
}


################################################################################

benchmark_blasso_pcg_lambda2_col_va <- function(
    vy,
    mX,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }
  if (!is.numeric(trials) || length(trials) != 1L || trials < 1) {
    stop("'trials' must be a positive scalar.", call. = FALSE)
  }

  # Run the Gibbs sampler trials times
  mStat = NULL
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = blasso_pcg_lambda2_va(vy,
                                           mX,
                                           a,
                                           b,
                                           u,
                                           v,
                                           nsamples,
                                           lambda_init,
                                           sigma2_init,
                                           verbose=nsamples/5)
    })[3]
    # print(time_val)

    inds_use = (nburn + 1):nsamples

    # Calculate summary statistics of efficiencies and mixing rates
    stats = mcmc_stats(res_mcmc$mBeta, res_mcmc$vsigma2, res_mcmc$vlambda2, time_val, inds_use)
    # print(stats)

    mStat = rbind(mStat,stats)
  }

  return(mcmc_diagnostics(res_mcmc$mBeta[inds_use,], res_mcmc$vsigma2[inds_use], res_mcmc$vlambda2[inds_use], beta_inds, mStat))

}

################################################################################

benchmark_blasso_pcg_sigma2_col_va <- function(
    vy,
    mX,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }
  if (!is.numeric(trials) || length(trials) != 1L || trials < 1) {
    stop("'trials' must be a positive scalar.", call. = FALSE)
  }


  # Run the Gibbs sampler trials times
  mStat = NULL
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = blasso_pcg_sigma2_va(vy,
                                                mX,
                                                a,
                                                b,
                                                u,
                                                v,
                                                nsamples,
                                                lambda_init,
                                                sigma2_init,
                                                va_init = rep(1,p),
                                                verbose = nsamples/5,
                                                lower = 1.0E-12,
                                                upper = 5000)
      })[3]
    # print(time_val)

    inds_use = (nburn + 1):nsamples

    # Calculate summary statistics of efficiencies and mixing rates
    stats = mcmc_stats(res_mcmc$mBeta, res_mcmc$vsigma2, res_mcmc$vlambda2, time_val, inds_use)
    # print(stats)

    mStat = rbind(mStat,stats)
  }

  return(mcmc_diagnostics(res_mcmc$mBeta[inds_use,], res_mcmc$vsigma2[inds_use], res_mcmc$vlambda2[inds_use], beta_inds, mStat))
}




################################################################################

benchmark_pcg_vbeta_col_sigma2 <- function(
    vy,
    mX,
    penalty_type,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }
  if (!is.numeric(trials) || length(trials) != 1L || trials < 1) {
    stop("'trials' must be a positive scalar.", call. = FALSE)
  }


  # Run the Gibbs sampler trials times
  mStat = NULL
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = penalized_pcg_beta_sigma2(vy,
                                                   mX,
                                                   penalty_type,
                                                   a,
                                                   b,
                                                   u,
                                                   v,
                                                   nsamples,
                                                   lambda_init,
                                                   sigma2_init,
                                                   verbose=nsamples/5)
    })[3]
    # print(time_val)

    inds_use = (nburn + 1):nsamples

    # Calculate summary statistics of efficiencies and mixing rates
    stats = mcmc_stats(res_mcmc$mBeta, res_mcmc$vsigma2, res_mcmc$vlambda2, time_val, inds_use)
    # print(stats)

    mStat = rbind(mStat,stats)
  }

  return(mcmc_diagnostics(res_mcmc$mBeta[inds_use,], res_mcmc$vsigma2[inds_use], res_mcmc$vlambda2[inds_use], beta_inds, mStat))
}

################################################################################

benchmark_pcg_sigma2_col_vbeta <- function(
    vy,
    mX,
    penalty_type,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }
  if (!is.numeric(trials) || length(trials) != 1L || trials < 1) {
    stop("'trials' must be a positive scalar.", call. = FALSE)
  }


  # Run the Gibbs sampler trials times
  mStat = NULL
  for (i in seq_len(trials)) {
    time_val = system.time({
      res_mcmc = penalized_pcg_sigma2_beta(vy,
                                                   mX,
                                                   penalty_type,
                                                   a,
                                                   b,
                                                   u,
                                                   v,
                                                   nsamples,
                                                   lambda_init,
                                                   va_init = rep(1,p),
                                                   verbose=nsamples/5,
                                                   lower = 1.0E-12,
                                                   upper = 5000)
    })[3]
    # print(time_val)

    inds_use = (nburn + 1):nsamples

    # Calculate summary statistics of efficiencies and mixing rates
    stats = mcmc_stats(res_mcmc$mBeta, res_mcmc$vsigma2, res_mcmc$vlambda2, time_val, inds_use)
    # print(stats)

    mStat = rbind(mStat,stats)
  }

  return(mcmc_diagnostics(res_mcmc$mBeta[inds_use,], res_mcmc$vsigma2[inds_use], res_mcmc$vlambda2[inds_use], beta_inds, mStat))
}


