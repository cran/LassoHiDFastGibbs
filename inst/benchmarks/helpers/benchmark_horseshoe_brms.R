

benchmark_horseshoe_brms <- function(
    vy,
    mX,
    lambda_init,
    sigma2_init,
    a, b, u, v,
    nburn, nsamples,
    trials,
    beta_inds=NA)
{

  if (!requireNamespace("brms", quietly = TRUE)) {
    stop(
      "The 'brms' package is required for this benchmark.\n",
      "Install it with install.packages('brms').",
      call. = FALSE
    )
  }

  if (!is.numeric(nsamples) || !is.numeric(nburn) || length(nsamples) != 1L || length(nburn) != 1L) {
    stop("'nsamples' and 'nburn' must be numeric scalars.", call. = FALSE)
  }
  if (nsamples <= nburn) {
    stop("'nsamples' must be greater than 'nburn'.", call. = FALSE)
  }



  vy <- as.numeric(vy)
  mX <- as.matrix(mX)
  p <- ncol(mX)

  prior_hs <- brms::prior(brms::horseshoe(df = 1), class = "b")

  warmup = nburn
  samples = nsamples - nburn

  fit <- brms::brm(
    formula = vy ~ -1+.,  # Include all predictors
    data = data.frame(vy=vy,mX=mX),
    family = gaussian(),
    prior = prior_hs,
    chains = 1,
    iter = nburn+1,
    warmup = nburn,
    cores = 1,
    control = list(adapt_delta = 0.95)  # To improve sampling efficiency
  )

  # Run the Gibbs sampler trials times

    # Initial fit

    time_val = system.time({
      res_mcmc <- brms::brm(
        formula = vy ~ -1+.,  # Include all predictors
        data = data.frame(vy=vy,mX=mX),
        family = gaussian(),
        prior = prior_hs,
        chains = 1,
        iter = samples,
        warmup = nburn,
        cores = 1,
        control = list(adapt_delta = 0.95)  # To improve sampling efficiency
      )
    })[3]
    # print(time_val)

    samples <- brms::as_draws_df(res_mcmc)
    beta_samples <- samples[, grepl("^b_", colnames(samples))]
    sigma_samples <- samples[, "sigma"]
    lambda_samples <- samples[, "hs_global"]

    sigma_samples = as.vector(sigma_samples$sigma)
    lambda_samples = as.vector(lambda_samples$hs_global)

    N = length(lambda_samples)

    vESS <- c()
    for(j in 1:p){
      vESS[j] <- FastGibbsSamplers:::effective_sample_size(beta_samples[[j]])
    }
    Ef = median(vESS)/time_val
    ESS_sigma2  = FastGibbsSamplers:::effective_sample_size(sigma_samples)
    Ef_sigma2 = ESS_sigma2/time_val
    ESS_lambda2  = FastGibbsSamplers:::effective_sample_size(lambda_samples)
    Ef_lambda2 = ESS_lambda2/time_val

    stats = c(
      100*median(vESS)/N,
      Ef,
      100*ESS_sigma2/N,2,
      Ef_sigma2,
      100*ESS_lambda2/N,
      Ef_lambda2,
      time_val)


  return(stats)
}

