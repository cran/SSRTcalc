# =============================================================================
# R/helpers.R
# Internal helpers shared by Monte Carlo and Stan extensions.
# None of these are exported.
# =============================================================================


# -- Ex-Gaussian parameter estimation ----------------------------------------

#' Fit ex-Gaussian parameters via maximum likelihood (Nelder-Mead),
#' starting from method-of-moments estimates when they are valid.
#'
#' Method-of-moments alone (tau = var/mean, mu = mean - tau,
#' sigma = sqrt(var - tau^2)) is badly biased whenever tau is not small
#' relative to sigma (it can underestimate tau by a factor of several and
#' correspondingly overestimate sigma and mu). MLE refinement is therefore
#' always performed; MoM merely supplies good starting values when valid.
#'
#' @noRd
.fit_exgaussian <- function(rt) {
  rt <- rt[!is.na(rt) & rt > 0]
  if (length(rt) < 5)
    stop("Need >= 5 valid RTs to fit ex-Gaussian.")

  m <- mean(rt);  v <- var(rt)

  # Method-of-moments starting values
  tau0   <- v / m
  mu0    <- m - tau0
  sigma0 <- sqrt(max(1e-6, v - tau0^2))

  # If MoM produced invalid (non-positive) values, use generic
  # rule-of-thumb starting values instead.
  if (!(mu0 > 0 && sigma0 > 0 && tau0 > 0)) {
    mu0    <- m * 0.8
    sigma0 <- sd(rt) * 0.5
    tau0   <- m * 0.2
  }

  neg_ll <- function(par) {
    mu_ <- par[1];  sigma_ <- exp(par[2]);  tau_ <- exp(par[3])
    lam  <- 1 / tau_

    # Outside this region the closed-form ex-Gaussian log-density suffers
    # catastrophic cancellation between two O(exp((sigma*lam)^2)) terms of
    # opposite sign, which can make wildly divergent (sigma >> tau)
    # parameter values look spuriously attractive to the optimiser. No RT
    # data plausibly needs sigma/tau outside this range, so it is excluded
    # outright rather than relying on is.finite() downstream.
    if (sigma_ * lam > 50 ||
        sigma_ < 1e-3 || sigma_ > 1e5 ||
        tau_   < 1e-3 || tau_   > 1e5) {
      return(1e10)
    }

    val <- -sum(log(lam) + lam * (mu_ + lam * sigma_^2 / 2 - rt) +
                   pnorm((rt - mu_ - lam * sigma_^2) / sigma_, log.p = TRUE))
    if (!is.finite(val)) return(1e10)
    val
  }

  fit <- tryCatch(
    optim(c(mu0, log(sigma0), log(tau0)),
          neg_ll, method = "Nelder-Mead",
          control = list(maxit = 3000)),
    error = function(e) NULL
  )

  if (!is.null(fit) && fit$convergence == 0) {
    mu    <- fit$par[1]
    sigma <- exp(fit$par[2])
    tau   <- exp(fit$par[3])
  } else {
    # MLE failed to converge -- fall back to the starting values
    mu <- mu0; sigma <- sigma0; tau <- tau0
  }

  list(mu = mu, sigma = sigma, tau = tau)
}


#' Sample n draws from ExGaussian(mu, sigma, tau)
#' @noRd
.rexgaussian <- function(n, mu, sigma, tau) {
  rnorm(n, mean = mu, sd = sigma) + rexp(n, rate = 1 / tau)
}


#' Ex-Gaussian CDF (consistent with Stan's exp_mod_normal)
#' F(x; mu, sigma, tau):  Phi((x-mu)/sigma)
#'   - exp(lambda*(mu + lambda*sigma^2/2 - x)) * Phi((x-mu-lambda*sigma^2)/sigma)
#' where lambda = 1/tau
#' @noRd
.pexgaussian <- function(x, mu, sigma, tau) {
  lam <- 1 / tau
  pnorm((x - mu) / sigma) -
    exp(lam * (mu + lam * sigma^2 / 2 - x)) *
    pnorm((x - mu - lam * sigma^2) / sigma)
}


#' Multivariate normal sampler (base-R, avoids MASS dependency by default)
#' Falls back to MASS::mvrnorm if available.
#' @noRd
.mvrnorm_simple <- function(n, mu, Sigma) {
  if (requireNamespace("MASS", quietly = TRUE)) {
    return(MASS::mvrnorm(n, mu, Sigma))
  }
  p <- length(mu)
  L <- chol(Sigma)
  matrix(rnorm(n * p), nrow = n, ncol = p) %*% L +
    matrix(mu, nrow = n, ncol = p, byrow = TRUE)
}


# -- Stop-trial splitter (shared by Stan helpers) -----------------------------

#' Split a trial data.frame into go / signal-respond / inhibited sub-frames
#' @noRd
.split_trials <- function(data, stop_col, rt_col, ssd_col) {
  go_trials  <- data[data[[stop_col]] == 0, ]
  stp_trials <- data[data[[stop_col]] == 1, ]

  sr_mask  <- !is.na(stp_trials[[rt_col]]) & stp_trials[[rt_col]] > 0
  list(
    go  = go_trials[ !is.na(go_trials[[rt_col]]) & go_trials[[rt_col]] > 0, ],
    sr  = stp_trials[ sr_mask, ],
    inh = stp_trials[!sr_mask, ]
  )
}


# -- Numerical P(inhibit | SSD) in R (mirrors Stan quadrature) ----------------

#' Compute P(inhibit | SSD) using midpoint quadrature (R implementation)
#' Used by ssrt_stan_inhibition_fn() and ssrt_robustness().
#' @noRd
.p_inhibit_r <- function(ssd, mu_go, sigma_go, tau_go,
                          mu_sp, sigma_sp, tau_sp,
                          n_quad = 100) {
  upper <- mu_sp + tau_sp + 5 * sqrt(sigma_sp^2 + tau_sp^2)
  dt    <- upper / n_quad
  t_seq <- seq(dt / 2, upper - dt / 2, by = dt)

  f_stop <- diff(c(0, .pexgaussian(t_seq, mu_sp, sigma_sp, tau_sp)))
  f_stop <- pmax(f_stop, 0)
  s_go   <- 1 - .pexgaussian(ssd + t_seq, mu_go, sigma_go, tau_go)
  s_go   <- pmax(pmin(s_go, 1), 0)

  sum(f_stop * s_go)
}


# -- Generic extractor of posterior draws -------------------------------------

#' Extract a flat vector of posterior draws for one parameter
#' Works with both cmdstanr and rstan objects.
#' @noRd
.extract_draws <- function(fit, backend, param) {
  if (backend == "cmdstanr") {
    m <- fit$draws(variables = param, format = "matrix")
    as.vector(m)
  } else {
    as.vector(rstan::extract(fit, pars = param)[[param]])
  }
}


#' Extract posterior summary as tidy data.frame
#' @noRd
.stan_summary_df <- function(fit, backend, params) {
  if (backend == "cmdstanr") {
    s <- fit$summary(variables = params)
    cols <- intersect(c("variable","mean","sd","q5","median","q95","rhat","ess_bulk"),
                      names(s))
    as.data.frame(s[, cols])
  } else {
    sm <- rstan::summary(fit, pars = params)$summary
    data.frame(
      variable = rownames(sm),
      mean     = sm[, "mean"],
      sd       = sm[, "sd"],
      q5       = sm[, "2.5%"],
      median   = sm[, "50%"],
      q95      = sm[, "97.5%"],
      rhat     = sm[, "Rhat"],
      ess_bulk = sm[, "n_eff"],
      row.names = NULL
    )
  }
}
