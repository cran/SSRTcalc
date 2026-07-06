# =============================================================================
# R/stan_extensions.R  --  Bayesian SSRT estimation via Stan
#
# Requires cmdstanr (recommended) or rstan.
# Stan model files live in inst/stan/ and are located at runtime via
# system.file("stan", "*.stan", package = "SSRTcalc").
#
# ssrt_stan()               -- fit single or hierarchical model
# print / plot / summary /  -- S3 methods
# coef / ranef
# ssrt_stan_pp_check()      -- posterior predictive checks
# ssrt_stan_inhibition_fn() -- posterior inhibition function
# ssrt_stan_compare()       -- compare base vs. trigger-failure model
# ssrt_stan_loo()           -- LOO-CV (requires log_lik in model)
# =============================================================================


# -- Backend detection --------------------------------------------------------

.detect_stan_backend <- function() {
  if (requireNamespace("cmdstanr", quietly = TRUE)) return("cmdstanr")
  if (requireNamespace("rstan",    quietly = TRUE)) return("rstan")
  stop(
    "No Stan backend found.\n\n",
    "Install cmdstanr (recommended):\n",
    "  install.packages('cmdstanr',\n",
    "    repos = c('https://mc-stan.org/r-packages/', getOption('repos')))\n",
    "  cmdstanr::install_cmdstan()\n\n",
    "Or rstan:\n",
    "  install.packages('rstan')"
  )
}


# -- Locate a bundled Stan model file ----------------------------------------

.stan_model_path <- function(name) {
  path <- system.file("stan", name, package = "SSRTcalc")
  if (nchar(path) == 0)
    stop("Stan model '", name, "' not found in package installation. ",
         "Try reinstalling SSRTcalc.")
  path
}


# -- Data preparation ---------------------------------------------------------

.prepare_single <- function(data, stop_col, rt_col, ssd_col, n_quad, use_tf) {
  s <- .split_trials(data, stop_col, rt_col, ssd_col)
  if (nrow(s$sr)  == 0) warning("No signal-respond trials -- stop params weakly identified.")
  if (nrow(s$inh) == 0) warning("No inhibited trials -- stop params weakly identified.")
  list(
    N_go    = nrow(s$go),
    go_rt   = s$go[[rt_col]],
    N_sr    = nrow(s$sr),
    sr_rt   = if (nrow(s$sr)  > 0) s$sr[[rt_col]]  else numeric(0),
    sr_ssd  = if (nrow(s$sr)  > 0) s$sr[[ssd_col]] else numeric(0),
    N_inh   = nrow(s$inh),
    inh_ssd = if (nrow(s$inh) > 0) s$inh[[ssd_col]] else numeric(0),
    n_quad  = as.integer(n_quad),
    use_tf  = as.integer(use_tf)
  )
}

.prepare_hier <- function(data, subject_col, stop_col, rt_col,
                           ssd_col, n_quad, use_tf) {
  ids     <- sort(unique(data[[subject_col]]))
  J       <- length(ids)
  map     <- setNames(seq_along(ids), as.character(ids))
  .si     <- function(df) as.integer(map[as.character(df[[subject_col]])])
  s       <- .split_trials(data, stop_col, rt_col, ssd_col)
  list(
    J         = J,
    N_go      = nrow(s$go),
    go_subj   = .si(s$go),
    go_rt     = s$go[[rt_col]],
    N_sr      = nrow(s$sr),
    sr_subj   = if (nrow(s$sr)  > 0) .si(s$sr)           else integer(0),
    sr_rt     = if (nrow(s$sr)  > 0) s$sr[[rt_col]]       else numeric(0),
    sr_ssd    = if (nrow(s$sr)  > 0) s$sr[[ssd_col]]      else numeric(0),
    N_inh     = nrow(s$inh),
    inh_subj  = if (nrow(s$inh) > 0) .si(s$inh)           else integer(0),
    inh_ssd   = if (nrow(s$inh) > 0) s$inh[[ssd_col]]     else numeric(0),
    n_quad    = as.integer(n_quad),
    use_tf    = as.integer(use_tf),
    .ids      = ids
  )
}


# -- Data-driven initial values ----------------------------------------------
# The hierarchical model carries population parameters on the log scale
# (mu_go_pop_log, ...). Stan's default initialisation draws unconstrained
# parameters uniformly from (-2, 2); for a log-mean that should sit near
# log(400) ~ 6 this starts mu_go at exp(0) = 1 ms, making every go-trial
# density log(0) (chains rejected at init) and, once the sd * z term is added,
# overflowing exp() to Inf ("Location parameter is inf" / "Inv_scale parameter
# is 0"). Supplying sensible starting values derived from the data puts every
# chain in the right region and removes the failure. Values are in milliseconds
# to match the rest of the package.
.make_init <- function(model, stan_data, trigger_failure, chains, seed) {
  set.seed(seed)
  go <- stan_data$go_rt
  go <- go[is.finite(go) & go > 0]
  m  <- mean(go)
  s  <- if (length(go) > 1) sd(go) else m * 0.3

  # Rough, robust starting values (ms), clamped well inside the hard parameter
  # bounds declared in the Stan models (with margin to absorb the per-chain
  # jitter applied below) so the backend never rejects a supplied init.
  clamp     <- function(x, lo, hi) min(max(x, lo), hi)
  mu_go0    <- clamp(m - s,     55,  1380)   # Gaussian component of the go RT
  sigma_go0 <- clamp(s * 0.5,   3.3,  460)
  tau_go0   <- clamp(s * 0.8,   3.3,  740)    # exponential tail of the go RT
  mu_st0    <- clamp(0.45 * m,  22,   740)    # SSRT ~ 45% of mean go RT
  sigma_st0 <- clamp(s * 0.25,  2.2,  370)
  tau_st0   <- clamp(s * 0.30,  2.2,  460)

  jit <- function(x, lo = 0.95, hi = 1.05) x * stats::runif(1, lo, hi)

  one_single <- function() {
    init <- list(
      mu_go    = jit(mu_go0),  sigma_go   = jit(sigma_go0), tau_go   = jit(tau_go0),
      mu_stop  = jit(mu_st0),  sigma_stop = jit(sigma_st0), tau_stop = jit(tau_st0)
    )
    if (trigger_failure) init$p_tf <- as.array(0.05)
    init
  }

  one_hier <- function() {
    J <- stan_data$J
    init <- list(
      mu_go_pop_log      = log(jit(mu_go0, 0.97, 1.03)),
      sigma_go_pop_log   = log(jit(sigma_go0, 0.97, 1.03)),
      tau_go_pop_log     = log(jit(tau_go0, 0.97, 1.03)),
      mu_stop_pop_log    = log(jit(mu_st0, 0.97, 1.03)),
      sigma_stop_pop_log = log(jit(sigma_st0, 0.97, 1.03)),
      tau_stop_pop_log   = log(jit(tau_st0, 0.97, 1.03)),
      sd_mu_go = 0.1, sd_sigma_go = 0.1, sd_tau_go = 0.1,
      sd_mu_stop = 0.1, sd_sigma_stop = 0.1, sd_tau_stop = 0.1,
      z_mu_go    = stats::rnorm(J, 0, 0.1),
      z_sigma_go = stats::rnorm(J, 0, 0.1),
      z_tau_go   = stats::rnorm(J, 0, 0.1),
      z_mu_stop    = stats::rnorm(J, 0, 0.1),
      z_sigma_stop = stats::rnorm(J, 0, 0.1),
      z_tau_stop   = stats::rnorm(J, 0, 0.1)
    )
    if (trigger_failure) init$p_tf <- rep(0.05, J)
    init
  }

  gen <- if (model == "single") one_single else one_hier
  lapply(seq_len(chains), function(i) gen())
}


# -- Fit wrapper --------------------------------------------------------------

.fit_stan <- function(stan_file, stan_data, backend,
                       chains, iter, warmup, cores, adapt_delta, seed,
                       init = NULL, ...) {
  if (backend == "cmdstanr") {
    # Compile from a temporary copy of the .stan file so that the compiled
    # executable (and any intermediate .hpp) is written to a writable temp
    # directory, never into the installed package's inst/stan directory.
    # Writing there would leave an undeclared executable in the package (an
    # R CMD check WARNING) and can pollute source/tarball rebuilds. Within a
    # session the temp copy persists, so the model is compiled once and reused.
    tmp_stan <- file.path(tempdir(), basename(stan_file))
    if (!file.exists(tmp_stan) ||
        file.mtime(tmp_stan) < file.mtime(stan_file)) {
      file.copy(stan_file, tmp_stan, overwrite = TRUE)
    }
    mod  <- cmdstanr::cmdstan_model(tmp_stan)
    args <- list(
      data            = stan_data,
      chains          = chains,
      iter_warmup     = warmup,
      iter_sampling   = iter - warmup,
      parallel_chains = min(cores, chains),
      seed            = seed,
      adapt_delta     = adapt_delta,
      refresh         = max(1L, as.integer((iter - warmup) / 10L)),
      ...
    )
    if (!is.null(init)) args$init <- init
    do.call(mod$sample, args)
  } else {
    args <- list(
      file    = stan_file,
      data    = stan_data,
      chains  = chains,
      iter    = iter,
      warmup  = warmup,
      cores   = min(cores, chains),
      seed    = seed,
      control = list(adapt_delta = adapt_delta),
      ...
    )
    if (!is.null(init)) args$init <- init
    do.call(rstan::stan, args)
  }
}


# =============================================================================
# ssrt_stan()  --  main fitting function
# =============================================================================

#' Bayesian SSRT estimation via Stan
#'
#' Fits the independent horse-race model with ex-Gaussian go- and
#' stop-process distributions using Hamiltonian Monte Carlo via Stan.
#' Supports single-subject and multi-subject hierarchical designs, with an
#' optional trigger-failure parameter following Matzke et al. (2013).
#'
#' @param data            data.frame in SSRTcalc long format.
#' @param model           \code{"single"} (default) or \code{"hierarchical"}.
#' @param subject_col     Column with subject IDs (required for hierarchical).
#' @param trigger_failure Logical: include trigger-failure parameter. Default FALSE.
#' @param n_quad          Quadrature resolution for P(inhibit) integral. Default 100.
#' @param chains          Number of MCMC chains. Default 4.
#' @param iter            Total iterations per chain including warmup. Default 2000.
#' @param warmup          Warmup iterations. Default 1000.
#' @param cores           Parallel cores. Default 4.
#' @param backend         \code{"cmdstanr"}, \code{"rstan"}, or \code{"auto"}.
#' @param adapt_delta     HMC target acceptance rate. Increase to 0.99 if
#'   divergences appear. Default 0.95.
#' @param stop_col,rt_col,acc_col,ssd_col Column names.
#' @param seed            Random seed. Default 42.
#' @param ...             Additional arguments passed to the Stan sampler.
#'
#' @return An object of class \code{ssrt_stan}.
#'
#' @references
#' Matzke, D., et al. (2013). Bayesian parametric estimation of stop-signal
#' reaction time distributions. \emph{Journal of Experimental Psychology:
#' General}, 142(4), 1047--1073.
#'
#' @examples
#' \dontrun{
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' fit <- ssrt_stan(d, chains = 4, iter = 2000)
#' print(fit)
#' plot(fit)
#'
#' # Hierarchical
#' fit_h <- ssrt_stan(adaptive, model = "hierarchical",
#'                    subject_col = "SubjID", chains = 4)
#' ranef(fit_h)
#' }
#'
#' @export
ssrt_stan <- function(data,
                      model           = c("single", "hierarchical"),
                      subject_col     = NULL,
                      trigger_failure = FALSE,
                      n_quad          = 100,
                      chains          = 4,
                      iter            = 2000,
                      warmup          = 1000,
                      cores           = 4,
                      backend         = "auto",
                      adapt_delta     = 0.95,
                      stop_col        = "vol",
                      rt_col          = "RT_exp",
                      acc_col         = "correct",
                      ssd_col         = "soa",
                      seed            = 42,
                      ...) {
  model   <- match.arg(model)
  backend <- if (identical(backend, "auto")) .detect_stan_backend() else backend

  cat(sprintf("SSRTcalc Stan  |  model: %s  |  backend: %s  |  TF: %s\n",
              model, backend, trigger_failure))

  if (model == "single") {
    stan_data  <- .prepare_single(data, stop_col, rt_col, ssd_col, n_quad, trigger_failure)
    stan_file  <- .stan_model_path("ssrt_single.stan")
    subj_ids   <- NULL
    key_params <- c("mu_go","sigma_go","tau_go","mu_stop","sigma_stop","tau_stop",
                    "mean_ssrt","mean_go_rt")
    if (trigger_failure) key_params <- c(key_params, "p_tf[1]")
  } else {
    if (is.null(subject_col))
      stop("`subject_col` must be specified for hierarchical model.")
    sd       <- .prepare_hier(data, subject_col, stop_col, rt_col, ssd_col,
                               n_quad, trigger_failure)
    subj_ids <- sd$.ids;  sd$.ids <- NULL
    stan_data <- sd
    stan_file  <- .stan_model_path("ssrt_hier.stan")
    key_params <- c("pop_mean_ssrt",
                    "mu_go_pop_log","sigma_go_pop_log","tau_go_pop_log",
                    "mu_stop_pop_log","sigma_stop_pop_log","tau_stop_pop_log",
                    "sd_mu_go","sd_sigma_go","sd_tau_go",
                    "sd_mu_stop","sd_sigma_stop","sd_tau_stop")
  }

  cat(sprintf("  Data: %d go | %d signal-respond | %d inhibited\n",
              stan_data$N_go, stan_data$N_sr, stan_data$N_inh))
  cat(sprintf("  %d chains x %d iter (%d warmup)  adapt_delta=%.3f\n\n",
              chains, iter, warmup, adapt_delta))

  dots      <- list(...)
  init_vals <- if ("init" %in% names(dots)) dots[["init"]]
               else .make_init(model, stan_data, trigger_failure, chains, seed)
  dots[["init"]] <- NULL

  fit <- do.call(.fit_stan,
                 c(list(stan_file, stan_data, backend, chains, iter, warmup,
                        cores, adapt_delta, seed, init = init_vals), dots))

  summary_df <- .stan_summary_df(fit, backend, key_params)

  structure(
    list(fit=fit, summary=summary_df, backend=backend, model_type=model,
         stan_data=stan_data, subject_ids=subj_ids,
         trigger_failure=trigger_failure, n_quad=n_quad,
         key_params=key_params, call=match.call()),
    class = c("ssrt_stan","list")
  )
}


# =============================================================================
# S3 methods
# =============================================================================

#' @export
print.ssrt_stan <- function(x, ...) {
  cat("== SSRTcalc Stan Fit ===========================================\n")
  cat(sprintf("  Model: %s  |  Backend: %s  |  TF: %s\n",
              x$model_type, x$backend, x$trigger_failure))
  cat(sprintf("  Data:  %d go | %d SR | %d inhibited\n\n",
              x$stan_data$N_go, x$stan_data$N_sr, x$stan_data$N_inh))

  s <- x$summary
  if (x$model_type == "single") {
    show <- c("mu_go","sigma_go","tau_go","mu_stop","sigma_stop","tau_stop","mean_ssrt")
    if (x$trigger_failure) show <- c(show, "p_tf[1]")
  } else {
    show <- c("pop_mean_ssrt","sd_mu_go","sd_sigma_go","sd_tau_go",
              "sd_mu_stop","sd_sigma_stop","sd_tau_stop")
  }

  cat(sprintf("  %-22s  %8s  %6s  %18s  %6s\n","Parameter","Mean","SD","90% CI","Rhat"))
  cat(sprintf("  %s\n", paste(rep("-",62),collapse="")))
  for (p in show) {
    r <- s[s$variable == p, ]
    if (nrow(r) == 0) next
    cat(sprintf("  %-22s  %8.3f  %6.3f  [%7.2f, %7.2f]  %6.3f\n",
                p, r$mean, r$sd, r$q5, r$q95, r$rhat))
  }

  max_rhat <- max(s$rhat,     na.rm=TRUE)
  min_ess  <- min(s$ess_bulk, na.rm=TRUE)
  cat(sprintf("\n  Max Rhat: %.4f %s\n", max_rhat,
              if (max_rhat < 1.01) "OK" else "WARN: >1.01"))
  cat(sprintf("  Min ESS:  %.0f %s\n", min_ess,
              if (min_ess > 400) "OK" else "WARN: <400"))

  if (x$backend == "cmdstanr") {
    divs <- tryCatch(sum(x$fit$diagnostic_summary()$num_divergent), error=function(e) NA)
    if (!is.na(divs))
      cat(sprintf("  Divergences: %d %s\n", divs,
                  if (divs==0) "OK" else "WARN: try adapt_delta=0.99"))
  }
  invisible(x)
}

#' @export
summary.ssrt_stan <- function(object, ...) {
  cat("Full posterior summary:\n\n")
  print(object$summary, digits=4, row.names=FALSE)
  invisible(object$summary)
}

#' @export
coef.ssrt_stan <- function(object, ...) {
  s <- object$summary
  setNames(s$mean, s$variable)
}

#' @export
plot.ssrt_stan <- function(x, type=c("posterior","trace","pairs"), ...) {
  type <- match.arg(type)

  stop_pars <- if (x$model_type=="single")
    c("mu_stop","sigma_stop","tau_stop","mean_ssrt")
  else
    c("pop_mean_ssrt","sd_mu_stop","sd_sigma_stop","sd_tau_stop")

  get_draws_arr <- function() {
    if (x$backend=="cmdstanr") x$fit$draws(format="array")
    else as.array(x$fit)
  }

  if (type == "posterior") {
    if (x$model_type == "single") {
      ssrt_d <- .extract_draws(x$fit, x$backend, "mean_ssrt")
      go_d   <- .extract_draws(x$fit, x$backend, "mean_go_rt")
      op <- par(mfrow=c(1,2), mar=c(4,4,3,1));  on.exit(par(op))
      for (draws in list(ssrt_d, go_d)) {
        lbl <- if (identical(draws, ssrt_d)) "Mean SSRT" else "Mean Go RT"
        hist(draws, breaks=50, col=if(identical(draws,ssrt_d)) "#59A14F" else "#4E79A7",
             border="white", freq=FALSE, main=paste("Posterior:", lbl), xlab="ms")
        lines(density(draws), lwd=2)
        abline(v=quantile(draws, c(0.05,0.5,0.95)),
               col=c("darkorange","firebrick","darkorange"), lty=c(2,1,2), lwd=2)
      }
    } else {
      # Hierarchical: forest plot
      ssrt_pars <- paste0("mean_ssrt[", seq_along(x$subject_ids), "]")
      ss <- .stan_summary_df(x$fit, x$backend, ssrt_pars)
      J  <- nrow(ss)
      op <- par(mar=c(4,7,3,2));  on.exit(par(op))
      plot(ss$mean, J:1, xlim=range(c(ss$q5, ss$q95)),
           xlab="Mean SSRT (ms)", ylab="", main="Per-Subject SSRT Posteriors",
           yaxt="n", pch=19, col="#59A14F")
      segments(ss$q5, J:1, ss$q95, J:1, col="#59A14F", lwd=2)
      axis(2, at=J:1, labels=x$subject_ids, las=1, cex.axis=0.8)
      abline(v=mean(ss$mean), lty=2, col="firebrick")
    }

  } else if (type == "trace") {
    if (requireNamespace("bayesplot", quietly=TRUE)) {
      print(bayesplot::mcmc_trace(get_draws_arr(), pars=stop_pars))
    } else if (x$backend=="rstan") {
      rstan::traceplot(x$fit, pars=stop_pars)
    } else {
      message("Install bayesplot for trace plots: install.packages('bayesplot')")
    }

  } else if (type == "pairs") {
    if (requireNamespace("bayesplot", quietly=TRUE)) {
      print(bayesplot::mcmc_pairs(get_draws_arr(), pars=stop_pars))
    } else {
      message("Install bayesplot for pairs plots: install.packages('bayesplot')")
    }
  }
  invisible(x)
}

#' Extract per-subject random effects from a hierarchical fit
#'
#' @param object An \code{ssrt_stan} object from a hierarchical model.
#' @param ... Unused.
#' @return data.frame with per-subject posterior summaries.
#' @export
ranef.ssrt_stan <- function(object, ...) {
  if (object$model_type != "hierarchical")
    stop("ranef() is only applicable to hierarchical fits.")
  J    <- object$stan_data$J
  pars <- c(paste0("mean_ssrt[", 1:J, "]"),
            paste0("mu_stop[",   1:J, "]"),
            paste0("tau_stop[",  1:J, "]"))
  df          <- .stan_summary_df(object$fit, object$backend, pars)
  df$subject  <- rep(object$subject_ids, 3)
  df[, c("subject", setdiff(names(df), "subject"))]
}


# =============================================================================
# Diagnostic and analysis functions
# =============================================================================

#' Posterior predictive checks for an ssrt_stan fit
#'
#' Overlays posterior-predictive RT distributions over observed data.
#'
#' @param fit       An \code{ssrt_stan} object (single-subject only).
#' @param n_samples Posterior draws to use. Default 200.
#'
#' @return Invisibly, a list with \code{pp_go_rt} and \code{pp_ssrt}.
#'
#' @examples
#' \dontrun{
#' data(adaptive)
#' fit <- ssrt_stan(adaptive[adaptive$SubjID == 1, ])
#' ssrt_stan_pp_check(fit)
#' }
#'
#' @export
ssrt_stan_pp_check <- function(fit, n_samples = 200) {
  if (fit$model_type != "single")
    stop("PP check currently supported for single-subject model only.")

  mu_go_d    <- .extract_draws(fit$fit, fit$backend, "mu_go")
  sigma_go_d <- .extract_draws(fit$fit, fit$backend, "sigma_go")
  tau_go_d   <- .extract_draws(fit$fit, fit$backend, "tau_go")
  ssrt_d     <- .extract_draws(fit$fit, fit$backend, "ssrt_pred")

  idx    <- sample(length(mu_go_d), min(n_samples, length(mu_go_d)))
  pp_go  <- mapply(function(mu, sg, tau) rnorm(1, mu, sg) + rexp(1, 1/tau),
                   mu_go_d[idx], sigma_go_d[idx], tau_go_d[idx])
  obs_go <- fit$stan_data$go_rt

  op <- par(mfrow=c(1,2), mar=c(4,4,3,1));  on.exit(par(op))
  xlim <- range(c(obs_go, pp_go))
  d_obs <- density(obs_go, from=xlim[1], to=xlim[2])
  d_pp  <- density(pp_go,  from=xlim[1], to=xlim[2])
  plot(d_obs, col="black", lwd=2,
       ylim=c(0, max(d_obs$y, d_pp$y)*1.15),
       main="PP Check: Go RT", xlab="RT (ms)")
  lines(d_pp, col="#4E79A7", lwd=2, lty=2)
  legend("topright", legend=c("Observed","PP median"),
         col=c("black","#4E79A7"), lwd=2, lty=c(1,2), bty="n")

  hist(ssrt_d[idx], breaks=40, col="#59A14F", border="white", freq=FALSE,
       main="Posterior SSRT", xlab="SSRT (ms)")
  lines(density(ssrt_d[idx]), lwd=2)
  abline(v=quantile(ssrt_d[idx], c(0.05,0.5,0.95)),
         col=c("darkorange","firebrick","darkorange"), lty=c(2,1,2), lwd=2)

  invisible(list(pp_go_rt=pp_go, pp_ssrt=ssrt_d[idx]))
}


#' Posterior inhibition function
#'
#' Computes P(inhibit | SSD) across a range of SSD values using posterior
#' parameter samples, with a 90% credible interval band.
#'
#' @param fit       An \code{ssrt_stan} object (single-subject only).
#' @param ssd_range SSD values to evaluate in ms. Default \code{seq(0, 600, 10)}.
#' @param n_draws   Posterior draws to use. Default 400.
#' @param plot      Produce the plot? Default TRUE.
#'
#' @return Invisibly, a data.frame with columns: ssd, mean, lo90, hi90.
#'
#' @export
ssrt_stan_inhibition_fn <- function(fit,
                                    ssd_range = seq(0, 600, by = 10),
                                    n_draws   = 400,
                                    plot      = TRUE) {
  if (fit$model_type != "single")
    stop("inhibition_fn currently supported for single-subject model only.")

  mu_go_d    <- .extract_draws(fit$fit, fit$backend, "mu_go")
  sigma_go_d <- .extract_draws(fit$fit, fit$backend, "sigma_go")
  tau_go_d   <- .extract_draws(fit$fit, fit$backend, "tau_go")
  mu_sp_d    <- .extract_draws(fit$fit, fit$backend, "mu_stop")
  sigma_sp_d <- .extract_draws(fit$fit, fit$backend, "sigma_stop")
  tau_sp_d   <- .extract_draws(fit$fit, fit$backend, "tau_stop")

  n   <- min(n_draws, length(mu_go_d))
  idx <- sample(length(mu_go_d), n)
  cat("Computing inhibition function across", n, "posterior draws...\n")

  p_mat <- vapply(idx, function(i)
    vapply(ssd_range, function(ssd)
      .p_inhibit_r(ssd, mu_go_d[i], sigma_go_d[i], tau_go_d[i],
                   mu_sp_d[i], sigma_sp_d[i], tau_sp_d[i]),
      numeric(1)),
    numeric(length(ssd_range)))

  df <- data.frame(
    ssd  = ssd_range,
    mean = rowMeans(p_mat),
    lo90 = apply(p_mat, 1, quantile, 0.05),
    hi90 = apply(p_mat, 1, quantile, 0.95)
  )

  if (plot) {
    plot(df$ssd, df$mean, type="l", lwd=2, col="#E15759",
         ylim=c(0,1), xlab="Stop-Signal Delay (ms)", ylab="P(inhibit | SSD)",
         main="Posterior Inhibition Function")
    polygon(c(df$ssd, rev(df$ssd)), c(df$hi90, rev(df$lo90)),
            col=adjustcolor("#E15759",0.2), border=NA)
    abline(h=0.5, lty=2, col="grey50")
    legend("topright", legend=c("Posterior mean","90% CI"),
           col=c("#E15759",adjustcolor("#E15759",0.4)), lwd=c(2,8), bty="n")
  }
  invisible(df)
}


#' Compare base model vs. trigger-failure model
#'
#' Fits both models and returns side-by-side SSRT posterior summaries.
#'
#' @param data  data.frame in SSRTcalc long format.
#' @param ...   Arguments forwarded to \code{ssrt_stan()}.
#'
#' @return Invisibly, a list with \code{base}, \code{tf}, and \code{comparison}.
#'
#' @export
ssrt_stan_compare <- function(data, ...) {
  cat("=== Model comparison: base vs. trigger-failure ===\n\n")
  cat("Step 1/2: Base model (no trigger failures)...\n")
  fit_base <- ssrt_stan(data, trigger_failure=FALSE, ...)
  cat("\nStep 2/2: Trigger-failure model...\n")
  fit_tf   <- ssrt_stan(data, trigger_failure=TRUE,  ...)

  ssrt_b  <- .extract_draws(fit_base$fit, fit_base$backend, "mean_ssrt")
  ssrt_t  <- .extract_draws(fit_tf$fit,   fit_tf$backend,   "mean_ssrt")
  p_tf_d  <- .extract_draws(fit_tf$fit,   fit_tf$backend,   "p_tf[1]")

  comp <- data.frame(
    model     = c("Base (no TF)","Trigger-failure"),
    mean_ssrt = c(mean(ssrt_b), mean(ssrt_t)),
    sd_ssrt   = c(sd(ssrt_b),   sd(ssrt_t)),
    lo90      = c(quantile(ssrt_b,.05), quantile(ssrt_t,.05)),
    hi90      = c(quantile(ssrt_b,.95), quantile(ssrt_t,.95)),
    mean_p_tf = c(NA, mean(p_tf_d)),
    max_rhat  = c(max(fit_base$summary$rhat, na.rm=TRUE),
                  max(fit_tf$summary$rhat, na.rm=TRUE))
  )
  cat("\n-- Comparison Summary --\n")
  print(comp, digits=3, row.names=FALSE)
  invisible(list(base=fit_base, tf=fit_tf, comparison=comp))
}


#' Leave-one-out cross-validation for an ssrt_stan fit
#'
#' Requires the Stan model to include a \code{log_lik} vector in
#' generated quantities. If absent, an informative error explains how to add it.
#'
#' @param fit An \code{ssrt_stan} object.
#'
#' @return A \code{loo} object.
#'
#' @export
ssrt_stan_loo <- function(fit) {
  if (!requireNamespace("loo", quietly=TRUE))
    stop("loo package required: install.packages('loo')")

  ll <- tryCatch({
    if (fit$backend == "cmdstanr")
      fit$fit$draws(variables="log_lik", format="matrix")
    else {
      # Resolved dynamically (rather than via rstan::extract_log_lik) so
      # that R CMD check does not require rstan's namespace to verify the
      # symbol when rstan is an unmet Suggests dependency.
      extract_log_lik_fn <- getExportedValue("rstan", "extract_log_lik")
      extract_log_lik_fn(fit$fit, parameter_name="log_lik", merge_chains=FALSE)
    }
  }, error=function(e) NULL)

  if (is.null(ll)) stop(
    "`log_lik` not found in the fitted model.\n\n",
    "To enable LOO, add to generated quantities in the Stan model:\n\n",
    "  vector[N_go + N_sr + N_inh] log_lik;\n",
    "  for (i in 1:N_go)\n",
    "    log_lik[i] = exp_mod_normal_lpdf(go_rt[i] | mu_go, sigma_go, lambda_go);\n",
    "  for (i in 1:N_sr)\n",
    "    log_lik[N_go+i] = exp_mod_normal_lpdf(sr_rt[i] | mu_go, sigma_go, lambda_go)\n",
    "                    + exp_mod_normal_lccdf(sr_rt[i]-sr_ssd[i] | mu_stop, sigma_stop, lambda_stop);\n",
    "  for (i in 1:N_inh)\n",
    "    log_lik[N_go+N_sr+i] = log_p_inhibit(inh_ssd[i],\n",
    "                             mu_go, sigma_go, lambda_go,\n",
    "                             mu_stop, sigma_stop, lambda_stop, n_quad);\n"
  )

  r_eff   <- loo::relative_eff(exp(ll))
  loo_obj <- loo::loo(ll, r_eff=r_eff)
  cat("-- LOO-CV results --\n");  print(loo_obj)
  invisible(loo_obj)
}
