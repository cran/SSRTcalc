# =============================================================================
# R/mc_extensions.R  --  Monte Carlo extensions for SSRT estimation
# ssrt_boot()        -- nonparametric bootstrap confidence intervals
# ssrt_simulate()    -- parametric ex-Gaussian Monte Carlo
# ssrt_power()       -- minimum trial count / power analysis
# ssrt_robustness()  -- sensitivity to horse-race assumption violations
# run_all_mc()       -- convenience wrapper
# =============================================================================


# =============================================================================
# 1.  ssrt_boot()
# =============================================================================

#' Bootstrap confidence intervals for SSRT
#'
#' Resamples trials with replacement \code{n_iter} times and applies the
#' chosen SSRT estimation function to each resample.
#'
#' @param data      data.frame in SSRTcalc long format.
#' @param ssrt_fn   SSRT function name. Default \code{"integration_adaptiveSSD"}.
#' @param n_iter    Bootstrap resamples. Default 2000.
#' @param conf      Confidence level. Default 0.95.
#' @param stop_col,rt_col,acc_col,ssd_col Column names.
#' @param seed      Random seed. Default 42.
#' @param parallel  Use parallel::mclapply? (Unix/macOS only). Default FALSE.
#' @param n_cores   Cores when parallel=TRUE. Default 2.
#'
#' @return Object of class \code{ssrt_boot}.
#'
#' @examples
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' b <- ssrt_boot(d, n_iter = 500)
#' print(b)
#'
#' @export
ssrt_boot <- function(data,
                      ssrt_fn  = "integration_adaptiveSSD",
                      n_iter   = 2000,
                      conf     = 0.95,
                      stop_col = "vol",
                      rt_col   = "RT_exp",
                      acc_col  = "correct",
                      ssd_col  = "soa",
                      seed     = 42,
                      parallel = FALSE,
                      n_cores  = 2) {
  set.seed(seed)
  fn       <- match.fun(ssrt_fn)
  ssrt_obs <- fn(data, stop_col=stop_col, rt_col=rt_col,
                 acc_col=acc_col, ssd_col=ssd_col)
  n        <- nrow(data)

  one_boot <- function(i) {
    tryCatch(
      fn(data[sample(n, replace=TRUE), ],
         stop_col=stop_col, rt_col=rt_col, acc_col=acc_col, ssd_col=ssd_col),
      error = function(e) NA_real_)
  }

  if (parallel) {
    if (!requireNamespace("parallel", quietly=TRUE))
      stop("parallel package required when parallel=TRUE.")
    dist <- unlist(parallel::mclapply(seq_len(n_iter), one_boot, mc.cores=n_cores))
  } else {
    dist <- vapply(seq_len(n_iter), one_boot, numeric(1))
  }

  dist  <- dist[!is.na(dist)]
  alpha <- 1 - conf
  ci    <- quantile(dist, probs=c(alpha/2, 1-alpha/2))
  names(ci) <- c(paste0("lower_",round(conf*100),"%"),
                 paste0("upper_",round(conf*100),"%"))

  structure(list(ssrt_obs=ssrt_obs, mean=mean(dist), bias=mean(dist)-ssrt_obs,
                 se=sd(dist), ci=ci, distribution=dist, n_valid=length(dist),
                 n_iter=n_iter, conf=conf, method=ssrt_fn),
            class=c("ssrt_boot","list"))
}

#' @export
print.ssrt_boot <- function(x, ...) {
  cat("-- SSRTcalc Bootstrap Results ----------------------------------\n")
  cat(sprintf("  Method:        %s\n", x$method))
  cat(sprintf("  Iterations:    %d  (valid: %d)\n", x$n_iter, x$n_valid))
  cat(sprintf("  Observed SSRT: %7.2f ms\n", x$ssrt_obs))
  cat(sprintf("  Bootstrap M:   %7.2f ms  (bias: %+.2f ms)\n", x$mean, x$bias))
  cat(sprintf("  SE:            %7.2f ms\n", x$se))
  cat(sprintf("  %d%% CI:        [%.2f, %.2f] ms\n",
              round(x$conf*100), x$ci[1], x$ci[2]))
  invisible(x)
}

#' @export
plot.ssrt_boot <- function(x, ...) {
  hist(x$distribution, breaks=40, col="#4E79A7", border="white",
       main=paste("Bootstrap SSRT Distribution --", x$method),
       xlab="SSRT (ms)", ylab="Frequency")
  abline(v=x$ssrt_obs, col="firebrick",  lwd=2, lty=1)
  abline(v=x$ci,       col="darkorange", lwd=2, lty=2)
  legend("topright", legend=c("Observed", paste0(round(x$conf*100),"% CI")),
         col=c("firebrick","darkorange"), lwd=2, lty=c(1,2), bty="n")
}


# =============================================================================
# 2.  ssrt_simulate()
# =============================================================================

#' Parametric Monte Carlo SSRT estimation via ex-Gaussian simulation
#'
#' Fits an ex-Gaussian to observed go-RT data, simulates synthetic datasets
#' under the horse-race model, and estimates SSRT on each.
#'
#' @param data      data.frame in SSRTcalc long format.
#' @param n_iter    MC iterations. Default 2000.
#' @param n_trials  Trials per simulated dataset. NULL uses nrow(data).
#' @param p_stop    Stop-trial proportion. NULL uses observed proportion.
#' @param ssrt_true Known true SSRT for parameter recovery. Default NULL.
#' @param conf      Confidence level. Default 0.95.
#' @param stop_col,rt_col,acc_col,ssd_col Column names.
#' @param seed      Random seed. Default 42.
#'
#' @return Object of class \code{ssrt_simulate}.
#'
#' @examples
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' s <- ssrt_simulate(d, n_iter=500)
#' print(s)
#'
#' @export
ssrt_simulate <- function(data,
                          n_iter    = 2000,
                          n_trials  = NULL,
                          p_stop    = NULL,
                          ssrt_true = NULL,
                          conf      = 0.95,
                          stop_col  = "vol",
                          rt_col    = "RT_exp",
                          acc_col   = "correct",
                          ssd_col   = "soa",
                          seed      = 42) {
  set.seed(seed)

  go_rts   <- data[[rt_col]][data[[stop_col]] == 0]
  go_rts   <- go_rts[!is.na(go_rts) & go_rts > 0]
  exg      <- .fit_exgaussian(go_rts)
  cat(sprintf("Ex-Gaussian fit: mu=%.1f  sigma=%.1f  tau=%.1f ms\n",
              exg$mu, exg$sigma, exg$tau))

  ssrt_obs <- integration_adaptiveSSD(data, stop_col, rt_col, acc_col, ssd_col)
  ssrt_ref <- if (!is.null(ssrt_true)) ssrt_true else ssrt_obs
  N        <- if (is.null(n_trials)) nrow(data) else n_trials
  p_stp    <- if (is.null(p_stop)) mean(data[[stop_col]]==1) else p_stop
  mean_ssd <- mean(data[[ssd_col]][data[[stop_col]]==1], na.rm=TRUE)

  sim_one <- function() {
    n_stop  <- round(N * p_stp);  n_go <- N - n_stop
    go_all  <- .rexgaussian(N, exg$mu, exg$sigma, exg$tau)
    stp_rt  <- .rexgaussian(n_stop, ssrt_ref, exg$sigma*0.5, exg$tau*0.3)
    ssd_s   <- pmax(50, rnorm(n_stop, mean_ssd, 50))
    inh     <- go_all[seq_len(n_stop)] > (ssd_s + stp_rt)
    df <- data.frame(
      vol    = c(rep(1,n_stop), rep(0,n_go)),
      RT_exp = c(ifelse(inh, NA, go_all[seq_len(n_stop)]), go_all[(n_stop+1):N]),
      correct= c(as.integer(inh), rbinom(n_go, 1, 0.95)),
      soa    = c(ssd_s, rep(NA_real_, n_go)))
    tryCatch(
      integration_adaptiveSSD(df, "vol", "RT_exp", "correct", "soa"),
      error=function(e) NA_real_)
  }

  dist  <- vapply(seq_len(n_iter), function(i) sim_one(), numeric(1))
  dist  <- dist[!is.na(dist)]
  alpha <- 1 - conf
  ci    <- quantile(dist, probs=c(alpha/2, 1-alpha/2))
  names(ci) <- c(paste0("lower_",round(conf*100),"%"),
                 paste0("upper_",round(conf*100),"%"))

  structure(list(ssrt_obs=ssrt_obs, exg_params=exg, mean=mean(dist),
                 se=sd(dist), ci=ci, distribution=dist, n_valid=length(dist),
                 n_iter=n_iter, conf=conf,
                 recovery_bias=if(!is.null(ssrt_true)) mean(dist)-ssrt_true else NULL),
            class=c("ssrt_simulate","list"))
}

#' @export
print.ssrt_simulate <- function(x, ...) {
  cat("-- SSRTcalc Parametric MC Simulation ---------------------------\n")
  cat(sprintf("  Ex-Gaussian (go): mu=%.1f  sigma=%.1f  tau=%.1f ms\n",
              x$exg_params$mu, x$exg_params$sigma, x$exg_params$tau))
  cat(sprintf("  Iterations:       %d  (valid: %d)\n", x$n_iter, x$n_valid))
  cat(sprintf("  Observed SSRT:    %.2f ms\n", x$ssrt_obs))
  cat(sprintf("  Simulated M:      %.2f ms  (SE: %.2f ms)\n", x$mean, x$se))
  cat(sprintf("  %d%% CI:           [%.2f, %.2f] ms\n",
              round(x$conf*100), x$ci[1], x$ci[2]))
  if (!is.null(x$recovery_bias))
    cat(sprintf("  Recovery bias:    %+.2f ms\n", x$recovery_bias))
  invisible(x)
}

#' @export
plot.ssrt_simulate <- function(x, ...) {
  hist(x$distribution, breaks=40, col="#59A14F", border="white",
       main="Parametric MC SSRT Distribution (ex-Gaussian)",
       xlab="SSRT (ms)", ylab="Frequency")
  abline(v=x$ssrt_obs, col="firebrick",  lwd=2, lty=1)
  abline(v=x$ci,       col="darkorange", lwd=2, lty=2)
  legend("topright", legend=c("Observed", paste0(round(x$conf*100),"% CI")),
         col=c("firebrick","darkorange"), lwd=2, lty=c(1,2), bty="n")
}


# =============================================================================
# 3.  ssrt_power()
# =============================================================================

#' Minimum trial count analysis via Monte Carlo
#'
#' Simulates datasets of increasing size and computes SSRT variance as a
#' function of trial count (power curve).
#'
#' @param data         data.frame for calibration.
#' @param trial_counts Stop-trial counts to evaluate.
#' @param n_iter       MC iterations per count. Default 500.
#' @param target_se    Target SE in ms. Default 10.
#' @param stop_col,rt_col,acc_col,ssd_col Column names.
#' @param seed         Random seed. Default 42.
#'
#' @return Object of class \code{ssrt_power}.
#'
#' @examples
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' p <- ssrt_power(d, trial_counts=c(10,30,50), n_iter=200)
#' print(p)
#'
#' @export
ssrt_power <- function(data,
                       trial_counts = c(10, 20, 30, 50, 75, 100, 150, 200),
                       n_iter       = 500,
                       target_se    = 10,
                       stop_col     = "vol",
                       rt_col       = "RT_exp",
                       acc_col      = "correct",
                       ssd_col      = "soa",
                       seed         = 42) {
  set.seed(seed)
  go_rts   <- data[[rt_col]][data[[stop_col]]==0]
  go_rts   <- go_rts[!is.na(go_rts) & go_rts > 0]
  exg      <- .fit_exgaussian(go_rts)
  ssrt_ref <- integration_adaptiveSSD(data, stop_col, rt_col, acc_col, ssd_col)
  mean_ssd <- mean(data[[ssd_col]][data[[stop_col]]==1], na.rm=TRUE)

  sim_n <- function(n_stop) {
    n_go   <- n_stop * 3L;  N <- n_go + n_stop
    go_all <- .rexgaussian(N, exg$mu, exg$sigma, exg$tau)
    stp_rt <- .rexgaussian(n_stop, ssrt_ref, exg$sigma*0.5, exg$tau*0.3)
    ssd_s  <- pmax(50, rnorm(n_stop, mean_ssd, 50))
    inh    <- go_all[seq_len(n_stop)] > (ssd_s + stp_rt)
    df <- data.frame(
      vol    = c(rep(1,n_stop), rep(0,n_go)),
      RT_exp = c(ifelse(inh, NA, go_all[seq_len(n_stop)]), go_all[(n_stop+1):N]),
      correct= c(as.integer(inh), rbinom(n_go, 1, 0.95)),
      soa    = c(ssd_s, rep(NA_real_, n_go)))
    tryCatch(
      integration_adaptiveSSD(df, "vol", "RT_exp", "correct", "soa"),
      error=function(e) NA_real_)
  }

  rows <- lapply(trial_counts, function(ns) {
    cat(sprintf("  Simulating n_stop = %d ...\n", ns))
    vals <- replicate(n_iter, sim_n(ns));  vals <- vals[!is.na(vals)]
    data.frame(n_stop_trials=ns, n_total=ns*4L,
               mean_ssrt=mean(vals), se=sd(vals),
               ci_width_95=diff(quantile(vals, c(0.025,0.975))),
               rmse=sqrt(mean((vals-ssrt_ref)^2)), n_valid=length(vals))
  })
  pt  <- do.call(rbind, rows)
  mn  <- pt$n_stop_trials[which(pt$se <= target_se)[1]]

  structure(list(power_table=pt, ssrt_true=ssrt_ref, exg_params=exg,
                 target_se=target_se, min_n=mn, n_iter=n_iter),
            class=c("ssrt_power","list"))
}

#' @export
print.ssrt_power <- function(x, ...) {
  cat("-- SSRTcalc Power Analysis -------------------------------------\n")
  cat(sprintf("  Reference SSRT:  %.2f ms\n", x$ssrt_true))
  cat(sprintf("  Target SE:       %.1f ms\n", x$target_se))
  if (!is.na(x$min_n))
    cat(sprintf("  Min stop trials: %d  (~%d total) to reach SE <= %.1f ms\n",
                x$min_n, x$min_n*4, x$target_se))
  else
    cat("  Target SE not reached within tested trial counts.\n")
  cat("\n");  print(x$power_table, digits=3, row.names=FALSE)
  invisible(x)
}

#' @export
plot.ssrt_power <- function(x, ...) {
  pt <- x$power_table
  op <- par(mfrow=c(1,2), mar=c(4,4,3,1));  on.exit(par(op))
  plot(pt$n_stop_trials, pt$se, type="b", pch=19, col="#4E79A7",
       xlab="Stop Trials", ylab="SE of SSRT (ms)", main="Precision")
  abline(h=x$target_se, col="firebrick", lty=2)
  if (!is.na(x$min_n)) abline(v=x$min_n, col="darkorange", lty=2)
  plot(pt$n_stop_trials, pt$rmse, type="b", pch=19, col="#E15759",
       xlab="Stop Trials", ylab="RMSE (ms)", main="Accuracy (RMSE)")
  abline(h=0, lty=3, col="grey50")
}


# =============================================================================
# 4.  ssrt_robustness()
# =============================================================================

#' Sensitivity of SSRT estimates to horse-race assumption violations
#'
#' Sweeps three violation types: go/stop process correlation, trigger failure
#' rate, and go-RT shift on stop trials.
#'
#' @param data          data.frame for calibration.
#' @param violation     "all", "correlation", "trigger_failure", or "go_shift".
#' @param n_iter        MC iterations per condition. Default 500.
#' @param corr_range    Correlations to test.
#' @param trigger_range Trigger failure rates to test.
#' @param shift_range   Go-RT shifts in ms to test.
#' @param stop_col,rt_col,acc_col,ssd_col Column names.
#' @param seed          Random seed. Default 42.
#'
#' @return Object of class \code{ssrt_robustness}.
#'
#' @examples
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' r <- ssrt_robustness(d, violation="trigger_failure", n_iter=20,
#'                      trigger_range=seq(0, 0.2, 0.1))
#' print(r)
#'
#' @export
ssrt_robustness <- function(data,
                            violation     = "all",
                            n_iter        = 500,
                            corr_range    = seq(0, 0.8, by=0.1),
                            trigger_range = seq(0, 0.4, by=0.05),
                            shift_range   = seq(0, 100, by=10),
                            stop_col      = "vol",
                            rt_col        = "RT_exp",
                            acc_col       = "correct",
                            ssd_col       = "soa",
                            seed          = 42) {
  set.seed(seed)
  violation <- match.arg(violation, c("all","correlation","trigger_failure","go_shift"))

  go_rts   <- data[[rt_col]][data[[stop_col]]==0]
  go_rts   <- go_rts[!is.na(go_rts) & go_rts > 0]
  exg      <- .fit_exgaussian(go_rts)
  ssrt_ref <- integration_adaptiveSSD(data, stop_col, rt_col, acc_col, ssd_col)
  mean_ssd <- mean(data[[ssd_col]][data[[stop_col]]==1], na.rm=TRUE)
  n_stop   <- sum(data[[stop_col]]==1)
  n_go     <- sum(data[[stop_col]]==0)

  run_sim <- function(corr=0, trigger_fail=0, go_shift=0) {
    bvn  <- .mvrnorm_simple(n_stop, c(0,0), matrix(c(1,corr,corr,1),2))
    u_go  <- pnorm(bvn[,1]);  u_stop <- pnorm(bvn[,2])

    qexg <- function(probs, mu, sigma, tau) {
      vapply(probs, function(p) {
        if (p<=0) return(mu - 3*sigma)
        if (p>=1) return(mu + tau + 3*sigma)
        tryCatch(
          uniroot(function(x) .pexgaussian(x,mu,sigma,tau) - p,
                  c(mu-4*sigma, mu+tau+4*sigma))$root,
          error=function(e) { s <- sort(.rexgaussian(1000,mu,sigma,tau)); s[max(1,round(p*1000))] })
      }, numeric(1))
    }

    go_rt_c  <- qexg(u_go,  exg$mu+go_shift, exg$sigma, exg$tau)
    stop_rt_c <- qexg(u_stop, ssrt_ref, exg$sigma*0.5, exg$tau*0.3)
    is_tf    <- runif(n_stop) < trigger_fail
    stp_eff  <- ifelse(is_tf, Inf, stop_rt_c)
    ssd_s    <- pmax(50, rnorm(n_stop, mean_ssd, 50))
    inh      <- go_rt_c > (ssd_s + stp_eff)
    go_rt_go <- .rexgaussian(n_go, exg$mu, exg$sigma, exg$tau)
    df <- data.frame(
      vol    = c(rep(1,n_stop), rep(0,n_go)),
      RT_exp = c(ifelse(inh, NA, go_rt_c), go_rt_go),
      correct= c(as.integer(inh), rbinom(n_go,1,0.95)),
      soa    = c(ssd_s, rep(NA_real_, n_go)))
    tryCatch(
      integration_adaptiveSSD(df, "vol", "RT_exp", "correct", "soa"),
      error=function(e) NA_real_)
  }

  cat("  Computing baseline...\n")
  base_mean <- mean(replicate(n_iter, run_sim()), na.rm=TRUE)

  sweep_one <- function(param_range, slot, label) {
    cat(sprintf("  Sweeping %s...\n", label))
    do.call(rbind, lapply(param_range, function(val) {
      args <- list(corr=0, trigger_fail=0, go_shift=0);  args[[slot]] <- val
      vals <- replicate(n_iter, do.call(run_sim, args));  vals <- vals[!is.na(vals)]
      data.frame(violation_value=val, mean_ssrt=mean(vals), se=sd(vals),
                 bias_from_baseline=mean(vals)-base_mean, n_valid=length(vals))
    }))
  }

  res <- list(baseline_ssrt=base_mean, ssrt_true=ssrt_ref, n_iter=n_iter)
  if (violation %in% c("all","correlation"))
    res$correlation    <- sweep_one(corr_range,    "corr",        "go/stop correlation")
  if (violation %in% c("all","trigger_failure"))
    res$trigger_failure <- sweep_one(trigger_range, "trigger_fail", "trigger failure rate")
  if (violation %in% c("all","go_shift"))
    res$go_shift       <- sweep_one(shift_range,   "go_shift",     "go-RT shift (ms)")

  structure(res, class=c("ssrt_robustness","list"))
}

#' @export
print.ssrt_robustness <- function(x, ...) {
  cat("-- SSRTcalc Robustness Analysis --------------------------------\n")
  cat(sprintf("  Reference SSRT:         %.2f ms\n", x$ssrt_true))
  cat(sprintf("  Baseline MC SSRT:       %.2f ms\n", x$baseline_ssrt))
  cat(sprintf("  Iterations/condition:   %d\n\n", x$n_iter))
  for (nm in c("correlation","trigger_failure","go_shift")) {
    if (!is.null(x[[nm]])) { cat(sprintf("  %s:\n", nm)); print(x[[nm]], digits=3, row.names=FALSE); cat("\n") }
  }
  invisible(x)
}

#' @export
plot.ssrt_robustness <- function(x, ...) {
  active <- intersect(c("correlation","trigger_failure","go_shift"), names(x))
  if (!length(active)) { message("No results to plot."); return(invisible(x)) }
  op <- par(mfrow=c(1,length(active)), mar=c(4,4,3,1));  on.exit(par(op))
  xlabs <- c(correlation="go/stop Correlation", trigger_failure="Trigger Failure Rate",
             go_shift="go-RT Shift (ms)")
  cols  <- c(correlation="#4E79A7", trigger_failure="#E15759", go_shift="#F28E2B")
  for (nm in active) {
    d <- x[[nm]]
    plot(d$violation_value, d$bias_from_baseline, type="b", pch=19, col=cols[nm],
         xlab=xlabs[nm], ylab="Bias in SSRT (ms)", main=paste("Robustness:", nm))
    abline(h=0, lty=2, col="grey50")
    polygon(c(d$violation_value, rev(d$violation_value)),
            c(d$bias_from_baseline+d$se, rev(d$bias_from_baseline-d$se)),
            col=adjustcolor(cols[nm],0.2), border=NA)
  }
}


# =============================================================================
# 5.  run_all_mc()
# =============================================================================

#' Run all four Monte Carlo analyses in one call
#'
#' @param data    data.frame in SSRTcalc long format.
#' @param n_iter  Iterations (shared). Default 1000.
#' @param seed    Random seed. Default 42.
#' @param stop_col,rt_col,acc_col,ssd_col Column names.
#'
#' @return Named list: bootstrap, simulation, power, robustness.
#'
#' @examples
#' \dontrun{
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' res <- run_all_mc(d, n_iter=500)
#' }
#'
#' @export
run_all_mc <- function(data, n_iter=1000, seed=42,
                       stop_col="vol", rt_col="RT_exp",
                       acc_col="correct", ssd_col="soa") {
  args <- list(data=data, n_iter=n_iter, seed=seed,
               stop_col=stop_col, rt_col=rt_col, acc_col=acc_col, ssd_col=ssd_col)
  # args without n_iter, for calls that need a different (smaller) value
  args_no_n <- args[setdiff(names(args), "n_iter")]

  cat("=== SSRTcalc MC Pipeline ===\n\n")
  cat("1/4  Bootstrap CIs...\n")
  boot <- do.call(ssrt_boot, args)
  cat("\n2/4  Parametric simulation...\n")
  sim  <- do.call(ssrt_simulate, args)
  cat("\n3/4  Power analysis...\n")
  pwr  <- do.call(ssrt_power,  c(args_no_n, list(n_iter=min(n_iter,500))))
  cat("\n4/4  Robustness analysis...\n")
  rob  <- do.call(ssrt_robustness, c(args_no_n, list(n_iter=min(n_iter,300))))
  cat("\n=== Done ===\n")
  list(bootstrap=boot, simulation=sim, power=pwr, robustness=rob)
}
