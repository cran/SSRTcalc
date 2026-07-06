# =============================================================================
# R/core_ssrt.R
# Core SSRT estimation functions following Verbruggen et al. (2019)
# doi: 10.7554/eLife.46323
# =============================================================================

# -- Input validation helper --------------------------------------------------
.check_data <- function(data, stop_col, rt_col, ssd_col) {
  if (!is.data.frame(data))
    stop("`data` must be a data.frame.")
  for (col in c(stop_col, rt_col, ssd_col)) {
    if (!col %in% names(data))
      stop(sprintf("Column '%s' not found in data.", col))
  }
  if (!all(data[[stop_col]] %in% c(0, 1, NA)))
    stop(sprintf("`%s` must contain only 0 (go) and 1 (stop).", stop_col))
}

# -- Go-RT extractor ----------------------------------------------------------
.get_go_rts <- function(data, stop_col, rt_col, min_rt = 50) {
  go  <- data[data[[stop_col]] == 0, ]
  rts <- go[[rt_col]]
  rts <- rts[!is.na(rts) & rts >= min_rt]
  if (length(rts) < 5)
    stop("Fewer than 5 valid go RTs -- cannot estimate SSRT.")
  n_omit <- sum(is.na(go[[rt_col]]))
  pct    <- 100 * n_omit / nrow(go)
  if (pct > 10)
    warning(sprintf("%.1f%% go-trial omissions detected. SSRT may be underestimated.", pct))
  rts
}

# -- p(respond|signal) --------------------------------------------------------
.p_respond <- function(stp_trials, rt_col) {
  if (nrow(stp_trials) == 0) stop("No stop trials found.")
  n_sr <- sum(!is.na(stp_trials[[rt_col]]) & stp_trials[[rt_col]] > 0)
  n_sr / nrow(stp_trials)
}


#' Estimate SSRT via the integration method (adaptive / staircase SSD design)
#'
#' Implements the recommended integration method from Verbruggen et al. (2019).
#' For each dataset:
#' \enumerate{
#'   \item Compute p(respond|signal) from all stop trials.
#'   \item Find the nth percentile of the go-RT distribution (n = p_respond).
#'   \item Subtract the mean SSD: SSRT = nth_percentile_RT - mean(SSD).
#' }
#'
#' @param data     A data.frame with one row per trial.
#' @param stop_col Column name for the stop-trial indicator (1 = stop, 0 = go).
#'   Default \code{"vol"}.
#' @param rt_col   Column name for reaction time in ms. Default \code{"RT_exp"}.
#' @param acc_col  Column name for accuracy (1 = correct). Default \code{"correct"}.
#' @param ssd_col  Column name for stop-signal delay in ms. Default \code{"soa"}.
#' @param min_rt   Minimum valid RT in ms; shorter responses are excluded as
#'   anticipations. Default 50.
#'
#' @return A single numeric value: the estimated SSRT in milliseconds.
#'
#' @references
#' Verbruggen, F., et al. (2019). A consensus guide to capturing the ability
#' to inhibit actions and impulses: the stop-signal task. \emph{eLife}, 8,
#' e46323. \doi{10.7554/eLife.46323}
#'
#' @examples
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' integration_adaptiveSSD(d)
#'
#' @export
integration_adaptiveSSD <- function(data,
                                    stop_col = "vol",
                                    rt_col   = "RT_exp",
                                    acc_col  = "correct",
                                    ssd_col  = "soa",
                                    min_rt   = 50) {
  .check_data(data, stop_col, rt_col, ssd_col)

  go_rts     <- .get_go_rts(data, stop_col, rt_col, min_rt)
  stp        <- data[data[[stop_col]] == 1, ]
  p_resp     <- .p_respond(stp, rt_col)
  pctile     <- quantile(go_rts, probs = p_resp)
  mean_ssd   <- mean(stp[[ssd_col]], na.rm = TRUE)

  as.numeric(pctile - mean_ssd)
}


#' Estimate SSRT via the integration method (fixed SSD design)
#'
#' Identical to \code{\link{integration_adaptiveSSD}} in computation, but
#' intended for experiments using a fixed (constant) stop-signal delay. When
#' multiple fixed SSD values are used, SSRT is computed separately for each
#' SSD and the results are averaged (Verbruggen et al., 2019, Appendix).
#'
#' @inheritParams integration_adaptiveSSD
#'
#' @return A single numeric value: the estimated SSRT in milliseconds.
#'
#' @references
#' Verbruggen, F., et al. (2019). A consensus guide to capturing the ability
#' to inhibit actions and impulses: the stop-signal task. \emph{eLife}, 8,
#' e46323. \doi{10.7554/eLife.46323}
#'
#' @examples
#' data(fixed)
#' d <- fixed[fixed$SubjID == 1, ]
#' integration_fixedSSD(d)
#'
#' @export
integration_fixedSSD <- function(data,
                                 stop_col = "vol",
                                 rt_col   = "RT_exp",
                                 acc_col  = "correct",
                                 ssd_col  = "soa",
                                 min_rt   = 50) {
  .check_data(data, stop_col, rt_col, ssd_col)

  go_rts <- .get_go_rts(data, stop_col, rt_col, min_rt)
  stp    <- data[data[[stop_col]] == 1, ]

  unique_ssds <- sort(unique(stp[[ssd_col]]))

  if (length(unique_ssds) == 1) {
    # Single fixed SSD -- standard integration
    p_resp   <- .p_respond(stp, rt_col)
    pctile   <- quantile(go_rts, probs = p_resp)
    return(as.numeric(pctile - unique_ssds))
  }

  # Multiple fixed SSDs -- compute per-SSD and average
  ssrt_vec <- vapply(unique_ssds, function(ssd_val) {
    stp_ssd <- stp[stp[[ssd_col]] == ssd_val, ]
    if (nrow(stp_ssd) < 2) return(NA_real_)
    p_resp  <- .p_respond(stp_ssd, rt_col)
    pctile  <- quantile(go_rts, probs = p_resp)
    as.numeric(pctile - ssd_val)
  }, numeric(1))

  mean(ssrt_vec, na.rm = TRUE)
}


#' Estimate SSRT via the mean method (adaptive SSD design)
#'
#' Computes SSRT as the difference between mean go RT and mean SSD:
#' \deqn{SSRT = \bar{RT}_{go} - \bar{SSD}}
#' This method is less accurate than the integration method but is included for
#' comparison and historical compatibility.
#'
#' @inheritParams integration_adaptiveSSD
#'
#' @return A single numeric value: the estimated SSRT in milliseconds.
#'
#' @examples
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' mean_adaptiveSSD(d)
#'
#' @export
mean_adaptiveSSD <- function(data,
                             stop_col = "vol",
                             rt_col   = "RT_exp",
                             acc_col  = "correct",
                             ssd_col  = "soa",
                             min_rt   = 50) {
  .check_data(data, stop_col, rt_col, ssd_col)
  go_rts   <- .get_go_rts(data, stop_col, rt_col, min_rt)
  stp      <- data[data[[stop_col]] == 1, ]
  if (nrow(stp) == 0) stop("No stop trials found.")
  mean(go_rts) - mean(stp[[ssd_col]], na.rm = TRUE)
}


#' Estimate SSRT via the mean method (fixed SSD design)
#'
#' Computes SSRT as \eqn{\bar{RT}_{go} - \bar{SSD}}. For multiple fixed SSD
#' values the mean is taken across all stop trials.
#'
#' @inheritParams integration_adaptiveSSD
#'
#' @return A single numeric value: the estimated SSRT in milliseconds.
#'
#' @examples
#' data(fixed)
#' d <- fixed[fixed$SubjID == 1, ]
#' mean_fixedSSD(d)
#'
#' @export
mean_fixedSSD <- function(data,
                          stop_col = "vol",
                          rt_col   = "RT_exp",
                          acc_col  = "correct",
                          ssd_col  = "soa",
                          min_rt   = 50) {
  mean_adaptiveSSD(data, stop_col, rt_col, acc_col, ssd_col, min_rt)
}
