#' Example adaptive stop-signal task data
#'
#' A long-format dataset from a stop-signal task using an adaptive
#' (staircase-tracking) stop-signal delay, suitable for demonstrating every
#' function in SSRTcalc.
#'
#' @format A data frame with 4000 rows (20 subjects x 200 trials) and 6
#'   columns:
#' \describe{
#'   \item{SubjID}{Integer subject identifier (1-20).}
#'   \item{trial}{Integer trial number within subject (1-200).}
#'   \item{vol}{Trial type: \code{0} = go trial, \code{1} = stop trial.}
#'   \item{RT_exp}{Observed reaction time in milliseconds. \code{NA} on
#'     successfully inhibited stop trials.}
#'   \item{correct}{Accuracy indicator (\code{1} = correct response on go
#'     trials, or successful inhibition on stop trials; \code{0} otherwise).}
#'   \item{soa}{Stop-signal delay (SOA) in milliseconds for stop trials;
#'     \code{NA} on go trials. Adjusted trial-by-trial by a staircase
#'     algorithm to track ~50\% inhibition.}
#' }
#'
#' @examples
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#' integration_adaptiveSSD(d)
#'
#' @source Simulated data distributed with the SSRTcalc package, generated to
#'   resemble a typical adaptive stop-signal task following
#'   Verbruggen et al. (2019) <doi:10.7554/eLife.46323>.
"adaptive"


#' Fixed-SSD stop-signal task data (motion discrimination)
#'
#' A long-format dataset from a stop-signal variant of a random-dot
#' motion-discrimination task in which the stop-signal delay (SOA) was set to
#' one of six \emph{fixed} values rather than being tracked by a staircase.
#' It complements the \code{\link{adaptive}} dataset and is intended for
#' demonstrating the fixed-SSD estimators \code{\link{integration_fixedSSD}}
#' and \code{\link{mean_fixedSSD}}.
#'
#' All timing variables are in \strong{milliseconds}, matching the rest of the
#' package and the \code{adaptive} dataset.
#'
#' @format A data frame with 28,799 rows (50 subjects) and 8 columns:
#' \describe{
#'   \item{SubjID}{Integer subject identifier. Subject 2 from the original
#'     recording is excluded (see Details), so the IDs run 1 and 3-51.}
#'   \item{trial}{Trial number within subject (1-576; a subject with one
#'     removed trial has 575).}
#'   \item{vol}{Trial type: \code{0} = go trial, \code{1} = stop trial.}
#'   \item{RT_exp}{Observed reaction time in milliseconds. \code{NA} when no
#'     response was made (successfully inhibited stop trials and go-trial
#'     omissions).}
#'   \item{correct}{Accuracy indicator (\code{1} = correct go response or
#'     successful inhibition; \code{0} otherwise).}
#'   \item{soa}{Stop-signal delay (SOA) in milliseconds for stop trials, one of
#'     six fixed values (100, 200, 300, 400, 500, 600 ms); \code{NA} on go
#'     trials, matching the \code{adaptive} dataset.}
#'   \item{coh}{Motion coherence for the go discrimination (0.1, 0.5, or 0.8),
#'     an experimental difficulty manipulation.}
#'   \item{response}{Response direction (\code{"left"}, \code{"right"}, or
#'     \code{NA} when no response was made).}
#' }
#'
#' @details
#' During preparation, subject 2 (who had no recorded responses on any trial)
#' and a single stop trial with a recorded reaction time of 0 ms were removed.
#' Subject identifiers otherwise follow the original recording, so identifier 2
#' is absent.
#'
#' @examples
#' data(fixed)
#' d <- fixed[fixed$SubjID == 1, ]
#' integration_fixedSSD(d)
#' mean_fixedSSD(d)
#'
#' @source Experimental data distributed with the SSRTcalc package, from a
#'   fixed-SSD motion-discrimination stop-signal task.
"fixed"
