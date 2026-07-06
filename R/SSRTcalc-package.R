#' SSRTcalc: Stop-Signal Reaction Time Calculator with Monte Carlo and
#' Bayesian Extensions
#'
#' Estimates stop-signal reaction time (SSRT) using the integration and mean
#' methods of Verbruggen et al. (2019), and extends these point estimates
#' with three families of tools:
#'
#' \itemize{
#'   \item \strong{Monte Carlo} (\code{\link{ssrt_boot}}, \code{\link{ssrt_simulate}},
#'     \code{\link{ssrt_power}}, \code{\link{ssrt_robustness}}): nonparametric
#'     bootstrap confidence intervals, parametric ex-Gaussian simulation,
#'     minimum-trial-count / power analysis, and sensitivity of SSRT to
#'     violations of the horse-race assumptions.
#'   \item \strong{Bayesian / Stan} (\code{\link{ssrt_stan}} and friends):
#'     single-subject and hierarchical ex-Gaussian horse-race models fit via
#'     Hamiltonian Monte Carlo, with an optional trigger-failure parameter
#'     following Matzke et al. (2013), posterior inhibition functions, and
#'     posterior predictive checks. Works with either the \pkg{cmdstanr} or
#'     \pkg{rstan} backend.
#'   \item \code{\link{run_all_mc}} and \code{\link{ssrt_stan_compare}}:
#'     convenience wrappers that run a full battery of analyses in one call.
#' }
#'
#' @section Typical workflow:
#' \preformatted{
#' data(adaptive)
#' d <- adaptive[adaptive$SubjID == 1, ]
#'
#' integration_adaptiveSSD(d)      # point estimate
#' ssrt_boot(d)                    # + bootstrap CI
#' ssrt_stan(d)                    # + full posterior (requires cmdstanr/rstan)
#' }
#'
#' @references
#' Verbruggen, F., Aron, A. R., Band, G. P. H., Beste, C., Bissett, P. G.,
#' Brockett, A. T., ... Boehler, C. N. (2019). A consensus guide to capturing
#' the ability to inhibit actions and impulses: the stop-signal task.
#' \emph{eLife}, 8, e46323. \doi{10.7554/eLife.46323}
#'
#' Matzke, D., Dolan, C. V., Logan, G. D., Brown, S. D., & Wagenmakers, E.-J.
#' (2013). Bayesian parametric estimation of stop-signal reaction time
#' distributions. \emph{Journal of Experimental Psychology: General}, 142(4),
#' 1047--1073. \doi{10.1037/a0030543}
#'
#' @keywords internal
#'
#' @importFrom stats density optim pnorm quantile rbinom rexp rnorm runif sd
#'   setNames uniroot var coef
#' @importFrom graphics abline axis hist legend lines par plot polygon segments
#' @importFrom grDevices adjustcolor
"_PACKAGE"


# =============================================================================
# `ranef` generic
#
# `ranef` is conventionally provided by nlme/lme4, not by base R. SSRTcalc
# defines its own minimal generic so that `ranef.ssrt_stan()` (see
# stan_extensions.R) dispatches correctly without depending on those packages.
# If nlme or lme4 is also loaded, R's S3 dispatch uses whichever generic is
# encountered first on the search path; both forward to `UseMethod("ranef")`
# with the same signature, so this is harmless either way.
# =============================================================================

#' Extract per-group random-effect summaries
#'
#' A generic function for extracting per-group (e.g. per-subject) parameter
#' summaries from a fitted model. SSRTcalc provides
#' \code{\link{ranef.ssrt_stan}} for hierarchical Stan fits.
#'
#' @param object A fitted model object.
#' @param ... Additional arguments passed to methods.
#'
#' @return A method-specific object (typically a data.frame of per-group
#'   summaries).
#'
#' @seealso \code{\link{ranef.ssrt_stan}}
#' @export
ranef <- function(object, ...) {
  UseMethod("ranef")
}
