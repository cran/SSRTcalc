
#'
#' Estimating SSRT using mean method for studies that use fixed (randomly chosen on each trial from a pre-determined set) stop-signal delays
#'
#' @export
#' @param df  Dataframe with response time, accuracy, indication whether trial is stop or go, and delays for a given trial.
#' @param stop_col Name of the column in the dataframe \code{df} that indicates whether a given trial is a "stop" or a "go" trial ( 0 = go, 1 = stop)
#' @param rt_col Name of the column in the dataframe \code{df} that contains response time in seconds
#' @param acc_col Name of the column in the dataframe \code{df} that contains accuracy of inhibition ( 0 = incorrect, 1 = correct)
#' @param ssd_col Name of the column in the dataframe \code{df} that contains stop-signal delays
#' @param ssd_list List of stop-signal delays used in the experiment
#' @return  Stop-signal reaction time corresponding roughly to 50 percent inhibition accuracy.
#' @examples
#' data(fixed)
#' sapply(split(fixed, fixed$new_id), mean_fixedSSD, stop_col = 'vol',acc_col ='acc',
#' rt_col = 'RT_exp', ssd_col = 'soa',ssd_list = c(0.1, 0.2,0.3, 0.4, 0.5, 0.6))
#'
#'
#'
#'

mean_fixedSSD <- function(df, stop_col, rt_col, acc_col, ssd_col, ssd_list) {
try({
  stop_trials <- df[ which(df[,stop_col]==1),]
  go_trials <- df[ which(df[,stop_col] == 0),]
  acc_soas <- data.frame()
  soas = ssd_list
  for (soa in soas){
    stp_tr <- stop_trials[ which(stop_trials[,ssd_col]==soa), ]
    total_stops <- sum(stp_tr[,acc_col])
    total_trials <- nrow(stp_tr)
    ssd_acc <- total_stops/total_trials
    acc_soas = rbind(acc_soas, data.frame(soa, ssd_acc))
  }

  f_of_ssd <- stats::splinefun(soa,ssd_acc)

  z <- f_of_ssd(0.5, deriv = 0) #evaluate the interpolating cubic spline (deriv = 0)
  meanRTGO <- mean(go_trials[, rt_col], na.rm = TRUE)
  ssrt_raw = meanRTGO - z
  if(isTRUE(ssrt_raw <= 0)){
    ssrt <- NA
  } else {
    ssrt <- ssrt_raw
  }

  return(ssrt)
})
}
