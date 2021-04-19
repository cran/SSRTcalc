#' SSRT using integration method for studies with "fixed" method of setting SSD
#'
#' Estimating SSRT using integration method for studies that use fixed (randomly chosen on each trial from a pre-determined set) stop-signal delays.
#' @export
#' @param df  Dataframe with response time, accuracy, indication whether trial is stop or go, and delays for a given trial.
#' @param stop_col Name of the column in the dataframe \code{df} that indicates whether a given trial is a "stop" or a "go" trial ( 0 = go, 1 = stop)
#' @param rt_col Name of the column in the dataframe \code{df} that contains response time in seconds
#' @param acc_col Name of the column in the dataframe \code{df} that contains accuracy of inhibition ( 0 = incorrect, 1 = correct)
#' @param ssd_col Name of the column in the dataframe \code{df} that contains stop-signal delays
#' @param ssd_list List of stop-signal delays used in the experiment
#' @return SSRT corresponding to the \code{nth rt - ssd; n = p(respond|signal)*number of goRTs}
#' @examples
#' data(fixed)
#' sapply(split(fixed, fixed$new_id), integration_fixedSSD, stop_col = 'vol',acc_col ='acc',
#' rt_col = 'RT_exp', ssd_col = 'soa',ssd_list = c(0.1, 0.2,0.3, 0.4, 0.5, 0.6))





# Integration_fixed SSD
# The integration method assumes that the finishing time of the stop process
# corresponds to the nth RT,
# with n equal to the number of RTs in the RT distribution
# multiplied by the overall p(respond|signal)
# (Logan, 1981); SSRT can then be estimated by subtracting the mean
# SSD from the nth RT
# (taken from Verbruggen 2013)


integration_fixedSSD <- function(df, stop_col, rt_col, acc_col, ssd_col, ssd_list) {

  ssrt_list = c()
  for (ssd in ssd_list){
  stop_trials <- df[ which(df[,stop_col]==1), ]
  stop_trials = stop_trials[which (stop_trials[,ssd_col] == ssd),]

  stop_correct <- sum(stop_trials[,acc_col])
  overall_prob = 1 - stop_correct/nrow(stop_trials)
  go_trials = df[ which(df[,stop_col] == 0),]

  df1 <- go_trials[order(go_trials[,rt_col], na.last = NA) , ]
  nrt <- length(df1[,rt_col])
  nthindex = as.integer(round(nrt*overall_prob))
  nthrt <- df1[,rt_col][nthindex]
  ssrt_raw = nthrt - ssd
  if(isTRUE(ssrt_raw <= 0)){
    ssrt <- NA
  } else {
    ssrt <- ssrt_raw
  }

  ssrt_list <- append(ssrt_list, ssrt)
  }
  ssrt_final <- mean(ssrt_list, na.rm= TRUE)
  return (ssrt_final)
}


