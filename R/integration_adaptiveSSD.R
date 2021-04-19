#' SSRT using integration method for studies with "adaptive" method of setting SSD
#'
#' Estimating SSRT using integration method for studies that use adaptive (increasing/decreasing by a given increment) stop-signal delays.
#' @export
#' @param df  Dataframe with response time, accuracy, indication whether trial is stop or go, and delays for a given trial.
#' @param stop_col Name of the column in the dataframe \code{df} that indicates whether a given trial is a "stop" or a "go" trial ( 0 = go, 1 = stop)
#' @param rt_col Name of the column in the dataframe \code{df} that contains response time in seconds
#' @param acc_col Name of the column in the dataframe \code{df} that contains accuracy of inhibition ( 0 = incorrect, 1 = correct)
#' @param ssd_col Name of the column in the dataframe \code{df} that contains stop-signal delays
#' @return SSRT corresponding to the \code{ nth rt - ssd; n = p(respond|signal)*number of goRTs}
#' @examples
#' data(adaptive)
#' sapply(split(adaptive, adaptive$new_id), integration_adaptiveSSD, stop_col = 'vol',
#' ssd_col = 'soa', rt_col = 'RT_exp', acc_col = 'correct')





# Integration_adaptiveSSD
# The integration method assumes that the finishing time of the stop process
# corresponds to the nth RT,
# with n equal to the number of RTs in the RT distribution
# multiplied by the overall p(respond|signal)
# (Logan, 1981); SSRT can then be estimated by subtracting the mean
# SSD from the nth RT
# (taken from Verbruggen 2013)



integration_adaptiveSSD <- function(df, stop_col, rt_col, acc_col, ssd_col) {
  go_trials = df[ which(df[,stop_col] == 0),]
  stop_trials <- df[ which(df[,stop_col]==1), ]
  stop_count <- sum(stop_trials[,acc_col])
  overall_prob = 1 - stop_count/nrow(stop_trials)
  df1 <- go_trials[order(go_trials[,rt_col], na.last = NA) , ]
  nrt <- length(df1[,rt_col])
  nthindex = as.integer(round(nrt*overall_prob))
  meanssd = mean(stop_trials[, ssd_col], na.rm =TRUE)
  nthrt <- df1[,rt_col][nthindex]
  ssrt_raw <- nthrt - meanssd

  if(isTRUE(ssrt_raw <= 0)){
    ssrt = NA
  } else {
    ssrt = ssrt_raw
  }
  return(ssrt)
}





