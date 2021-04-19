#' SSRT using mean method for studies with "adaptive" method of setting SSD
#'
#' Estimating SSRT using mean method for studies that use adaptive (increasing/decreasing by a given increment) stop-signal delays
#' @export
#' @param df  Dataframe with response time, accuracy, indication whether trial is stop or go, and delays for a given trial.
#' @param stop_col Name of the column in the dataframe \code{df} that indicates whether a given trial is a "stop" or a "go" trial ( 0 = go, 1 = stop)
#' @param rt_col Name of the column in the dataframe \code{df} that contains response time in seconds
#' @param ssd_col Name of the column in the dataframe \code{df} that contains stop-signal delays
#' @return Spline-interpolated stop-signal reaction time corresponding roughly to 50% inhibition accuracy.
#' @examples
#' data(adaptive)
#' sapply(split(adaptive, adaptive$new_id), mean_adaptiveSSD, stop_col = 'vol',
#' ssd_col = 'soa', rt_col = 'RT_exp')



# "The mean method uses the mean of the inhibition function,
# which corresponds to the average SSD obtained with
# the tracking procedure when p(respond|signal) = .50.
# In other words, the mean method assumes that the
# mean RT equals SSRT plus the mean SSD, so
# SSRT can be estimated easily by subtracting the mean SSD
#  from the mean RT" (Verbruggen 2013)



mean_adaptiveSSD <- function (df, rt_col, ssd_col, stop_col)
{
  stop_trials <- df[ which(df[,stop_col]==1), ]
  go_trials <- df[ which(df[,stop_col] == 0),]
  meanRT <-mean(go_trials[,rt_col], na.rm = TRUE)
  meanSSD <- mean(stop_trials[,ssd_col])
  ssrt_raw <- meanRT - meanSSD
  if(isTRUE(ssrt_raw <= 0)){
    ssrt <- NA
  } else {
    ssrt <- ssrt_raw
  }
  return(ssrt)

}
