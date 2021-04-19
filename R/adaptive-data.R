#' Adaptive-SSD dataset for stop-signal task
#'
#' Data from a mouse movement-based stop-signal experiment, using dynamically set stop-signal delays, with
#' random dot kinematogram as the "go" task, collected from 63 participants
#'
#' @docType data
#'
#' @usage data(adaptive)
#'
#' @format A data frame with 36288 rows and 7 variables
#' \describe{
#'  \item{new_id}{Subject identifier}
#'  \item{soa}{stop-signal delay (ignore for "go" trials)}
#'  \item{vol}{stop (1) or go (0) trial}
#'  \item{coh}{Percent coherent dots in the kinematogram go task}
#'  \item{RT_exp}{Response time in seconds, NA if no response was made}
#'  \item{correct}{Did the participant correctly respond in "go" trials/omit response in "stop" trials (1) or not(0)?}
#' }
#'
#' @keywords datasets
#'
#' @references Leontyev and Yamauchi (2019) PLoS One
#' (\doi{10.1371/journal.pone.0225437})
#'
#' @source \href{https://osf.io/g3fxs/}{OSF archive}
#'
#' @examples
#' data(adaptive)
#' head(adaptive)
"adaptive"
