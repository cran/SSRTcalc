#' Fixed-SSD dataset for stop-signal task
#'
#' Data from a mouse movement-based stop-signal experiment, using preset stop-signal delays, with
#' random dot kinematogram as the "go" task, collected from 51 participants
#'
#' @docType data
#'
#' @usage data(fixed)
#'
#' @format A data frame with 29376 rows and 7 variables
#' \describe{
#'  \item{new_id}{Subject identifier}
#'  \item{soa}{stop-signal delay (ignore for "go" trials)}
#'  \item{vol}{stop (1) or go (0) trial}
#'  \item{coh}{Percent coherent dots in the kinematogram go task}
#'  \item{RT_exp}{Response time in seconds, NA if no response was made}
#'  \item{response}{Which button did the participant click?}
#'  \item{acc}{Did the participant respond in "go" trials/omit response in "stop" trials (1) or not(0)?}
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
#' data(fixed)
#' head(fixed)
"fixed"
