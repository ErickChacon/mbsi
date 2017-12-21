#' @title Simulated rainfall data
#'
#' @description
#' Simulated rainfall data to be used in the examples of \code{mbsi} package. This
#' data has a seasonal component and a serially correlated component.
#'
#' @docType data
#'
#' @usage data(simrain)
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' data(simrain)
#' plot(simrain$time, simrain$rain, type = "b")
#' lines(simrain$time, exp(-2 + simrain$season), col = 2)
#'
"simrain"
