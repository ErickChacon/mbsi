#' @title mbsi: model-based standardized index
#'
#' @docType package
#' @aliases mbsi-package
#'
#' @description
#' The \code{mbsi} package provide tools to compute and visualize extreme
#' hydro-climatic events using the standardized precipitation index (SPI) and the
#' model-based standardized index (MBSI). The difference with between the MBSPI and
#' the classical SPI index is that it consider the association between continuous
#' times using cycle P-splines \code{pbc}. The package can also with work with
#' precipitation series containing missing values (NA), 0 or only non-zero values.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @useDynLib mbsi, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom magrittr %>%
"_PACKAGE"

if(getRversion() >= "2.15.1") utils::globalVariables(
  c(".",
    "y", "time",
    "mu", "sigma", "pzero",
    "ecdf", "spi",
    "ecdf0", "spi0",
    "q025", "q975",
    "data", "droughts", "floods", "event",
    "varname", "varvalue"
    )
  )

