
#' @title Detect floods and droughts based on Mckee et al. (1993)
#'
#' @description
#' \code{find_flood_drough} detect events of floods and droughts based on the
#' definition proposed by Mckee et al. (1993) using a time series of standardized
#' precipitation values. A drought (flood) is defined  as a period of time in
#' which the SPI is continuously negative (positive) reaching at least
#' one value lower (higher) or equal to -1 (1).
#'
#' @details
#' Although the definition of Mckee et al. (1993) uses a threshold of 1 (-1). Other
#' values can be used for this threshold.
#'
#' @param spi A vector representing the time series of the standardized values.
#' @param threshold An numeric value used to detect the extreme events, 1 (-1) is
#' used by default as proposed by Mckee et al. (1993).
#' @param labels A character vector providing the labels used to return the three
#' types of events \code{normal}, \code{drought} and \code{flood}.
#'
#' @return A factor vector indicating the time of event the standardized
#' precipitation value corresponds.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' 
#'
#' @export
find_flood_drought <- function (spi, threshold = qnorm(1-0.05/2) ,
                                labels = c("normal", "drought", "flood")) {
  # Convert spi to string to detect extreme events.
  spi <- (spi > -threshold) + (spi > 0) + (spi > threshold)
  spi[is.na(spi)] <- 9
  spichar <- paste(spi, collapse = "")
  # Detect flood and drought.
  flood_exp <- gregexpr("[2-3]*3+[2-3]*", spichar)
  drought_exp <- gregexpr("[0-1]*0+[0-1]*", spichar)
  flood <- function(n) strrep("8", n)
  drought <- function(n) strrep("7", n)
  regmatches(spichar, flood_exp) <-
    Map(flood, lapply(regmatches(spichar, flood_exp), nchar))
  regmatches(spichar, drought_exp) <-
    Map(drought, lapply(regmatches(spichar, drought_exp), nchar))
  # Convert to vector.
  spichar <- substring(spichar, 1:nchar(spichar), 1:(nchar(spichar)))
  spi <- as.numeric(spichar)
  spi[spi == 9] <- NA
  spi[!(spi %in% c(7:8, NA))] <- 0
  spi <- factor(spi, c(0, 7, 8), labels)
  return(spi)
}

