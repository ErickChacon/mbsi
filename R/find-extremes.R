
#' @title Identify floods and droughts based on SPI.
#'
#' @description
#' \code{find_flood_drought} identifies floods and droughts for a vector series of
#' standardized precipitation values.
#'
#' @details
#' Floods (droughts) are identified by a subset of negative values where at least one value
#' exceeds -1 (1).
#'
#' @param spi A vector series of standardized precipitation values.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @export
find_flood_drought <- function (spi) {
  # Convert spi to string to detect extreme events.
  spi <- (spi > -1) + (spi > 0) + (spi > 1)
  spi[is.na(spi)] <- 9
  spichar <- spi %>% paste(collapse = "")
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
  spi <- factor(spi, c(0, 7, 8), c("normal", "drought", "flood"))
  return(spi)
}
