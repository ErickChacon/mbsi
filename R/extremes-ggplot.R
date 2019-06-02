
# #' @param threshold Threshold to plot the extreme events.
StatEvents <- ggplot2::ggproto("StatEvents", ggplot2::Stat,
  required_aes = c("x", "y"),
  compute_group = function(data, scales, threshold = 0) {
    transition <- data %>%
      # dplyr::select(-event) %>%
      dplyr::arrange(x) %>%
      dplyr::mutate(
        ini = c(abs(diff(y > threshold)), NA), ini = cumsum_na(ini) * ini,
        end = c(NA, abs(diff(y > threshold))), end = cumsum_na(end) * end) %>%
      tidyr::gather(ini_end, change_id,  ini:end) %>%
      dplyr::filter(change_id > 0) %>%
      dplyr::group_by(change_id) %>%
      dplyr::do(data.frame(x = predict(lm(x ~ y, .), data.frame(y = threshold)))) %>%
      dplyr::mutate(y = threshold) %>% dplyr::ungroup() %>% dplyr::select(- change_id)

    data$y[data$event == 0] <- NA
    # data$y[data$y < threshold] <- NA
    data <- data.frame(x = c(data$x, transition$x),
                       y = c(data$y, transition$y),
                       ymin = threshold,
                       ymax = c(data$y, transition$y))
  }
)

#' @title Extreme events.
#'
#' @description
#' \code{stat_events} computes a data frame to adequatly draw areas
#' corresponding to extreme events in a time series.
#'
#' @param mapping Aesthetic mapping created by \code{aes} or \code{aes}.
#' @param data Dataset to use.
#' @param geom Geometry.
#' @param position Position.
#' @param na.rm Logical value to remove missing values.
#' @param show.legend Show legend.
#' @param inherit.aes Inherit aesthetics.
#' @param ... Additional arguments passed to \code{layer}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' # Obtain a zero-mean time series
#' library(ggplot2)
#' economics <- transform(economics,
#'   unemploy_zero = 2 * (unemploy - mean(unemploy)) / sd(unemploy))
#'
#' # Plot time series with appropiate shade for negative and positive values
#' ggplot(economics, aes(date, unemploy_zero)) +
#'   geom_line() +
#'   geom_point(color = 2, alpha = 0.3) +
#'   stat_events(alpha = 0.3)
#'
#' # Plot time series with appropiate shade for a given threshold
#' ggplot(economics, aes(date, unemploy_zero)) +
#'   geom_line() +
#'   geom_point(color = 2, alpha = 0.3) +
#'   stat_events(threshold = 0.5, alpha = 0.3)
#'
#' # Plot time series with appropiate shade for a given threshold
#' ggplot(economics, aes(date, unemploy_zero)) +
#'   geom_line() +
#'   stat_events(aes(event = I(1 * (unemploy_zero > 1)), fill = "positive peak"),
#'               threshold = 1, alpha = 0.3) +
#'   stat_events(aes(event = I(1 * (unemploy_zero < -1)), fill = "negative peak"),
#'               threshold = -1, alpha = 0.3)
#'
#' @importFrom ggplot2 layer
#'
#' @export
stat_events <- function(mapping = NULL, data = NULL, geom = "ribbon",
                        position = "identity", na.rm = FALSE, show.legend = NA,
                        inherit.aes = TRUE, ...) {
  ggplot2::layer(
    stat = StatEvents, data = data, mapping = mapping, geom = geom,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}


#' @title Cumsum vector with NA values.
#'
#' @description
#' \code{cumsum_na} Cumulative sum for vectors with NA values.
#'
#' @details
#' When missing values are present, \code{cumsum_na} replace missing values for 0 and
#' compute the usual \code{cumsum}.
#'
#' @param x The vector for which the cumulative sum is desired.
#' @param ... Additional arguments for the \code{cumsum} function.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' x <- 1:10
#' x[3] <- NA
#' cumsum_na(x)
#' @export
cumsum_na <- function (x, ...) {
  x[is.na(x)] <- 0
  cumsum(x, ...)
}


#' @title Plot floods and droughts episodes using standardized precipitation values
#'
#' @description
#' \code{plot_extremes} Plot floods and droughts episodes using the function
#' \code{find_flood_drought}. This function employs the
#' definition proposed by Mckee et al. (1993) using a time series of standardized
#' precipitation values. A drought (flood) is defined  as a period of time in
#' which the SPI is continuously negative (positive) reaching at least
#' one value lower (higher) or equal to -1 (1).
#'
#' @details
#' Although the definition of Mckee et al. (1993) uses a threshold of 1 (-1). Other
#' values can be used for this threshold.
#'
#' @param data A \code{mbsi} object returned by \code{spi} or \code{mbsi}, which
#' contains the standardized precipitation values.
#' @param threshold An numeric value used to detect the extreme events, 1 (-1) is
#' used by default as proposed by Mckee et al. (1993).
#'
#' @return A \code{ggplot} object. The graph is shown when this object is printed.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' data(simrain)
#' spi_rain <- mbsi(simrain$rain, simrain$time)
#'
#' # Visualize extreme events
#' plot_extremes(spi_rain, threshold = 2)
#' plot_extremes(spi_rain, threshold = 1.5)
#' plot_extremes(spi_rain, threshold = 1)
#'
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_line geom_hline scale_fill_brewer theme
#' @importFrom ggplot2 element_blank
#'
#' @export
plot_extremes <- function (data, threshold = 1, timevar = time) {

  # Detect extreme events
  data <- data %>% dplyr::mutate(
    event = find_flood_drought(spi, threshold),
    floods = 1 * (event == "flood"),
    droughts = 1 * (event == "drought")
    )

  gg_extreme <- ggplot2::ggplot(data, ggplot2::aes(!!ensym(timevar), spi)) +
    ggplot2::geom_line(size = 0.3) +
    ggplot2::geom_hline(yintercept = c(-1, 1) * threshold, linetype = 2, alpha = 0.3) +
    stat_events(aes(event = floods, fill = "Flood"), alpha = 0.6) +
    stat_events(aes(event = droughts, fill = "Drought     "), alpha = 0.7) +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::labs(y = "standardized value") +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank())

  return(gg_extreme)
}
