
#' @title Compute the model-based standardized index (MBSI)
#'
#' @description
#' \code{mbsi} computes the model-based standardized index using \code{gamlss}
#' models. The difference with the classical \code{spi} index is that it consider the
#' association between continuous times using cycle P-splines \code{pbc}.
#' It can also with work with precipitation series containing NA, 0 or
#' only non-zero values.
#'
#' @details
#' details.
#'
#' @param y Precipitation level or other variable under study.
#' @param time Time associated with the variable \code{y}. This variable can be in
#' weeks, months, etc.
#' @param tscale Time-scale to compute the MBSI in \code{time} units. This argument
#' is used to return the standardized values at this scale. Usually defined to
#' identify different types of droughts of floods.
#' @param period Numeric value representing the period to define seasonality. For
#' example, 53 weeks, 12 months, 365 days. As it can be seen, it depends of the units
#' of the argument \code{time}.
#'
#' @return A dataframe consisting of \code{y}, \code{time}, \code{season}, \code{mu},
#' \code{sigma}, \code{pzero}, \code{ecdf} and \code{spi}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#'
#' data(simrain)
#' spi_rain <- mbsi(simrain$rain, simrain$time)
#'
#' # Visualize model fitting
#' plot(spi_rain)
#' # Visualize distribution of empirical cumulative density function
#' plot(spi_rain, which = "ecdf", binwidth = 0.05)
#' # Visualize extreme events
#' plot_extremes(spi_rain, threshold = 2)
#'
#' @importFrom utils getFromNamespace
#' @importFrom gamlss gamlss
#' @importFrom dplyr mutate
#' @importFrom gamlss.dist pZAGA ZAGA GA
#' @importFrom gamlss gamlss
#' @importFrom stats pnorm
#'
#' @export
mbsi <- function(y, time, tscale = 1, period = 365 / 7) {

  # Compute moving average of y based on time-scale
  y <- runmean(y, tscale)

  # Create data frame and add season
  data <- data.frame(y, time = 1:length(y))
  data <- transform(data, season = time %% period)
  data <- within(data, season[season == 0] <- period)

  # Fit gamlss
  # The optimization does not work when using gamlss::pbc, weird!!
  pbc <- getFromNamespace("pbc", "gamlss")
  form_mu <- y ~ pbc(season)
  form <- ~ pbc(season)
  data0 <<- na.omit(data)
  data0 <- na.omit(data)
  ZA <- sum(y == 0, na.rm = TRUE) > 0
  if (ZA) {
    lss <- gamlss::gamlss(form_mu, sigma.formula = form, nu.formula = form,
                          data = data0, family = gamlss.dist::ZAGA, trace = TRUE)
  } else {
    lss <- gamlss::gamlss(form_mu, sigma.formula = form,
                          data = data0, family = gamlss.dist::GA, trace = TRUE)
  }
  predict.gamlss <- getFromNamespace("predict.gamlss", "gamlss")
  mu <- predict.gamlss(lss, what = "mu", newdata = data, type = "response")
  sigma <- predict.gamlss(lss, what = "sigma", newdata = data, type = "response")
  if (ZA) {
    pzero <- predict.gamlss(lss, what = "nu", newdata = data, type = "response")
  } else {
    pzero <- 0
  }
  data$mu <- mu
  data$sigma <- 1 / sigma ^ 2
  data$pzero <- pzero
  # Notes:
  # 1) predict only works if data0 exists in the global env.
  # 2) stepGAIC only works if form_mu and form_sg exits in the global env.
  # Solution: use <<- instead of <- for global assignment.

  # Compute the spi and additional variables.
  data <- data %>% dplyr::mutate(
    q025 = gamlss.dist::qZAGA(rep(0.025, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
    q975 = gamlss.dist::qZAGA(rep(0.975, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
    mean = mu)

  data$ecdf <- NA
  ind <- !is.na(y)
  data <- data %>%
    within({
      ecdf[ind] <- gamlss.dist::pZAGA(y[ind], mu = mu[ind],
                                      sigma = 1 / sqrt(sigma[ind]), nu = pzero[ind])
      spi <- qnorm(ecdf)
    })

  # Compute the spi
  data <- data %>%
    dplyr::mutate(
      qq_emp = quantile(spi, ecdf, na.rm = TRUE),
      qq_the = qnorm(ecdf, mean(spi, na.rm = TRUE), sd(spi, na.rm = TRUE))
    )


  # # Compute the spi
  # data <- data %>% dplyr::mutate(
  #   ecdf = gamlss.dist::pZAGA(y, mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
  #   spi = qnorm(ecdf)
  #   # spi = runmean(spi0, width = tscale) * sqrt(tscale),
  #   # ecdf = pnorm(spi)
  #   )

  class(data) <- c("mbsi", class(data))

  return(data)
}


# #' @title Plot goodness of fit of the computation of the standardized precipitation
# #'
# #' @description
# #' \code{plot.mbsi} make graphs to evaluate the goodness of fit of the \code{gamlss}
# #' model. This is useful to compare the \code{spi} and the \code{mbsi}.
# #'
# #' @details
# #' Two options are provided. The default option is a graph of the mean and coverage
# #' interval obtained from the computation of \code{spi} or \code{mbsi}. This is
# #' useful to evaluate the seasonal behaviour and the parameter estimation. The second
# #' option is the histogram of the empirical cumulative density function to assess the
# #' probability integral transform.
# #'
# #' @param x A \code{mbsi} object returned by \code{spi} or \code{mbsi}.
# #' @param which A character value indicating which type of graph to return. The
# #' \strong{fit} option plots the mean and coverage interval of rainfall, while the
# #' \strong{ecdf} option plots an histogram of the empirical cumulative density
# #' function.
# #' @param binwidth The binwidth value for the histogram if \code{which == "ecdf"}.
# #' @param ... Additional arguments
# #'
# #' @return A \code{ggplot} object. The graph is shown when this object is printed.
# #'
# #' @author Erick A. Chacon-Montalvan
# #'
# #' @examples
# #'
# #' data(simrain)
# #' spi_rain <- mbsi(simrain$rain, simrain$time)
# #'
# #' # Visualize model fitting
# #' plot(spi_rain)
# #' # Visualize distribution of empirical cumulative density function
# #' plot(spi_rain, which = "ecdf", binwidth = 0.05)
# #'
# #' @importFrom dplyr mutate n
# #' @importFrom tidyr gather
# #' @importFrom gamlss.dist qZAGA
# #' @importFrom ggplot2 ggplot aes geom_ribbon geom_line scale_colour_brewer theme
# #' @importFrom ggplot2 element_blank
# #' @importFrom grDevices rgb
# #'
# #' @export
# plot.mbsi <- function (x, which = c("fit", "ecdf"), binwidth = 0.05, ...) {
#
#   data <- x
#   which <- which[1]
#
#   if (which == "fit") {
#     data <- data %>%
#       dplyr::mutate(
#         q025 = gamlss.dist::qZAGA(rep(0.025, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
#         q975 = gamlss.dist::qZAGA(rep(0.975, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero)
#         )
#
#     data_longer <- data %>%
#       tidyr::gather(varname, varvalue, y, mu) %>%
#       dplyr::mutate(
#         varname = factor(
#           varname,
#           c("y", "mu"),
#           c( "Moving average ", "Estimated mean"))
#         )
#
#     gg_mbsi <- data %>%
#       ggplot2::ggplot(ggplot2::aes(time, y)) +
#       ggplot2::geom_ribbon(
#         ggplot2::aes(ymin = q025, ymax = q975, fill = "95% Coverage interval"),
#         col = grDevices::rgb(1,0,0, 0.3), alpha = 0.1
#         ) +
#       ggplot2::geom_line(aes(y = varvalue, col = varname, linetype = varname),
#         data_longer) +
#       ggplot2::scale_colour_brewer(palette = "Set1", direction = -1) +
#       ggplot2::theme(legend.position = "bottom",
#                      legend.title = ggplot2::element_blank(),
#                      axis.title.x = ggplot2::element_blank())
#   } else if (which == "ecdf") {
#     gg_mbsi <- data %>%
#       ggplot2::ggplot(ggplot2::aes(ecdf)) +
#       ggplot2::geom_histogram(binwidth = binwidth, center = binwidth/2)
#   }
#
#   return(gg_mbsi)
# }
#
#
#
