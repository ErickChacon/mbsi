
#' @title Classic standardized precipitation index (SPI)
#'
#' @description
#' \code{spi} computes the spi index using \code{gamlss} models.
#' It can also with work with precipitation series containing NA, 0 or
#' only non-zero values.
#'
#' @details
#' details.
#'
#' @param y Usually rainfall from which to compute the SPI, but other variables like
#' precipitation can also be used.
#' @param time Time indexed to the values of the argument \code{y}, it can be in weeks,
#' months, etc.
#' @param tscale Time-scale to compute the SPI in time units. This argument is used
#' to perform a moving average to the rainfall series before computing the SPI.
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
#' spi_rain <- spi(simrain$rain, simrain$time)
#'
#' # Visualize model fitting
#' plot(spi_rain)
#' # Visualize distribution of empirical cumulative density function
#' plot(spi_rain, which = "ecdf", binwidth = 0.05)
#' # Visualize extreme events
#' plot_extremes(spi_rain, threshold = 2)
#'
#' @importFrom dplyr group_by mutate select arrange
#' @importFrom tidyr nest unnest
#' @importFrom gamlss.dist pZAGA ZAGA
#' @importFrom gamlss gamlss lpred
#' @importFrom stats na.omit
#'
#' @export
spi <- function(y, time, tscale = 1, period = 52) {

  # Compute moving average of y based on time-scale
  season <- time %% period
  season[season == 0] <- period
  y <- runmean(y, tscale)

  # Create dataframe
  y_data <- data.frame(time, y, season)

  # Nesting data
  y_data <- y_data %>%
    dplyr::group_by(season) %>%
    tidyr::nest() %>%
    dplyr::mutate(mod = lapply(data, fit_iid_ZAGA)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::arrange(time)

  # Compute the spi and additional variables.
  y_data <- y_data %>% dplyr::mutate(
    q025 = gamlss.dist::qZAGA(rep(0.025, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
    q975 = gamlss.dist::qZAGA(rep(0.975, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
    mean = mu)

  y_data$ecdf <- NA
  ind <- !is.na(y)
  y_data <- y_data %>%
    within({
      ecdf[ind] <- gamlss.dist::pZAGA(y[ind], mu = mu[ind],
                                      sigma = 1 / sqrt(sigma[ind]), nu = pzero[ind])
      spi <- qnorm(ecdf)
    })

    y_data$qq_emp <- quantile(y_data$spi, y_data$ecdf, na.rm = TRUE)
    y_data$qq_the <- qnorm(y_data$ecdf, mean(y_data$spi, na.rm = TRUE), sd(y_data$spi, na.rm = TRUE))


  # y_data <- y_data %>% dplyr::mutate(
  #   ecdf = gamlss.dist::pZAGA(y, mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
  #   spi = qnorm(ecdf)
  #   )

  class(y_data) <- c("mbsi", class(y_data))

  return(y_data)
}

fit_iid_ZAGA <- function (data) {
  data0 <- na.omit(data)
  gam1 <- gamlss::gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                         data = data0, family = gamlss.dist::ZAGA, trace = FALSE)
  mu <- gamlss::lpred(gam1, what = "mu", type = "response")[1]
  sigma <- gamlss::lpred(gam1, what = "sigma", type = "response")[1]
  nu <- gamlss::lpred(gam1, what = "nu", type = "response")[1]
  data$mu <- mu
  data$sigma <- 1 / sigma ^ 2
  data$pzero <- nu
  return(data)
}
