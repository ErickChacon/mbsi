
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
#' @return A dataframe consisting of \code{y}, \code{time}, \code{spi}, \code{mu},
#' \code{sigma} and \code{pzero}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @importFrom dplyr group_by mutate select arrange
#' @importFrom tidyr nest unnest
#' @importFrom gamlss.dist pZAGA ZAGA
#' @importFrom gamlss gamlss lpred
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
    ecdf = gamlss.dist::pZAGA(y, mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
    spi = qnorm(ecdf)
    )

  class(y_data) <- c("mbsi", class(y_data))

  return(y_data)
}

fit_iid_ZAGA <- function (data) {
  data0 <- na.omit(data)
  gam1 <- gamlss::gamlss(y ~ 1, sigma.formula = ~ 1, nu.formula = ~ 1,
                         data = data0, family = gamlss.dist::ZAGA)
  mu <- gamlss::lpred(gam1, what = "mu", type = "response")[1]
  sigma <- gamlss::lpred(gam1, what = "sigma", type = "response")[1]
  nu <- gamlss::lpred(gam1, what = "nu", type = "response")[1]
  data$mu <- mu
  data$sigma <- 1 / sigma ^ 2
  data$pzero <- nu
  return(data)
}

#' @title Plot fitting 
#'
#' @description
#' \code{function} description.
#'
#' @details
#' details.
#'
#' @param par.
#'
#' @return return.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' 
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom gamlss.dist pZAGA ZAGA
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line scale_colour_brewer theme
#' @importFrom ggplot2 element_blank
#'
#' @export
ggplot.mbsi <- function (data) {

  data <- data %>%
    dplyr::mutate(
      q025 = gamlss.dist::qZAGA(0.025, mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
      q975 = gamlss.dist::qZAGA(0.975, mu = mu, sigma = 1 / sqrt(sigma), nu = pzero)
      )

  data_longer <- data %>%
    tidyr::gather(varname, varvalue, y, mu) %>%
    dplyr::mutate(
      varname = factor(
        varname,
        c("y", "mu"),
        c( "Moving average ", "Estimated mean"))
      )

  gg_mbsi <- data %>%
    ggplot2::ggplot(ggplot2::aes(time, y)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = q025, ymax = q975, fill = "95% Coverage interval"),
      col = rgb(1,0,0, 0.3), alpha = 0.1
      ) +
    ggplot2::geom_line(aes(y = varvalue, col = varname, linetype = varname),
      data_longer) +
    ggplot2::scale_colour_brewer(palette = "Set1", direction = -1) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank())

  return(gg_mbsi)
}



