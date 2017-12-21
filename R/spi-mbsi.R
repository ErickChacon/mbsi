
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
#' @param time Time associated with the variable \code{y}.
#' @param tscale Time-scale to compute the MBSI in \code{time} units.
#' @param period Period (e.g. 53 weeks) defined to model the seasonal effect.
#'
#' @return A dataframe consisting of \code{y}, \code{time}, \code{season}, \code{mu},
#' \code{sigma}, \code{pzero}, \code{ecdf} and \code{spi}.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @importFrom utils getFromNamespace
#' @importFrom gamlss gamlss
#' @importFrom dplyr mutate
#' @importFrom gamlss.dist pZAGA ZAGA GA
#' @importFrom gamlss gamlss
#'
#' @export
mbsi <- function(y, time, tscale = 1, period = 365 / 7) {

  # # Compute moving average of y based on time-scale
  # y <- runmean(y, tscale)

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
  if (sum(y == 0) > 0) {
    lss <- gamlss::gamlss(form_mu, sigma.formula = form, nu.formula = form,
                          data = data0, family = gamlss.dist::ZAGA, trace = TRUE)
  } else {
    lss <- gamlss::gamlss(form_mu, sigma.formula = form,
                          data = data0, family = gamlss.dist::GA, trace = TRUE)
  }
  predict.gamlss <- getFromNamespace("predict.gamlss", "gamlss")
  mu <- predict.gamlss(lss, what = "mu", newdata = data, type = "response")
  sigma <- predict.gamlss(lss, what = "sigma", newdata = data, type = "response")
  if (sum(y == 0) > 0) {
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

  # Compute the spi
  data <- data %>% dplyr::mutate(
    ecdf0 = gamlss.dist::pZAGA(y, mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
    spi0 = qnorm(ecdf0),
    spi = runmean(spi0, width = tscale) * sqrt(tscale),
    ecdf = pnorm(spi)
    )

  class(data) <- c("mbsi", class(data))

  return(data)
}


#' @title Plot fitting when computing the SPI
#'
#' @description
#' \code{plot.mbsi} make a graph of the mean and coverage interval obtained when
#' computing the \code{spi}. This is useful to evaluate the seasonal behaviour and
#' the parameter estimation made by the classical SPI.
#'
#' @details
#' details.
#'
#' @param data An \code{mbsi} object returned by \code{spi}.
#'
#' @return A \code{ggplot} object. The graph is shown when this object is printed.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @examples
#' 
#'
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom gamlss.dist qZAGA
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line scale_colour_brewer theme
#' @importFrom ggplot2 element_blank
#'
#' @export
plot.mbsi <- function (x, which = c("fit", "ecdf"), binwidth = 0.05, ...) {

  data <- x
  which <- which[1]

  if (which == "fit") {
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
      ggplot2::theme(legend.position = "bottom",
                     legend.title = ggplot2::element_blank(),
                     axis.title.x = ggplot2::element_blank())
  } else if (which == "ecdf") {
    gg_mbsi <- data %>%
      ggplot2::ggplot(ggplot2::aes(ecdf)) +
      ggplot2::geom_histogram(binwidth = binwidth, center = binwidth/2)
  }

  return(gg_mbsi)
}



