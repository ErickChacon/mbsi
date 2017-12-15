
#' @title Compute the Standardized Precipitation Index.
#'
#' @description
#' \code{spi_week} computes the weekly spi index with gamlss models.
#' It can also with work with precipitation series containing NA, 0 or
#' only non-zero values.
#'
#' @details
#' details.
#'
#' @param rain Precipitation level.
#' @param tscale Time-scale to compute the SPI in week units.
#' @param period Period (e.g. 53 weeks) defined to model the seasonal effect.
#' @param package Package to use. \code{gamlss} for backfitting and \code{bamlss} for baysian gamlss models.
#' @param plot Logical value indicating if a plot should be draw or not.
#'
#' @return dataframe consisting of spi, mu, sigma and pzero.
#'
#' @author Erick A. Chacon-Montalvan
#'
#' @importFrom mgcv gam predict.gam
#' @importFrom gamlss pbc gamlss
#' @importFrom utils getFromNamespace
#' @importFrom graphics par points lines
#'
#' @export
spi_week <- function(rain, tscale = 1, period = 365 / 7, package = "gamlss", plot = FALSE) {

  # require(mgcv)

  # Compute moving average of rain based on time-scale.
  rain <- runmean(rain, tscale)

  # Create dataframe.
  data <- data.frame(rain, weeks = 1:length(rain))
  data <- transform(data, no_rain = (rain == 0) * 1, rain_level = rain)
  data <- within(data, rain_level[rain == 0] <- NA)
  data <- transform(data, weeks2 = weeks %% period) # for seasonal trend

  # Modelling probability of no rain.
  gam0 <- mgcv::gam(no_rain ~ s(weeks2, bs = "cc"), binomial("logit"), data)
  data$pzero <- mgcv::predict.gam(gam0, data, type = "response")

  if (package == "gamlss") {
    # require(gamlss)
    predict.gamlss <- getFromNamespace("predict.gamlss", "gamlss")
    form_mu <- rain_level ~ gamlss::pbc(weeks2)
    form_sg <- ~ gamlss::pbc(weeks2)
    data0 <<- na.omit(data)
    lss <- gamlss::gamlss(form_mu, sigma.formula = form_sg, data = data0, family = GA,
                  trace = FALSE)
    mu <- predict.gamlss(lss, what = "mu", newdata = data, type = "response")
    sigma <- predict.gamlss(lss, what = "sigma", newdata = data, type = "response")
    data$mu <- mu
    data$sigma <- 1 / sigma ^ 2
    # Notes:
    # 1) predict only works if data0 exists in the global env.
    # 2) stepGAIC only works if form_mu and form_sg exits in the global env.
    # Solution: use <<- instead of <- for global assignment.
  } else {
  }

  # Compute the spi and additional variables.
  data <- transform(data, cdf_gamma = pgamma(rain, sigma, sigma / mu))
  data <- transform(data, cdf_mix = pzero + (1-pzero) * cdf_gamma)
  data <- transform(data, spi = qnorm(data$cdf_mix),
                    q025 = qgamma((0.025 - pzero) / (1-pzero), sigma,
                                  sigma / mu),
                    q975 = qgamma((0.975 - pzero) / (1-pzero), sigma,
                                  sigma / mu))
  data <- within(data, {q025[is.na(q025)] <- 0;
                        q975[is.na(q975)] <- 0})



  if (plot) {
    opar <- par(mfrow = c(1, 2))
    plot(data$weeks, data$spi)
    points(data$weeks[abs(data$spi) > qnorm(0.975)],
           data$spi[abs(data$spi) > qnorm(0.975)], col = 2, lwd = 2)
    plot(data$weeks, data$rain)
    points(data$weeks[abs(data$spi) > qnorm(0.975)],
           data$rain[abs(data$spi) > qnorm(0.975)], col = 2, lwd = 2)
    lines(data$weeks, data$mu)
    lines(data$weeks, data$q025)
    lines(data$weeks, data$q975)
    par(opar)
  }

  return(subset(data, select = c(spi, mu, sigma, pzero)))
}

