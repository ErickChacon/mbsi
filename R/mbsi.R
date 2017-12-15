
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

  # Compute moving average of y based on time-scale
  y <- runmean(y, tscale)

  # Create data frame and add season
  data <- data.frame(y, time = 1:length(y))
  data <- transform(data, season = time %% period)
  data <- within(data, season[season == 0] <- period)

  # Fit gamlss
  # The optimization does not work when using gamlss::pbc, weird!!
  predict.gamlss <- getFromNamespace("pbc", "gamlss")
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
    ecdf = gamlss.dist::pZAGA(y, mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
    spi = qnorm(ecdf)
    )

  class(data) <- c("mbsi", class(data))

  return(data)
}

