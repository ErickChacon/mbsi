mbsi_fit <- function(y, time, period = 365 / 7, size = NULL) {

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
  ZA <- sum(y == 0) > 0 # check if zero augmented
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

  # Compute the spi
  data <- data %>%
    dplyr::mutate(
      q025 = gamlss.dist::qZAGA(rep(0.025, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
      q975 = gamlss.dist::qZAGA(rep(0.975, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
      mean = mu,
      ecdf = gamlss.dist::pZAGA(y, mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
      spi = qnorm(ecdf),
      qq_emp = quantile(spi, ecdf, na.rm = TRUE),
      qq_the = qnorm(ecdf, mean(spi, na.rm = TRUE), sd(spi, na.rm = TRUE))
      # spi = runmean(spi0, width = tscale) * sqrt(tscale),
      # ecdf = pnorm(spi)
    )

  # samples
  if (!is.null(size)) {
    samples <- with(data,
                    replicate(size,
                              rZAGA(length(mu), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero))
                    )
  } else {
    samples <- NULL
  }

  # add atributes
  attr(data, "model") <- lss
  attr(data, "samples") <- samples

  class(data) <- c("mbsi", class(data))

  return(data)
}

mbsi_sample <- function (x, size = 1000) {
  samples <- with(x,
              replicate(size,
                        rZAGA(length(mu), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero))
            )
  attr(x, "samples") <- samples
  return(x)
}

mbsi_ma <- function (x, tscale = 1) {

  samples_ma <- apply(attr(x, "samples"), 2, function (x) runmean(x, tscale))
  x$y <- runmean(x$y, tscale)

  lprob <- function (i) {
    sum(samples_ma[i, ] < x$y[i]) / length(samples_ma[i, ])
  }

  x$q025 <- apply(samples_ma, 1, function (x) quantile(x, 0.025, na.rm = TRUE))
  x$q975 <- apply(samples_ma, 1, function (x) quantile(x, 0.975, na.rm = TRUE))
  x$mean <- apply(samples_ma, 1, function (x) mean(x))
  x$ecdf <- sapply(1:nrow(x), lprob)
  if (min(x$ecdf, na.rm = TRUE) == 0) {
    x$ecdf[x$ecdf == 0] <- min(x$ecdf[x$ecdf != 0], na.rm = TRUE)
  }
  if (min(x$ecdf, na.rm = TRUE) == 1) {
    x$ecdf[x$ecdf == 1] <- max(x$ecdf[x$ecdf != 1], na.rm = TRUE)
  }
  x$spi <- qnorm(x$ecdf)

  # Compute the spi
  x <- x %>%
    dplyr::mutate(
      qq_emp = quantile(spi, ecdf, na.rm = TRUE),
      qq_the = qnorm(ecdf, mean(spi, na.rm = TRUE), sd(spi, na.rm = TRUE))
    )

  attr(x, "samples") <- samples_ma

  class(x) <- c("mbsi", class(x))

  return(x)
}


#' @title Plot goodness of fit of the computation of the standardized precipitation
#'
#' @description
#' \code{plot.mbsi} make graphs to evaluate the goodness of fit of the \code{gamlss}
#' model. This is useful to compare the \code{spi} and the \code{mbsi}.
#'
#' @details
#' Two options are provided. The default option is a graph of the mean and coverage
#' interval obtained from the computation of \code{spi} or \code{mbsi}. This is
#' useful to evaluate the seasonal behaviour and the parameter estimation. The second
#' option is the histogram of the empirical cumulative density function to assess the
#' probability integral transform.
#'
#' @param x A \code{mbsi} object returned by \code{spi} or \code{mbsi}.
#' @param which A character value indicating which type of graph to return. The
#' \strong{fit} option plots the mean and coverage interval of rainfall, while the
#' \strong{ecdf} option plots an histogram of the empirical cumulative density
#' function.
#' @param binwidth The binwidth value for the histogram if \code{which == "ecdf"}.
#' @param ... Additional arguments
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
#' # Visualize model fitting
#' plot(spi_rain)
#' # Visualize distribution of empirical cumulative density function
#' plot(spi_rain, which = "ecdf", binwidth = 0.05)
#'
#' @importFrom dplyr mutate n
#' @importFrom tidyr gather
#' @importFrom gamlss.dist qZAGA
#' @importFrom ggplot2 ggplot aes geom_ribbon geom_line scale_colour_brewer theme
#' @importFrom ggplot2 element_blank
#' @importFrom grDevices rgb
#'
#' @export
plot.mbsi <- function (x, which = c("fit", "ecdf", "qq"), binwidth = 0.05, timevar = time, ...) {

  data <- x
  which <- which[1]

  if (which == "fit") {
    # data <- data %>%
    #   dplyr::mutate(
    #     q025 = gamlss.dist::qZAGA(rep(0.025, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero),
    #     q975 = gamlss.dist::qZAGA(rep(0.975, n()), mu = mu, sigma = 1 / sqrt(sigma), nu = pzero)
    #     )

    data_longer <- data %>%
      tidyr::gather(varname, varvalue, y, mean) %>%
      dplyr::mutate(
        varname = factor(
          varname,
          c("y", "mean"),
          c( "Moving average ", "Estimated mean"))
        )

    gg_mbsi <- data %>%
      ggplot2::ggplot(ggplot2::aes(!!ensym(timevar), y)) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = q025, ymax = q975, fill = "95% Coverage interval"),
        col = grDevices::rgb(1,0,0, 0.3), alpha = 0.1
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
  } else if (which == "qq") {
    # data$qq_emp <- quantile(data$spi, data$ecdf, na.rm = TRUE)
    # data$qq_the <- qnorm(data$ecdf, mean(data$spi, na.rm = TRUE), sd(data$spi, na.rm = TRUE))
    gg_mbsi <- data %>%
      ggplot2::ggplot(ggplot2::aes(qq_the, qq_emp)) +
      ggplot2::geom_abline(intercept = 0, slope = 1, col = 2, size = 1) +
      ggplot2::geom_point(size = 0.5)
  }

  return(gg_mbsi)
}



