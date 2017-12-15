
# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())
# library(day2day)
#
# CUSTOM FUNCTIONS -------------------------------------------------------------

gp1 <- function (t1, cov.model = NULL, cov.params = NULL) {
  coords <- cbind(t1)
  n <- nrow(coords)
  distance <- as.matrix(dist(coords))
  varcov <- do.call(cov.model, c(list(dis = distance), cov.params))
  right <- chol(varcov)
  output <- as.numeric(crossprod(right, rnorm(n)))
}

exp_cov <- function (dis, phi, sigma2) {
  sigma2 * exp(-dis/phi)
}

# SIMULATE RAINFALL ------------------------------------------------------------

weeks <- 1:52
years <- 1:10
data <- expand.grid(weeks = weeks, years = years)
data <- transform(data, time = (years - 1) * 52 + weeks)
data <- transform(data,
  season = sin(weeks / 52 * 2 * pi),
  gp = gp1(time, "exp_cov", list(phi = 1, sigma2 = 0.03)))
data <- transform(data, rain = day2day::rgamma_mu(520, exp(-2 + season + gp), exp(5)))
plot(data$time, data$rain, type = "b")
lines(data$time, exp(-2 + data$season), col = 2)
databla <- data

# COMPUTE CLASSIC SPI ----------------------------------------------------------

getwd()
setwd("..")

Rcpp::compileAttributes(getwd())
roxygen2::roxygenize(getwd(), roclets = c("collate", "namespace", "rd"))
devtools::install_local(getwd())

library(mbsi)
# library(tidyverse)
# library(dplyr)
bla <- mbsi::spi(data$rain, data$time)
ggbla <- plot(bla)
ggbla

library(mbsi)
bla <- mbsi::mbsi(data$rain, data$time, period = 52)
ggbla <- plot(bla)
ggbla

plot(bla$time, bla$y)
lines(bla$time, bla$mu)

plot(bla$season, bla$y)
