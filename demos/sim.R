
# CLEAN WORKSPACE AND LOAD PACKAGES --------------------------------------------

rm(list = ls())
setwd("..")
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


set.seed(1)
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
simrain <- data
save(simrain, file = file.path("data", "simrain.RData"))


# COMPUTE CLASSIC SPI ----------------------------------------------------------

getwd()
setwd("..")

Rcpp::compileAttributes(getwd())
roxygen2::roxygenize(getwd(), roclets = c("collate", "namespace", "rd"))
devtools::install_local(getwd())

devtools::check()

library(mbsi)
bla <- mbsi::spi(data$rain, data$time)
plot(bla)
plot(bla, which = "ecdf", binwidth = 0.05)
plot_extremes(bla, 2)

bla <- mbsi::mbsi(data$rain, data$time, period = 52, tscale = 8)
plot(bla)
plot(bla, which = "ecdf")
plot_extremes(bla, 2)

tools::showNonASCII(readLines("man/find_flood_drought.Rd"))


t.data.frame <- function(x) {
    x <- as.matrix(x)
    NextMethod("t")
}

ok <- data.frame(x = 1:5)

ok <- as.character(1:5)
sum.character <- function(x) {
  x <- as.numeric(x)
  NextMethod("sum")
}
sum.character(ok)

sum(ok)

x <- structure(1, class = letters)
bar <- function(x) UseMethod("bar", x)
bar.z <- function(x) "z"
bar(x)


baz <- function(x) UseMethod("baz", x)
baz.A <- function(x) "A"
baz.B <- function(x) "B"

ab <- structure(1, class = c("A", "B"))
ba <- structure(1, class = c("B", "A"))
baz(ab)
baz(ba)

baz.C <- function(x) c("C", NextMethod())
ca <- structure(1, class = c("C", "A"))
cb <- structure(1, class = c("C", "B"))
baz(ca)
baz(cb)


