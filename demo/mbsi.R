
data(simrain)

# Compute mbsi
spi_rain <- mbsi(simrain$rain, simrain$time)

# Visualize model fitting
plot(spi_rain)

# Visualize distribution of empirical cumulative density function
plot(spi_rain, which = "ecdf", binwidth = 0.05)

# Visualize extreme events
plot_extremes(spi_rain, threshold = 2)
