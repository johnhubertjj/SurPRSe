library(ggplot2)
library(grid)
library(dplyr)
library(lubridate)
library(reshape2)

#' Create some data to play with. Two time series with the same timestamp.
df <- data.frame(DateTime = ymd("2010-07-01") + c(0:8760) * hours(2), series1 = rnorm(8761), series2 = rnorm(8761, 100))

#' Create the two plots.
plot1 <- df %>%
  select(DateTime, series1) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series1), size = 0.5, alpha = 0.75) +
  ylab("Red dots / m") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

plot2 <- df %>%
  select(DateTime, series2) %>%
  na.omit() %>%
  ggplot() +
  geom_point(aes(x = DateTime, y = series2), size = 0.5, alpha = 0.75) +
  ylab("Blue drops / L") +
  theme_minimal() +
  theme(axis.title.x = element_blank())

grid.newpage()
grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "last"))

x <- seq(1, 5, length = 100)
y <- replicate(10, sin(2 * pi * x) + rnorm(100, 0, 0.3), "list")
z <- replicate(10, sin(2 * pi * x) + rnorm(100, 5, 0.3), "list")
y <- melt(y)
z <- melt(z)
df <- data.frame(x = y$Var1, rep = y$Var2, y = y$value, z = z$value)
dat <- melt(df, id = c("x", "rep"))