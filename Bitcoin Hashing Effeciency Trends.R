################################################################################
# Libraries
################################################################################
library(lubridate)
library(latex2exp)
library(ggplot2)


################################################################################
# Import Data
################################################################################

# Microprocessor trend data is available from:
# https://github.com/karlrupp/microprocessor-trend-data.git
# Specify the cloned directory in config.yml under default:

processor_trends <- config::get("processor_trends")
cores <- read.table(paste(processor_trends,"48yrs/cores.dat", sep = "/"))
frequency <- read.table(
  paste(processor_trends,"48yrs/frequency.dat", sep = "/"))
specint <- read.table(paste(processor_trends,"48yrs/specint.dat", sep = "/"))
transitors <- read.table(
  paste(processor_trends,"48yrs/transistors.dat", sep = "/"))
watts <- read.table(paste(processor_trends,"48yrs/watts.dat", sep = "/"))
mining_efficiency <- read.csv("data/mining hash rate.csv", comment.char = "#")


# Properly format the imported data
names(cores)[names(cores) == "V1"] <- "t"
names(frequency)[names(frequency) == "V1"] <- "t"
names(specint)[names(specint) == "V1"] <- "t"
names(transitors)[names(transitors) == "V1"] <- "t"
names(watts)[names(watts) == "V1"] <- "t"

names(cores)[names(cores) == "V2"] <- "y"
names(frequency)[names(frequency) == "V2"] <- "y"
names(specint)[names(specint) == "V2"] <- "y"
names(transitors)[names(transitors) == "V2"] <- "y"
names(watts)[names(watts) == "V2"] <- "y"

mining_efficiency[,4] <- as.Date(mining_efficiency[, 4], format = "%d/%m/%Y")


################################################################################
# Fit Models
################################################################################

# Models
sigmoid = function(params, x) {
  params[1] / (1 + exp(-params[2] * (x - params[3])))
}

gen_logistic = function(params, x) {
  params[1] - (params[1] - params[2]) /
    (1 + exp(-params[3] * (x - params[4])))
}

# General Processor Trends
params_spec <- nls(log(y) ~ log(a/(1 + exp(-b * (t - c)))),
                   data = specint,
                   start = list(a = 5.93e+04, b = 0.3695, c = 2008))
params_freq <- nls(log(y) ~ log(a/(1 + exp(-b * (t - c)))),
                   data = frequency,
                   start = list(a = 3251, b = 0.27, c = 2005))
params_watt <- nls(log(y) ~ log(a/(1 + exp(-b * (t - c)))),
                   data = watts,
                   start = list(a = 255, b = 0.17, c = 2010))
params_core <- nls(log(y) ~ log(1 + a / (1 + exp(-b * (t - c)))), data = cores,
                   start = list(a = 95, b = 0.4, c = 2023),
                   algorithm = "port")

# Bitcoin Mining Efficiency
hash = data.frame(t = decimal_date(mining_efficiency[, 4]),
                  y = 1 / mining_efficiency[, 3] * 1e6)

params_hash <- nls(log(y) ~ log(0.01 + (a - 0.01) / (1 + exp(-b * (t - c)))),
                   data = hash,
                   start = list( a = 42729, b = 1.25, c = 2019))


################################################################################
# Plot Data and Models
################################################################################

# General Processor Trends
plot(transitors[, 1], transitors[, 2], log = "y",
     xlab = "year", ylab = "Number of", ylim = c(10^(-1), 10^8))
points(specint[, 1], specint[, 2], pch = 2, col = 2)
points(frequency[, 1], frequency[, 2], pch = 3, col = 3)
points(watts[, 1], watts[, 2], pch = 4, col = 4)
points(cores[, 1], cores[, 2], pch = 5, col = 6)
lines(specint[, 1], sigmoid(coef(params_spec), specint[, 1]), col = 2)
lines(frequency[, 1], sigmoid(coef(params_freq), frequency[, 1]), col = 3)
lines(watts[, 1], sigmoid(coef(params_watt), watts[, 1]), col = 4)
lines(cores[, 1], sigmoid(coef(params_core), cores[, 1]) + 1, col = 6)

legend("topleft",
       legend = c("Transistors (thousands)", "Single Thread Performance",
                  "Frequency (MHz)", "Typical Power (W)", "Logical Cores"),
       col = c(1,2,3,4,6), pch = c(1,2,3,4,5), lty = c(1,1,1,1,1), ncol = 1)

# Bitcoin Mining Efficiency
hash$fit <- gen_logistic(c(0.01,coef(params_hash)), hash[, 1])
plot_miner_eff <- ggplot(hash, aes(x = t, y = y)) +
  geom_jitter(color = "#7895aa") +
  geom_line(aes(y = fit),
            na.rm = TRUE,
            size = 1,
            color = "#004c6d") +
  scale_y_continuous(trans = 'log10',
                     labels =
                       scales::trans_format('log10',
                                            scales::math_format(10^.x))) +
  labs(title = "Bitcoin Miner Performance",
       subtitle = TeX("Fitted Model: $\\nu(t)^{-1} = 0.01\\left[MH/J\\right] + \\frac{42730\\left[MH/J\\right]}{1 + e^{-1.258\\left[1/a\\right](t - 2019.2\\left[a\\right])}}$"),
       color = "Legend",
       x = "year",
       y = "Miner Performance [MH/J]")
plot_miner_eff

# General Processor Efficiency
t <- seq(2001, length.out = 22)
processor_eff <- data.frame(t)
for (i in t) {
  processor_eff$eff[i == t] <- sigmoid(coef(params_spec), i) * (sigmoid(coef(params_core), i) + 1) /
    sigmoid(coef(params_watt), i)
}
plot_processor_eff <- ggplot(processor_eff, aes(x = t, y = eff)) +
  geom_line(size = 1, color = "#004c6d") +
  labs(title = "Processor Efficiency Trends",
       subtitle = "2001 to 2022",
       x = "year",
       y = "Processor Performance [SPECINT/J]")
plot_processor_eff

