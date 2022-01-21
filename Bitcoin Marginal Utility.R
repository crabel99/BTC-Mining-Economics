# This file imports the necessary block header data into R for calculating
# bitcoin's marginal utility.
#
# The Bitcoin blockchain is parsed using Blockchain Postgres Import:
# https://github.com/blkchain/blkchain
#
# Bitcoin Price History: https://fred.stlouisfed.org/series/CBBTCUSD#0

################################################################################
# Libraries
################################################################################
library(dplyr)
library(lubridate)
library(invgamma)
library(fitdistrplus)
library(foreach)
library(doParallel)
library(zoo)
library(hrbrthemes)
library(ggplot2)
library(mice)
library(caTools)
library(Metrics)


################################################################################
# Utility Functions
################################################################################

target <- function(nbits) {
  # Compute the target from the nbits in the block header
  hex <- as.hexmode(nbits)
  shift <- strtoi(substr(hex, 1, 2), 16)
  base <- strtoi(substr(hex, 3, 8), 16)
  target <- base * 2 ^ (8 * (shift - 3))
  return(target)
}

difficulty <- function(nbits) {
  # Bitcoin mining difficulty
  MAX_TARGET <- target(486604799L)
  difficulty <- MAX_TARGET / target(nbits)
  return(max(difficulty, 1))
}

hash_rate <- function(nbits, delta_t) {
  # HR = D / T * 2 ^ 32, where D is the difficulty and T is the time between
  # blocks.
  # This function returns the hash rate as TH/s
  return(difficulty(nbits) / delta_t * 2 ^ 32 / 10 ^ 12)
}

sats_to_block <- function(height) {
  #Total bitcoin miner subsidy issued as a function of block height
  i <- floor(height / 210000)
  if (i <= 32) {
    total <- 0
    if ( i > 0) {
      for (j in 0:(i - 1)) {
        total <- total + floor(50e8 / 2 ^ j) * 210000
      }
    }
    total <- total + floor(50e8 / 2 ^ i) * (height + 1 - 210000 * i)
    return(total)
  }
  return(0)
}

gen_logistic = function(params, x) {
  # Generalized Logistic Function
  params[1] - (params[1] - params[2]) /
    (1 + exp(-params[3] * (x - params[4])))
}
################################################################################
# Gather the Data
################################################################################

# The connection is established using a DSN to the appropriate database. See:
# https://db.rstudio.com/best-practices/drivers/

conn <- DBI::dbConnect(odbc::odbc(), "Abe", timeout = 10)

latest_block <- DBI::dbGetQuery(conn,
  "SELECT height
  FROM blocks
  WHERE height IS NOT NULL AND
  orphan = false
  ORDER BY height DESC LIMIT 1") %>%
  as.integer()

headers <- DBI::dbGetQuery(conn,
  paste("
        SELECT
        	b.height AS height,
        	b.time AS time,
        	b.bits AS bits,
        	SUM(o.value) AS coinbase
        FROM blocks AS b
        LEFT JOIN block_txs AS bt ON b.id = bt.block_id
        LEFT JOIN txouts AS o ON bt.tx_id = o.tx_id",
        # "WHERE block_height < 200000 AND block_height > 195000 AND
        # "WHERE height < ", 30000, " AND
        "WHERE height < ", latest_block, " AND
          bt.n = 0 AND b.orphan = false
        GROUP BY bt.tx_id, b.height, b.time, b.bits
        ORDER BY b.height ASC"))



# Use the smallest data class possible to save on storage space
headers$height <- as.integer(headers$height)
headers$time <- as.integer(headers$time)
headers$bits <- as.integer(headers$bits)

DBI::dbDisconnect(conn = conn)
rm(conn)

################################################################################
# Mining Efficiency Model
################################################################################
mining_efficiency <- read.csv("data/mining hash rate.csv", comment.char = "#")
mining_efficiency[,4] <- as.Date(mining_efficiency[, 4], format = "%d/%m/%Y")

hash = data.frame(t = decimal_date(mining_efficiency[, 4]),
                  y = 1 / mining_efficiency[, 3] * 1e6)

params_hash <- nls(log(y) ~ log(0.01 + (a - 0.01) / (1 + exp(-b * (t - c)))),
                   data = hash,
                   start = list( a = 42729, b = 1.25, c = 2019))

plot(hash[,1], hash[,2], log = "y", xlab = "year", ylab = "MH/J")
lines(hash[,1], gen_logistic(c(0.01,coef(params_hash)), hash[, 1]), col = 1)


################################################################################
# Finalize Model Data and Estimates
################################################################################
headers[ , "total_sats"] <- numeric()
headers[ , "hash_rate"] <- numeric()
headers[ , "delta_t"] <- numeric()
headers[ , "marginal_utility"] <- numeric()
headers$year <- headers$time %>% as_datetime() %>% decimal_date()

# Intitialize parallelization
ncores = detectCores() - 3
cl = makeCluster(ncores)
registerDoParallel(cl)


headers <- foreach(
  block = iter(headers, by = "row"),
  .export = c("headers", "params_hash"),
  .combine = rbind) %dopar% {
    i <- which(headers$height == block$height)
    block$total_sats <- sats_to_block(block$height)

    if (i > 1) {
      # Median Time Past Rule time difference
      block$delta_t <- (block$time - median(
        headers$time[max((i - 11), 1):(i - 1)])) / min(6, i / 2)

      block$hash_rate <- hash_rate(block$bits, block$delta_t)
      block$marginal_utility <- block$hash_rate * 600 * 1e8 /
        gen_logistic(c(0.01,coef(params_hash)), block$year) / block$coinbase
    }

    block
  }

stopCluster(cl)

rm(cl, ncores)

# plot(headers$year, headers$hash_rate, log = "y",
#      xlab = "year", ylab = "MH/J")
#
# plot(headers$year, headers$marginal_utility, log = "y",
#      xlab = "year", ylab = "MJ/BTC")
#
# fit <- fitdist(headers$delta_t[-1:-11], distr = "gamma", method = "mle")
# fit <- fitdist(headers$marginal_utility[-1:-11], distr = "invgamma", method = "mle")
# summary(fit)
# plot(fit)

# Resolve any missing data
headers$marginal_utility[is.infinite(headers$marginal_utility)] <- NA
headers$height[is.na(headers$marginal_utility)][-1]
temp <- headers[,c(1,8)]
temp <- temp[-1,] # remove the first NA
tempData <- mice(temp, m = 5, maxit = 50, meth = 'pmm', seed = 500)
temp <- mice::complete(tempData)
headers$marginal_utility[-1] <- temp$marginal_utility

rm(temp, tempData)

# Calculate the moving averages, these need to be odd numbers so that the
# averages are symmetric. For the daily average use 143. For the bi-weekly,
# use the off-by-one bug value of 2015. And for the monthly, use 1/12th of an
# 8766-hour year.
headers <- headers %>%
  dplyr::mutate(util_01da = zoo::rollmean(marginal_utility, k = 143, fill = NA),
                util_14da = zoo::rollmean(marginal_utility,
                                          k = 2015,
                                          fill = NA),
                util_30da = zoo::rollmean(marginal_utility,
                                          k = 4383,
                                          fill = NA))



plot_util <- ggplot(headers, aes(x = year, y = marginal_utility)) +
  geom_density_2d(na.rm = TRUE, size = 0.25, color = "blue") +
  geom_line(mapping = aes(x = year, y = util_14da),
            na.rm = TRUE,
            size = 1,
            color = "black") +
  geom_vline(mapping = aes(xintercept = headers[210001,c("year")]),
             color = "red") +
  geom_vline(mapping = aes(xintercept = headers[420001,c("year")]),
             color = "red") +
  geom_vline(mapping = aes(xintercept = headers[630001,c("year")]),
             color = "red") +
  scale_y_continuous(trans = 'log10',
                     breaks = scales::trans_breaks('log10', function(x) 10^x),
                     labels =
                       scales::trans_format('log10',
                                            scales::math_format(10^.x))) +
  scale_color_discrete(name = "Legend",
                       # breaks = c("marginal_utility", "util_14da", "trt2"),
                       labels = c("Histogram",
                                "14-day MA",
                                "Halvenings")) +
  labs(title = "Marginal Utility of Bitcoin [MJ/BTC]",
       subtitle = "All Time: histogram with 14-day moving average",
       color = "Legend",
       x = "Year",
       y = "Marginal Utility [MJ/BTC]")


# Download Price History
CBBTCUSD <- readxl::read_excel("data/CBBTCUSD.xls", skip = 10)
CBBTCUSD$observation_date <- as.Date(CBBTCUSD$observation_date,
                                     format = "%Y-%m-%d")
CBBTCUSD$year <- decimal_date(CBBTCUSD$observation_date)
CBBTCUSD <- subset(CBBTCUSD, select = -observation_date)
CBBTCUSD <- CBBTCUSD[!(CBBTCUSD$CBBTCUSD %in% 0.00),]
CBBTCUSD$marginal_utility <- NA

closest <- function(vect,val){
  vect[which(abs(vect - val) == min(abs(vect - val)))]
}
for (i in 1:length(CBBTCUSD$year)) {
   CBBTCUSD$marginal_utility[i] <-
     headers$util_01da[headers$year %in%
                                closest(headers$year, CBBTCUSD$year[i])]
}
CBBTCUSD <- CBBTCUSD[!is.na(CBBTCUSD$marginal_utility), ]
rm(i)

# Fit the data to a power law relationship and plot
data_input <- subset(CBBTCUSD, select = -year)
split <- sample.split(data_input$marginal_utility, SplitRatio = 0.8)

# Linear regression of the log of the data
train <- subset(data_input, split == "TRUE") %>% log()
test <- subset(data_input, split == "FALSE")  %>% log()
fit <- lm(CBBTCUSD ~ marginal_utility, train)
pred <- predict(fit, test, se.fit = TRUE)
res <- resid(fit)
fit_norm <- fitdist(res, distr = "norm", method = "mle")

# produce a residual vs fitted plot for visualizing heteroscedasticity
#produce residual vs. fitted plot
plot(fitted(fit), res)
#add a horizontal line at 0
abline(0,0)

plot(fit_norm)




Y_test <- test$CBBTCUSD
sum_squared_error <- sse(Y_test, pred$fit)
mean_squared_error <- mse(Y_test, pred$fit)
root_mean_squared_error <- rmse(Y_test, pred$fit)
relative_squared_error <- rse(Y_test, pred$fit)
mean_abs_error <- mae(Y_test, pred$fit)
mean_abs_deviation <- mean_abs_error/length(Y_test)
relative_abs_error <- rae(Y_test, pred$fit)

error <- Y_test - pred$fit
R2 <- 1 - sum(error^2)/sum((Y_test - mean(Y_test))^2)
Adj_R2 <- 1 - (mean_squared_error/var(Y_test))

fitted_data <- data_input
fitted_data$CBBTCUSD <- exp(fit$coefficients[1]) *
  fitted_data$marginal_utility ^ fit$coefficients[2]

est_price <- exp(fit$coefficients[1]) *
  tail(headers$util_14da[!is.na(headers$util_14da)],
       n = 1) ^ fit$coefficients[2]
act_price <- 38000
dif_price <- log(act_price) - log(est_price)
dif_sigma <- dif_price/fit_norm$estimate[2]

plot_BTCUSD <- ggplot(CBBTCUSD, aes(x = marginal_utility, y = CBBTCUSD)) +
  geom_jitter() +
  geom_line(data = fitted_data,
            na.rm = TRUE,
            size = 1,
            color = "red") +
  scale_x_continuous(trans = 'log10') +
  scale_y_continuous(trans = 'log10') +
  labs(title = "Bitcoin Price to Marginal Utility",
       subtitle = "Dec 2014 - Present",
       color = "Legend",
       x = "Marginal Utility [MJ/BTC]",
       y = "Price [USD/BTC]")

# Outputs
plot_util
plot_BTCUSD
