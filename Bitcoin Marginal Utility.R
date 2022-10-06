# This file imports the necessary block header data into R for calculating
# bitcoin's marginal utility.
#
# The Bitcoin blockchain is parsed using Blockchain Postgres Import:
# https://github.com/blkchain/blkchain
#
#
# Thank you to @mutatrum for the BTCUSD API endpoint.
# Bitcoin Price History: https://community-api.coinmetrics.io/
# See:
# https://docs.coinmetrics.io/exchanges/all-exchanges
# https://docs.coinmetrics.io/methodologies/reference-rates/hourly-reference-rates-methodology
# https://docs.coinmetrics.io/api#timestamps
#
#
# The general approach to the model is the following:
# I. Estimate bitcoin's marginal utility
#     A. Download the block data
#     B. Model the miner efficiency
#     C. Compute the marginal utility
#     D. Compute the 1-day, 14-day, and 30-day moving averages
# II. Model Price Movements as a Function of Marginal Utility
#     A. Download price data
#     B. Estimate model parameters
#     c. Evaluate the predictive capacity of the model
#
# The price movements use the 1-day moving average of the marginal utility as
# the input parameter. The price data should be the daily moving average.
# The price used from FRED are not daily averages. Once a daily average price
# history is found, adjust the search time by 12-hr so that the moving average
# corresponds to the daily average. Make sure that the price history is in UTC.
#
# The included plotted data is up to block height 719186.

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
library(latex2exp)

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

closest <- function(vect,val) {
  vect[which(abs(vect - val) == min(abs(vect - val)))]
}
################################################################################
# Gather the Data
################################################################################

# External price history time series
api_url <- "https://community-api.coinmetrics.io/"
api_endpint <- "v4/timeseries/asset-metrics"
api_query <- "?page_size=10000&metrics=PriceUSD&assets=btc"
api_request <- paste(api_url, api_endpint, api_query, sep = "")
api_response <- httr::GET(url = api_request)
httr::http_status(api_response)
api_content <- httr::content(api_response)
api_data <- bind_rows(api_content)

# The 12 hour shift is so that the daily moving average of marginal utility
# corresponds to the average for the day.
api_data$year <- api_data$time %>% ymd_hms() %>% decimal_date() + 12/8766
api_data$PriceUSD <- api_data$PriceUSD %>% as.numeric()
api_data <- api_data %>% subset(select = -asset) %>% subset(select = -time)

rm(api_url, api_endpint, api_query, api_request, api_response, api_content)


# Local parsed blockchain
# The connection is established using a DSN to the appropriate database. See:
# https://db.rstudio.com/best-practices/drivers/

conn <- DBI::dbConnect(odbc::odbc(), "Abe", timeout = 10)

if (exists("latest_block")) {
  query_where_stmt <- paste("height >",latest_block, "AND")
} else {
  query_where_stmt <- ""
}

query_result <- DBI::dbGetQuery(conn,
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
        "WHERE ", query_where_stmt, "
          bt.n = 0 AND b.orphan = false
        GROUP BY bt.tx_id, b.height, b.time, b.bits
        ORDER BY b.height ASC"))

# Use the smallest data class possible to save on storage space
query_result$height <- as.integer(query_result$height)
query_result$time <- as.integer(query_result$time)
query_result$bits <- as.integer(query_result$bits)

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
# Allocate memory for values to be generated
query_result[ , "total_sats"] <- numeric()
query_result[ , "hash_rate"] <- numeric()
query_result[ , "delta_t"] <- numeric()
query_result[ , "marginal_utility"] <- numeric()
query_result$year <- query_result$time %>% as_datetime() %>% decimal_date()

# Update existing data
if (exists("latest_block")) {
  headers <- rbind(headers[, 1:9], query_result)
} else {
  headers <- query_result
}

# Initialize parallelization
ncores = detectCores() - 3
cl = makeCluster(ncores)
registerDoParallel(cl)

query_result <- foreach(
  block = iter(query_result, by = "row"),
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

headers[is.na(headers$total_sats),] <- query_result
latest_block <- tail(headers[, 1], n = 1)

rm(cl, ncores, query_result)

# Estimate the MTP Delta t for a given mining epoch
time_data <- headers[
  headers[
    closest(headers$time,1666298078) == headers$time, 3] == headers$bits, 7]
time_fit <- fitdist(time_data, distr = "gamma", method = "mle")
# fit <- fitdist(headers$marginal_utility[-1:-11], distr = "invgamma", method = "mle")
plot(time_fit)

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
  # stat_density_2d(aes(fill = factor(stat(level))), geom = "polygon") +
  geom_density_2d(na.rm = TRUE,
                  size = 0.25,
                  color = "#7895aa",
                  show.legend = TRUE) +
  geom_line(mapping = aes(x = year, y = util_14da),
            na.rm = TRUE,
            size = 1,
            color = "#004c6d") +
  geom_vline(mapping = aes(xintercept = headers[210001,c("year")]),
             color = "#ff6361") +
  geom_vline(mapping = aes(xintercept = headers[420001,c("year")]),
             color = "#ff6361") +
  geom_vline(mapping = aes(xintercept = headers[630001,c("year")]),
             color = "#ff6361") +
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
  labs(title = "Marginal Value of Bitcoin [MJ/BTC]",
       subtitle = "All Time: histogram with 14-day moving average",
       color = "Legend",
       x = "Year",
       y = "Marginal Value [MJ/BTC]")

################################################################################
# Model Price Movements as a Function of Marginal Utility
################################################################################

api_data$marginal_utility <- NA

for (i in 1:length(api_data$year)) {
   api_data$marginal_utility[i] <-
     headers$util_01da[headers$year %in%
                                closest(headers$year, api_data$year[i])]
}
api_data <- api_data[!is.na(api_data$marginal_utility), ]
rm(i)

#######################################
# Estimate Model Parameters
#######################################
# Fit the data to a power law relationship and plot
data_input <- subset(api_data, select = -year)
split <- sample.split(data_input$marginal_utility, SplitRatio = 0.8)

# Linear regression of the log of the data
train <- subset(data_input, split == "TRUE") %>% log()
test <- subset(data_input, split == "FALSE")  %>% log()
fit <- lm(PriceUSD ~ marginal_utility, train)
pred <- predict(fit, test, se.fit = TRUE)
res <- resid(fit)
fit_norm <- fitdist(res, distr = "norm", method = "mle")

# produce a residual vs fitted plot for visualizing heteroscedasticity
#produce residual vs. fitted plot
plot(fitted(fit), res)
#add a horizontal line at 0
abline(0,0)

plot(fit_norm)

# Model the bitcoin price using the power law model
btc_price <- function(marg_utl, const, exponent) {
  return(exp(const) * marg_utl ^ exponent)
}

#Evaluate the size of bitcoin price movements from the power law model
btc_num_sigma <- function(mu_price, act_price,sigma) {
  dif_price <- log(act_price) - log(mu_price)
  return(dif_price/sigma)
}

btc_marg_util <- function(year, params) {
  exp(params[1] + params[2]*year)
}

# Example of model using real world data
est_price <- btc_price(tail(headers$util_01da[!is.na(headers$util_01da)],
                            n = 1),
                       fit$coefficients[1],
                       fit$coefficients[2])
dif_sigma <- btc_num_sigma(est_price, 36000, fit_norm$estimate[2])

# Plot the btc price model vs the data
fitted_data <- data_input
fitted_data$PriceUSD <- exp(fit$coefficients[1]) *
  fitted_data$marginal_utility ^ fit$coefficients[2]

latest_point <- tail(api_data, n = 1)
latest_point <- data.frame("PriceUSD" = latest_point[1,1],
                           "marginal_utility" = latest_point[1,3])

plot_P_m <- headers %>% ggplot(aes(x = total_sats, y = marginal_utility)) +
  geom_jitter(color = "#7895aa") +
  # geom_density_2d(na.rm = TRUE,
  #                 size = 0.25,
  #                 color = "#7895aa",
  #                 show.legend = TRUE) +
  scale_x_continuous(trans = 'log10',
                     labels =
                       scales::trans_format('log10',
                                            scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10',
                     labels =
                       scales::trans_format('log10',
                                            scales::math_format(10^.x)))

ggsave(
  "BTC P-m.jpg",
  plot = plot_P_m,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "print"
)

plot_BTCUSD <- ggplot(api_data, aes(x = marginal_utility, y = PriceUSD)) +
  geom_jitter(color = "#7895aa") +
  geom_line(data = fitted_data,
            na.rm = TRUE,
            size = 1,
            color = "#004c6d") +
  geom_point(data = latest_point,
             size = 3,
             color = "red") +
  scale_x_continuous(trans = 'log10',
                     labels =
                       scales::trans_format('log10',
                                            scales::math_format(10^.x))) +
  scale_y_continuous(trans = 'log10',
                     labels =
                       scales::trans_format('log10',
                                            scales::math_format(10^.x))) +
  labs(title = "Bitcoin Price to Marginal Value",
       subtitle = TeX(paste("(18 July 2010 - Present) Fitted Model: $price =",
                            round(exp(fit$coefficients[1]), digits = 4),
                            " P^{",
                            round(fit$coefficients[2], digits = 4),
                            "}$",
                            sep = "")),
       color = "Legend",
       x = TeX("$P \\left[MJ/BTC\\right]$"),
       y = "Price [USD/BTC]")

#######################################
# Evaluate the Predictive Capacity of the Model
#######################################
Y_test <- test$PriceUSD
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
fit_norm_test <- fitdist(error, distr = "norm", method = "mle")
plot(fit_norm_test)

#######################################
# Projected Bitcoin Marginal Value
#######################################
proj_data <- headers[headers$year > 2014, ]
proj_data$year <- proj_data$year - 2014
proj_data <- subset(proj_data, select = c(height, year, marginal_utility))
proj_data$ln_util <- log(proj_data$marginal_utility)
proj_split <- sample.split(proj_data$marginal_utility, SplitRatio = 0.8)

# Linear regression of the log of the data
proj_train <- subset(proj_data, proj_split == "TRUE")
proj_test <- subset(proj_data, proj_split == "FALSE")
proj_fit <- lm(ln_util ~ year, proj_train)
proj_fit_norm <- fitdist(resid(proj_fit), distr = "norm", method = "mle")
proj_pred <- predict(proj_fit, proj_test, se.fit = TRUE)
# plot(proj_fit_norm)
proj_data$pred_util <- exp(proj_fit$coefficients[2]*proj_data$year +
                             proj_fit$coefficients[1])

Y_proj_test <- proj_test$ln_util
proj_mean_squared_error <- mse(Y_proj_test, proj_pred$fit)
proj_error <- Y_proj_test - proj_pred$fit
proj_R2 <- 1 - sum(proj_error^2)/sum((Y_proj_test - mean(Y_proj_test))^2)
proj_Adj_R2 <- 1 - (proj_mean_squared_error/var(Y_proj_test))
proj_fit_norm_test <- fitdist(proj_error, distr = "norm", method = "mle")
plot(proj_fit_norm_test)
proj_data$year <- proj_data$year + 2014

# Add new data series to plot
plot_util <- plot_util +
  geom_line(mapping = aes(x = year, y = pred_util),
            data = proj_data,
            na.rm = TRUE,
            size = 1.0,
            color = "#bc5090")

# Plot Outputs
plot_util
ggsave(
  "BTC Marginal Value.pdf",
  plot = plot_util,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "retina"
)
ggsave(
  "BTC Marginal Value.jpg",
  plot = plot_util,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "print"
)


plot_BTCUSD
ggsave(
  "BTC Marginal Value-Price.pdf",
  plot = plot_BTCUSD,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "retina"
)
ggsave(
  "BTC Marginal Value-Price.jpg",
  plot = plot_BTCUSD,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "print"
)

# Fit Summaries
print(summary(proj_fit))
print(summary(fit))
