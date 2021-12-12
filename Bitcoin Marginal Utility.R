# This file imports the necessary block header data into R for calculating
# bitcoin's marginal utility.
#
# The Bitcoin blockchain is parsed using a fork of Bitcoin-Abe:
# https://github.com/crabel99/bitcoin-abe/tree/python2-to-python3

################################################################################
# Libraries
################################################################################
library(dplyr)
library(lubridate)
library(invgamma)
library(fitdistrplus)
library(foreach)
library(doParallel)


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
  #Bitcoin Miner Subsidy as a function of block height
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

con <- DBI::dbConnect(odbc::odbc(), "Abe", timeout = 10)

latest_block <- DBI::dbGetQuery(con,
  "SELECT block_height
  FROM block
  WHERE block_height IS NOT NULL ORDER BY block_height DESC LIMIT 1") %>%
  as.integer()

headers <- DBI::dbGetQuery(con,
  paste("
        SELECT
          block_id,
          block_height,
          block_ntime,
          block_nbits,
          prev_block_id
        FROM block",
        # "WHERE block_height < 200000 AND block_height > 195000 ",
        "WHERE block_height < ", latest_block,
        "ORDER BY block_height ASC"))

txouts <- DBI::dbGetQuery(con,
  paste("
        SELECT
          block_id,
          in_longest,
          txout_pos,
          txout_value
        FROM txout_detail",
        # "WHERE block_height < 200000 AND block_height > 195000 ",
        "WHERE block_height < ", latest_block,
        "ORDER BY block_height, txout_pos ASC"))


# Use the smallest data class possible to save on storage space
headers$block_id <- as.integer(headers$block_id)
headers$block_height <- as.integer(headers$block_height)
headers$prev_block_id <- as.integer(headers$prev_block_id)
headers$block_nbits <- as.integer(headers$block_nbits)

txouts$block_id <- as.integer(txouts$block_id)
txouts$in_longest <- as.logical(txouts$in_longest)
txouts$txout_pos <- as.integer(txouts$txout_pos)

# Sort based upon longest chain and remove orphan/forks
n_blocks <- length(headers$block_id)
longest_id <- integer(n_blocks)
longest_id[1]  <- headers$block_id[
  !(headers$block_id %in% headers$prev_block_id)]
longest_id[1] <- which(headers$block_id == longest_id[1])

for (i in 2:n_blocks) {
  longest_id[i] <- which(
    headers$block_id == headers$prev_block_id[longest_id[i - 1]])
}

longest_id <- rev(longest_id)
headers <- headers[longest_id, ]

# rm(con, i, longest_id)


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
headers[ , "coinbase"] <- numeric()
headers[ , "hash_rate"] <- numeric()
headers[ , "delta_t"] <- numeric()
headers[ , "marginal_utility"] <- numeric()
headers$year <- headers$block_ntime %>% as_datetime() %>% decimal_date()

# Intitialize parallelization
ncores = detectCores() - 1
cl = makeCluster(ncores)
registerDoParallel(cl)


headers <- foreach(
  block = iter(headers, by = "row"),
  .export = c("headers", "txouts", "params_hash"),
  .combine = rbind) %dopar% {
    i <- which(headers$block_id == block$block_id)
    block$total_sats <- sats_to_block(block$block_height)
    block$coinbase <- sum(subset(
      txouts,txouts$block_id == block$block_id & txouts$txout_pos == 0)[,4])

    if (i > 1) {
      # Median Time Past Rule time difference
      block$delta_t <- (block$block_ntime - median(
        headers$block_ntime[max((i - 11), 1):(i - 1)])) / min(6, i / 2)

      block$hash_rate <- hash_rate(block$block_nbits, block$delta_t)
      block$marginal_utility <- block$hash_rate * 600 * 1e8 /
        gen_logistic(c(0.01,coef(params_hash)), block$year) / block$coinbase
    }

    block
  }

stopCluster(cl)

# plot(headers$year, headers$hash_rate, log = "y",
#      xlab = "year", ylab = "MH/J")
#
# plot(headers$year, headers$marginal_utility, log = "y",
#      xlab = "year", ylab = "MJ/BTC")
#
# fit <- fitdist(headers$marginal_utility[-1:-11], distr = "invgamma", method = "mle")
# summary(fit)
# plot(fit)

