# Estimate the dimensionless Plank constant for the bitcoin network
# \Delta \frac{1}{t} \Delta t \geq \frac{1}{2}

library(dplyr)
library(invgamma)
library(fitdistrplus)


n_diff <- as.integer(nrow(headers)/2016)
diff_periods <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(diff_periods) <- c("period","var_t","var_hr","hbar","alpha_t","beta_t")

for (i in 1:n_diff) {
  eval_set <- headers[2016 * (i - 1) + 1:2016 * i, ] %>%
    dplyr::select(bits, delta_t)
  eval_set <- eval_set[!is.na(eval_set$delta_t),]
  if (i == 1) {
    eval_set <- eval_set[-1,] %>% subset(delta_t < 3 * 3600)
  }
  eval_set$hash_rate <- 1 / eval_set$delta_t
  diff_periods[i, 1:3] <- c(i, var(eval_set$delta_t), var(eval_set$hash_rate))
  fit <- fitdist(eval_set$delta_t,
                 distr = "gamma",
                 method = "mle")
  diff_periods[i, 4:6] <- c(sqrt(4*diff_periods[i,2]*diff_periods[i,3]),
                            fit$estimate)
                            # fit_inv$estimate)
}

rm(i, eval_set, fit) #, fit_inv)

plot(diff_periods$hbar,log = "y")
