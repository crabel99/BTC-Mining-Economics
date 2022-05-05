# Estimate the dimensionless Plank constant for the bitcoin network
# \Delta \frac{1}{t} \Delta t \geq \frac{1}{4\pi}

library(dplyr)
library(ggplot2)


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
  # Dimensionless hashes
  eval_set$hash_rate <- mean(eval_set$delta_t) / eval_set$delta_t
  diff_periods[i, 1:3] <- c(i, var(eval_set$delta_t), var(eval_set$hash_rate))
  fit <- fitdist(eval_set$delta_t,
                 distr = "gamma",
                 method = "mle")
  diff_periods[i, 4:6] <- c(2*sqrt(diff_periods[i,2])*sqrt(diff_periods[i,3]),
                            fit$estimate)
                            # fit_inv$estimate)
}

rm(i, eval_set, fit) #, fit_inv)

plot(diff_periods$hbar,log = "y")

sigma_tsigma_H <- mean(diff_periods$hbar[41:length(diff_periods$hbar)])/2

plot_hbar <- ggplot(diff_periods, aes(x = period, y = hbar/2)) +
  geom_jitter(color = "#7895aa") +
  geom_hline(yintercept = sigma_tsigma_H,
            size = 1,
            color = "#004c6d") +

  scale_y_continuous(trans = 'log10',
                     labels =
                       scales::trans_format('log10',
                                            scales::math_format(10^.x))) +
  labs(title = "Mining Epoch Uncertainty Plot",
       subtitle = TeX(paste("Fitted Model: (Blocks ",
                            41*2016,
                            " - Present)  $\\sigma_t\\sigma_H =",
                            sigma_tsigma_H,
                            "$",
                            sep = "")),
       color = "Legend",
       x = "Mining Epoch",
       y = TeX("$\\sigma_t\\sigma_H \\left[\\textrm{s}\\right]$"))

plot_hbar

kd <- density(diff_periods$hbar[41:length(diff_periods$hbar)])
plot(kd)
