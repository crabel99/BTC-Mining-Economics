# Estimate the dimensionless Plank constant for the bitcoin network
# \Delta \frac{1}{t} \Delta t \geq \frac{1}{4\pi}

library(dplyr)
library(ggplot2)
library(latex2exp)
library(fitdistrplus)
library(metRology)


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


df <- data.frame("x" = diff_periods$hbar[41:length(diff_periods$hbar)])
fit_h_bar <- fitdist(df$x,"t.scaled",
                     start = list(df = 3, mean = mean(df$x), sd = sd(df$x)))

plot(fit_h_bar)
print(summary(fit_h_bar))
params <- as.list(fit_h_bar$estimate)

plot_hbar <- ggplot(diff_periods, aes(x = period, y = hbar)) +
  geom_jitter(color = "#7895aa") +
  geom_hline(yintercept = params$m,
            size = 1,
            color = "#004c6d") +

  scale_y_continuous(trans = 'log10',
                     labels =
                       scales::trans_format('log10',
                                            scales::math_format(10^.x))) +
  labs(title = "Mining Difficulty Period Uncertainty Plot",
       subtitle = TeX(paste("Fitted Model: (Blocks ",
                            40*2016,
                            " - Present)  $2\\sigma_t\\sigma_H =",
                            round(params$m, 3),
                            "$",
                            sep = "")),
       color = "Legend",
       x = "Mining Difficulty Period",
       y = TeX("$2\\sigma_t\\sigma_H$"))

plot_hbar_fit <- ggplot(data = df, aes(x)) +
  geom_histogram(aes(y = stat(density)),
                 color = "black",
                 fill = "white",
                 bins = 70) +
  geom_density(color = "black", fill = "grey", alpha = .7) +
  stat_function(fun = dt.scaled,
                args = params, color = "red") +
  labs(title = "Mining Plank's Constant Distribution",
       subtitle = TeX(paste("Fitted Model: Block ",
                            40*2016,
                            " - Present)  $h \\sim T(\\mu,\\sigma,\\nu)=",
                            "T(",round(params$m, 0),
                            ",", round(params$s,1),
                            ",", round(params$df,1),
                            ")$",
                            sep = "")),
       color = "Legend",
       x = TeX("$2 \\sigma_t \\, \\sigma_H$"),
       y = "density")

plot_hbar
ggsave(
  "hbar.pdf",
  plot = plot_hbar,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "retina"
)
ggsave(
  "hbar.jpg",
  plot = plot_hbar,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "print"
)

plot_hbar_fit
ggsave(
  "hbar Fit.pdf",
  plot = plot_hbar_fit,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "retina"
)
ggsave(
  "hbar Fit.jpg",
  plot = plot_hbar_fit,
  path = "images/",
  scale = 1,
  width = 10,
  height = 5.625,
  units = "in",
  dpi = "print"
)
