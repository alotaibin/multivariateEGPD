library("ggplot2")
library("ggh4x")
library("dplyr")
library("reshape2")
library("NeuralEstimators")

theta_test <- readRDS("intermediates/theta_test.rds")
theta_test_NBE <- readRDS("intermediates/theta_test_pointestimates.rds")
theta_test_samples <- readRDS("intermediates/theta_test_samples.rds")
theta_test_NPE <- sapply(theta_test_samples, function(x) apply(x, 1, median))
theta_test_naive <- readRDS("intermediates/theta_test_naiveestimates.rds")

parameter_labels <- c(
  "kappa" = expression(kappa), 
  "sigma" = expression(sigma), 
  "xi" = expression(xi), 
  "thL" = expression(theta[L]), 
  "thU" = expression(theta[U]), 
  "thw" = expression(theta[omega])
)

# Plots: marginals for multiple cases
K <- 3
theta <- theta_test[, 1:K, drop = F]
theta_samples <- theta_test_samples[1:K]
parameter_names <- rownames(theta)
plotlist <- lapply(1:K, function(k) {
  samples <- theta_samples[[k]]  # Extract samples for parameter vector k
  samples_df <- as.data.frame(t(samples))  # Convert to data frame
  names(samples_df) <- parameter_names
  samples_df <- reshape2::melt(samples_df, variable.name = "parameter")  # Convert to long format
  true_params_df <- data.frame(parameter = parameter_names, true_value = theta[, k])
  NBE_estimates_df <- data.frame(parameter = parameter_names, estimate = theta_test_NBE[, k])
  NPE_estimates_df <- data.frame(parameter = parameter_names, estimate = theta_test_NPE[, k])
  Naive_estimates_df <- data.frame(parameter = parameter_names, estimate = theta_test_naive[, k])
  
  Naive_estimates_df <- dplyr::mutate_at(Naive_estimates_df, .vars = "parameter", .funs = factor, levels = names(parameter_labels), labels = parameter_labels)
  NBE_estimates_df   <- dplyr::mutate_at(NBE_estimates_df, .vars = "parameter", .funs = factor, levels = names(parameter_labels), labels = parameter_labels)
  NPE_estimates_df   <- dplyr::mutate_at(NPE_estimates_df, .vars = "parameter", .funs = factor, levels = names(parameter_labels), labels = parameter_labels)
  true_params_df     <- dplyr::mutate_at(true_params_df, .vars = "parameter", .funs = factor, levels = names(parameter_labels), labels = parameter_labels)
  samples_df         <- dplyr::mutate_at(samples_df, .vars = "parameter", .funs = factor, levels = names(parameter_labels), labels = parameter_labels)
  
  gg <- ggplot(samples_df, aes(x = value)) +
    geom_density(alpha = 0.3, colour = "#619CFF") +
    facet_wrap(~parameter, scales = "free", nrow = 1, labeller = label_parsed) +
    theme_bw() +
    labs(x = "", y = "Density", colour = "", fill = "") +
    geom_vline(data = true_params_df, aes(xintercept = true_value), linetype = "dashed",colour="black") +
    # geom_vline(data = NBE_estimates_df, aes(xintercept = estimate), colour = "#00BA38") +
    # geom_vline(data = NPE_estimates_df, aes(xintercept = estimate), colour = "#619CFF") +
    # geom_vline(data = Naive_estimates_df, aes(xintercept = estimate), colour="#F8766D") + 
    theme(
      strip.background = element_blank(),
      panel.grid.major = element_blank(),  # removes major gridlines
      panel.grid.minor = element_blank(),   # removes minor gridlines
      strip.text = element_text(size = 12)  # Increase font size
      )
  if (k > 1) {
    gg <- gg + theme(strip.text = element_blank())
  }
  if (k != 2) {
    gg <- gg + labs(y = "")
  } else if (k == 3) {
    gg <- gg + labs(x = "Value")
  }
  return(gg)
})
figure <- egg::ggarrange(plots = plotlist, nrow = K)
figure

ggsave(
  plot = figure, 
  file = "posterior_and_estimates.pdf", 
  path = "Figures", device = "pdf", width = 12, height = 7
)

# Plots: all test cases 
df.naive <- data.frame(
  estimator = "Likelihood-moment", 
  parameter = rep(rownames(theta_test), ncol(theta_test)), 
  truth = as.vector(theta_test), 
  estimate = as.vector(theta_test_naive)
  )
df.nbe <- data.frame(
  estimator = "NBE", 
  parameter = rep(rownames(theta_test), ncol(theta_test)), 
  truth = as.vector(theta_test), 
  estimate = as.vector(theta_test_NBE)
)
df.npe <- data.frame(
  estimator = "NPE", 
  parameter = rep(rownames(theta_test), ncol(theta_test)), 
  truth = as.vector(theta_test), 
  estimate = as.vector(theta_test_NPE)
)
df <- rbind(df.naive, df.nbe, df.npe)
df <- dplyr::mutate_at(df, .vars = "parameter", .funs = factor, levels = names(parameter_labels), labels = parameter_labels)

figure <- ggplot2::ggplot(df) + 
  ggplot2::geom_point(ggplot2::aes(x=truth, y = estimate, colour  = estimator), alpha = 0.75) + 
  ggplot2::geom_abline(colour = "black", linetype = "dashed") +
  ggh4x::facet_grid2(estimator~parameter, scales = "free", independent = "y", labeller = label_parsed) + 
  ggplot2::labs(colour = "") + 
  ggplot2::theme_bw() +
  ggplot2::theme(strip.background = ggplot2::element_blank()) + 
  theme(
    strip.text.y = element_blank(),
    strip.background = element_blank(),
    # panel.grid.major = element_blank(),  # removes major gridlines
    # panel.grid.minor = element_blank(),   # removes minor gridlines
    strip.text = element_text(size = 12)  # Increase font size
    ) #+ 
  # theme(legend.position = "top")
  # theme(legend.position = "none")

ggsave(
  plot = figure, 
  file = "assessment.pdf", 
  path = "Figures", device = "pdf", width = 11, height = 5
)


### RMSEs 
rmse_NBE <- sqrt(apply((theta_test_NBE-theta_test)^2,1,mean))
rmse_NPE <- sqrt(apply((theta_test_NPE-theta_test)^2,1,mean))
rmse_naive <- sqrt(apply((theta_test_naive-theta_test)^2,1,mean))
rbind(rmse_naive, rmse_NBE, rmse_NPE)

### Trimmeds RMSEs
trimmed_rmse <- function(errors, trim = 0.1) {
  n <- length(errors)
  k <- floor(n * trim / 2)  # trim equally on both sides
  sorted <- sort(errors)
  trimmed <- sorted[(k + 1):(n - k)]
  sqrt(mean(trimmed))
}
rmse_trimmed_naive <- apply((theta_test_naive - theta_test)^2, 1, trimmed_rmse)
rmse_trimmed_NBE <- apply((theta_test_NBE - theta_test)^2, 1, trimmed_rmse)
rmse_trimmed_NPE <- apply((theta_test_NPE - theta_test)^2, 1, trimmed_rmse)
rbind(rmse_trimmed_naive, rmse_trimmed_NBE, rmse_trimmed_NPE)