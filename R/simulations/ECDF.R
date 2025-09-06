library("ggplot2")

parameter_labels <- c(
  "kappa" = expression(kappa), 
  "sigma" = expression(sigma), 
  "xi" = expression(xi), 
  "thL" = expression(theta[L]), 
  "thU" = expression(theta[U]), 
  "thw" = expression(theta[omega])
)

# See sections 1.2 and 1.3 of: https://betanalpha.github.io/assets/case_studies/principled_bayesian_workflow.html#1_Questioning_Authority
# https://hyunjimoon.github.io/SBC/articles/rank_visualizations.html
# devtools::install_github("hyunjimoon/SBC")
# library("SBC")
# ranks_matrix <- matrix(c(sample(1000), sample(1000), c(sample(500), sample(500))), ncol= 3, byrow = FALSE)
# SBC::plot_ecdf(ranks_matrix, max_rank = max(ranks_matrix))
# SBC::plot_ecdf_diff(ranks_matrix, max_rank = max(ranks_matrix))

compute_ranks_matrix <- function(posterior_samples, theta_true) {
  K <- length(posterior_samples)
  d <- nrow(theta_true)
  ranks_mat <- matrix(NA_integer_, nrow = d, ncol = K)
  
  for (k in seq_len(K)) {
    samples <- posterior_samples[[k]]  
    theta <- theta_true[, k]           
    
    for (j in seq_len(d)) {
      ranks_mat[j, k] <- sum(samples[j, ] < theta[j])
    }
  }
  
  rownames(ranks_mat) <- rownames(theta_true)
  
  return(ranks_mat)
}

# Compute simultaneous confidence intervals for one or more samples from the
# standard uniform distribution.
# N - sample length
# L - number of samples
# K - size of uniform partition defining the ECDF evaluation points.
# gamma - coverage parameter for the marginal distribution (binomial for
# one sample and hypergeometric for multiple rank transformed chains).
ecdf_intervals <- function(N, L, K, gamma) {
  lims <- list()
  z <- seq(0,1, length.out = K + 1)
  if (L == 1) {
    lims$lower <- qbinom(gamma / 2, N, z)
    lims$upper <- qbinom(1 - gamma / 2, N, z)
  } else {
    n <- N * (L - 1)
    k <- floor(z * L * N)
    lims$lower <- qhyper(gamma / 2, N, n, k)
    lims$upper <- qhyper(1 - gamma / 2, N, n, k)
  }
  lims$lower <- c(rep(lims$lower[1:K], each=2), lims$lower[K + 1])
  lims$upper <- c(rep(lims$upper[1:K], each=2), lims$upper[K + 1])
  lims
}

# Adjust coverage parameter to find simultaneous confidence intervals for the
# ECDF of a sample from the uniform distribution.
# N - length of samples
# K - number of equally spaced evaluation points, i.e. the right ends of the
# partition intervals.
adjust_gamma_optimize <- function(N, K, conf_level=0.95) {
  target <- function(gamma, conf_level, N, K) {
    z <- 1:(K - 1) / K
    z1 <- c(0,z)
    z2 <- c(z,1)
    
    # pre-compute quantiles and use symmetry for increased efficiency.
    x2_lower <- qbinom(gamma / 2, N, z2)
    x2_upper <- c(N - rev(x2_lower)[2:K], N)
    
    # Compute the total probability of trajectories inside the confidence
    # intervals. Initialize the set and corresponding probabilities known
    # to be 0 and 1 for the starting value z1 = 0.
    x1 <- 0
    p_int <- 1
    for (i in seq_along(z1)) {
      tmp <- p_interior(
        p_int, x1 = x1, x2 = x2_lower[i]: x2_upper[i],
        z1 = z1[i], z2 = z2[i], gamma = gamma, N = N
      )
      x1 <- tmp$x1
      p_int <- tmp$p_int
    }
    abs(conf_level - sum(p_int))
  }
  optimize(target, c(0, 1 - conf_level), conf_level, N = N, K = K)$minimum
}
p_interior <- function(p_int, x1, x2, z1, z2, gamma, N) {
  z_tilde <- (z2 - z1) / (1 - z1)
  N_tilde <- rep(N - x1, each = length(x2))
  p_int <- rep(p_int, each = length(x2))
  x_diff <- outer(x2, x1, "-")
  p_x2_int <- p_int * dbinom(x_diff, N_tilde, z_tilde)
  list(p_int = rowSums(p_x2_int), x1 = x2)
}

data_for_ecdf_plots <- function(ranks_matrix, prob = 0.99) {
  
  max_rank <- max(ranks_matrix)
  N <- nrow(ranks_matrix)
  K <- min(max_rank + 1, N)
  
  gamma <- adjust_gamma_optimize(N, K, conf_level = prob)
  
  z <- seq(0,1, length.out = K + 1)
  z_twice <- c(0, rep(z[2:(K + 1)], each = 2))
  
  limits_df <- as.data.frame(ecdf_intervals(
    N = N,
    L = 1,
    K = K,
    gamma = gamma))
  
  limits_df <- dplyr::mutate(limits_df,
                             x = z_twice,
                             lower = lower / N,
                             upper = upper / N,
                             # The uniform_val needs to be shifted w.r.t z_twice
                             uniform_val =  c(rep(z[1:K], each = 2), 1)
  )
  
    base_vals <- floor((0:K) * ((max_rank + 1) / K))
  ecdf_vals <- matrix(nrow = K + 1, ncol = ncol(ranks_matrix))
  colnames(ecdf_vals) <- colnames(ranks_matrix)
  for(i in 1:(K + 1)) {
    # Note: for pit calculations we would use (col + 1) / (max_rank + 1)
    # For ecdf we would use pit <= base_val / (max_rank + 1)
    # So the "+ 1" and "<=" can be subsumed in "<"
    ecdf_vals[i,] <- colMeans(ranks_matrix < base_vals[i])
  }
  
  
  ecdf_df <- as.data.frame(ecdf_vals)
  ecdf_df$..z <- z
  ecdf_df <- tidyr::pivot_longer(ecdf_df, -..z, names_to = "parameter", values_to = "ecdf")
  ecdf_df <- dplyr::rename(ecdf_df, z = ..z)
  
  list(limits_df = limits_df, ecdf_df = ecdf_df, K = K, N = N)
}


plot_ecdf <- function(theta, theta_samples, prob = 0.99, parameter_labels = NULL) {
  
  if (!is.null(parameter_labels)) {
    rownames(theta) <- names(parameter_labels)
  } 
  
  ranks_matrix <- compute_ranks_matrix(theta_samples, theta)
  ecdf_data <- data_for_ecdf_plots(t(ranks_matrix), prob = prob)
  
  N <- ecdf_data$N
  K <- ecdf_data$K
  
  ecdf_df <- dplyr::mutate(ecdf_data$ecdf_df, type = "sample ECDF")
  limits_df <- ecdf_data$limits_df
  limits_df$type <- "theoretical CDF"
  
  
  if (!is.null(parameter_labels)) {
    param_labeller <- ggplot2::label_parsed
    ecdf_df <- dplyr::mutate_at(ecdf_df, .vars = "parameter", .funs = factor, levels = names(parameter_labels), labels = parameter_labels)
  } else {
    param_labeller <- identity
  }
  
  # construct figure
  ggplot(ecdf_df, aes(color = type, fill = type)) +
    geom_ribbon(
      data = limits_df,
      aes(x = x, ymax = upper, ymin = lower),
      alpha = 0.33) +
    geom_step(
      aes(x = z, y = ecdf, group = parameter, alpha = 0.8)
    ) +
    scale_x_continuous(breaks = c(0, 0.5, 1), name = "Fractional rank of true parameter") +
    scale_y_continuous(breaks = c(0, 0.5, 1), name = "ECDF") +
    scale_color_manual(
      name = "",
      values = rlang::set_names(
        c("skyblue1", "black"),
        c("theoretical CDF", "sample ECDF")),
      # labels = c(
      #   "theoretical CDF" = expression(italic("theoretical CDF")),
      #   "sample ECDF" = expression(italic("sample ECDF"))
      # )
    ) +
    scale_fill_manual(
      name = "",
      values = c("theoretical CDF" = "skyblue",
                 "sample ECDF" = "transparent"),
      # labels = c(
      #   "theoretical CDF" = expression(italic("theoretical CDF")),
      #   "sample ECDF" = expression(italic("sample ECDF"))
      # )
    ) +
    scale_alpha_identity() +
    xlab(NULL) +
    ylab(NULL) +
    facet_wrap(~ parameter, nrow = 1, labeller = param_labeller)
}
  
theta_test <- readRDS("intermediates/theta_test.rds")
theta_test_samples <- readRDS("intermediates/theta_test_samples.rds")
figure <- plot_ecdf(theta_test, theta_test_samples, prob = 0.99, parameter_labels = parameter_labels) 
figure <- figure +
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12), 
    legend.position = "none"
  )

ggsave(
  plot = figure, 
  file = "ranks_ecdf.pdf", 
  path = "Figures", device = "pdf", width = 14, height = 2.5
)

figure_ecdf <- figure

# Posterior z-scores and contraction 
mu_post <- sapply(theta_test_samples, function(samples) rowMeans(samples))
sd_post <- sapply(theta_test_samples, function(samples) apply(samples, 1, sd))
z_score <- (mu_post - theta_test) / sd_post

var_prior <- apply(theta_test, 1, var)
var_post  <- sd_post^2
contraction <- 1 - var_post / var_prior
rownames(contraction) <- rownames(z_score)

df_z <- reshape2::melt(z_score, value.name = "z_score")
df_c <- reshape2::melt(contraction, value.name = "contraction")
df <- merge(df_z, df_c)
df <- dplyr::mutate_at(df, .vars = "Var1", .funs = factor, levels = names(parameter_labels), labels = parameter_labels)

#TODO could add true value as a colour (on a sequential, standardized colour scale showing relative values with respect to the value of the prior)
figure <- ggplot(df) + 
  geom_point(aes(x = contraction, y = z_score), alpha = 0.75) + 
  facet_wrap(Var1~., nrow = 1, labeller = label_parsed) + 
  labs(x = "Posterior contraction", y = "Posterior z-score") + 
  theme_bw() + 
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  )
ggsave(
  plot = figure, 
  file = "zscore_contraction.pdf", 
  path = "Figures", device = "pdf", width = 11, height = 2.5
)
figure_zscore <- figure

figure <- egg::ggarrange(
  figure_ecdf, 
  figure_zscore + theme(strip.text = element_blank()), 
  nrow = 2
  )

ggsave(
  plot = figure, 
  file = "ranks_ecdf_zscore_contraction.pdf", 
  path = "Figures", device = "pdf", width = 12, height = 5
)

# # ---- Toy examples ----
# 
# # Load packages
# library("NeuralEstimators")
# library("JuliaConnectoR")
# library("ggplot2")
# juliaEval('using NeuralEstimators, Flux, CUDA') 
# 
# # Sampling from the prior 
# # K: number of samples to draw from the prior
# sampler <- function(K) {
#   mu    <- rnorm(K)
#   sigma <- rgamma(K, 1)
#   Theta <- matrix(c(mu, sigma), byrow = TRUE, ncol = K)
#   return(Theta)
# }
# K <- 10000
# theta_train <- sampler(K)
# theta_val   <- sampler(K/10)
# 
# # Simulation from the statistical model
# # theta: a matrix of parameters drawn from the prior
# # m: number of conditionally independent replicates for each parameter vector
# simulate <- function(Theta, m) {
#   apply(Theta, 2, function(theta) t(rnorm(m, theta[1], theta[2])), simplify = FALSE)
# }
# m <- 100
# Z_train <- simulate(theta_train, m)
# Z_val   <- simulate(theta_val, m)
# 
# # Initialise the estimator
# estimator <- juliaEval('
#   n = 1    # dimension of each data replicate (univariate)
#   d = 2    # dimension of the parameter vector θ
#   w = 128  # width of each hidden layer 
#   
#   # Distribution used to approximate the posterior 
#   q = NormalisingFlow(d, d) 
# 
#   # Neural network (outputs d summary statistics)
#   psi = Chain(Dense(n, w, relu), Dense(w, w, relu))    
#   phi = Chain(Dense(w, w, relu), Dense(w, d))          
#   network = DeepSet(psi, phi)
# 
#   # Wrap the approximating distribution and neural network as a PosteriorEstimator
#   estimator = PosteriorEstimator(q, network)
# ')
# 
# # Train the estimator 
# estimator <- train(estimator,   
#                    theta_train = theta_train,
#                    theta_val   = theta_val,
#                    Z_train = Z_train,
#                    Z_val   = Z_val)
# 
# # Assess the estimator 
# theta_test <- sampler(1000)
# Z_test     <- simulate(theta_test, m)
# assessment <- assess(estimator, theta_test, Z_test)
# plotestimates(assessment)
# 
# theta_test_samples <- sampleposterior(estimator, Z_test)
# ranks_matrix <- compute_ranks_matrix(theta_test_samples, theta_test)
# plot_ecdf(ranks_matrix, max_rank = max(ranks_matrix))
# plot_ecdf_diff(ranks_matrix, max_rank = max(ranks_matrix))
# 
