# Boolean indicating whether (TRUE) or not (FALSE) to quickly establish that the code is working properly
quick <- identical(commandArgs(trailingOnly = TRUE)[1], "--quick")

library("NeuralEstimators")
library("JuliaConnectoR")
library("dplyr")
library("parallel")
library("ggplot2")
source("R/simulations/ModelSimulator.R")
# Helper functions for neural posterior estimation
sampleposterior <- juliaFun("sampleposterior")
logdensity <- juliaFun("logdensity")
Sys.setenv("JULIACONNECTOR_JULIAOPTS" = "--project=.") 
juliaEval('using NeuralEstimators, Flux, CUDA')
source("R/simulations/Architecture.R")

int_path <- "intermediates"
dir.create(int_path, recursive = TRUE, showWarnings = FALSE)
theta_scenarios <- readRDS("intermediates/theta_scenarios.rds")

# Training hyperparameters
epochs <- ifelse(quick, 10, 100)

# Sampler from the prior
sampler <- function(K) {
  kappa <- runif(K, 0.1, 10)
  sigma <- runif(K, 0.1, 3)
  xi    <- runif(K, 0.01, 0.5)
  thL   <- runif(K, 0.1, 20)
  thU   <- runif(K, 0.1, 20)
  thw   <- runif(K, 0.01, 0.49)
  theta <- matrix(c(kappa, sigma, xi, thL, thU, thw), byrow = TRUE, ncol = K)
  return(theta)
}

# Variance stabilizing transformation applied to the data 
signed_log <- function(x) (sign(x) * log1p(abs(x))) - 1

# Data simulator 
simulator <- function(Theta, m, mc.cores = detectCores() - 1) {
  Z <- mclapply(seq_len(ncol(Theta)), function(k) {
    m_k <- ifelse(length(m) > 1, sample(m, 1), m)
    theta_k <- Theta[, k]
    z <- t(rsim(m_k, theta_k)$Y)
    return(z)
  }, mc.cores = mc.cores)
  return(Z)
}

# Training/validation/test sets
set.seed(1)
K <- if(quick) 20000 else 100000
m <- if(quick) 500:1000 else 1000:4000
theta_train <- sampler(K)
theta_val   <- sampler(K/10)
theta_test  <- sampler(1000)
theta_test <- cbind(theta_scenarios, theta_test)
Z_train <- simulator(theta_train, m)
Z_val   <- simulator(theta_val, m)
Z_test  <- simulator(theta_test, max(m))
Z_test  %>% saveRDS(file = file.path(int_path, "Z_test.rds"))
theta_test  %>% saveRDS(file = file.path(int_path, "theta_test.rds"))

# Apply variance stabilizer to data 
Z_train <- lapply(Z_train, signed_log)
Z_val   <- lapply(Z_val, signed_log)
Z_test  <- lapply(Z_test, signed_log)

# Apply variance stabilizer to parameters 
theta_train <- log(theta_train)
theta_val   <- log(theta_val)
theta_test  <- log(theta_test)

# Training stage: NPE
NPE <- train(
  NPE, 
  theta_train = theta_train,
  theta_val = theta_val,
  Z_train = Z_train,
  Z_val = Z_val,
  epochs = epochs,
  stopping_epochs = 10,
  savepath = file.path(int_path, "NPE")
)
savestate(NPE, file.path(int_path, "NPE.bson"))

# Training stage: NBE
NBE <- train(
  NBE, 
  theta_train = theta_train,
  theta_val = theta_val,
  Z_train = Z_train,
  Z_val = Z_val,
  stopping_epochs = 10,
  epochs = epochs,  
  savepath = file.path(int_path, "NBE")
)
savestate(NBE, file.path(int_path, "NBE.bson"))

# Load the trained estimators
loadstate(NBE, file.path(int_path, "NBE.bson"))
loadstate(NPE, file.path(int_path, "NPE.bson"))

# Assessment with the test set
parameter_names <- rownames(theta_test)
assessment <- assess(list(NBE, NPE), theta_test, Z_test, estimator_names = c("NBE", "NPE"), parameter_names = parameter_names)
assessment$estimates$truth <- exp(assessment$estimates$truth)
assessment$estimates$estimate <- exp(assessment$estimates$estimate)

## Save figures into intermediates/Figures/Simulation
fig_path <- file.path("intermediates", "Figures", "Simulation")
dir.create(fig_path, recursive = TRUE, showWarnings = FALSE)

ggsave(
  plot = plotestimates(assessment), 
  file = "Neural_assessment1.pdf", 
  path = fig_path, device = "pdf", width = 18, height = 4
)

assessment_NBE <- assess(NBE, theta_test, Z_test, estimator_names = "NBE", parameter_names = parameter_names)
assessment_NPE <- assess(NPE, theta_test, Z_test, estimator_names = "NPE", parameter_names = parameter_names)
assessment_NBE$estimates$estimate <- exp(assessment_NBE$estimates$estimate)
assessment_NBE$estimates$truth <- exp(assessment_NBE$estimates$truth)
assessment_NPE$estimates$estimate <- exp(assessment_NPE$estimates$estimate)
assessment_NPE$estimates$truth <- exp(assessment_NPE$estimates$truth)

ggsave(
  plot = ggpubr::ggarrange(plotestimates(assessment_NBE), plotestimates(assessment_NPE), nrow = 2), 
  file = "Neural_assessment2.pdf", 
  path = fig_path, device = "pdf", width = 18, height = 6
)
 
# Inference over the test set (NBE)
inference_time <- system.time({
  thetahat <- exp(estimate(NBE, Z_test))
})
writeLines(as.character(inference_time["elapsed"]), file.path(int_path, "inference_time_NBE.txt"))
thetahat %>% saveRDS(file = file.path(int_path, "theta_test_pointestimates.rds"))

# Inference over the test set (NPE)
inference_time <- system.time({
  samples <- sampleposterior(NPE, Z_test)
  samples <- juliaLet('broadcast.(exp, samples)', samples = samples)
})
writeLines(as.character(inference_time["elapsed"]), file.path(int_path, "inference_time_NPE.txt"))
samples <- juliaGet(samples)
samples <- lapply(samples, function(x) {
  attributes(x)$JLTYPE <- NULL
  return(x)
})
saveRDS(samples, file = file.path(int_path, "theta_test_samples.rds"))
logdensities <- c(logdensity(NPE, theta_test, Z_test))
saveRDS(logdensities, file = file.path(int_path, "theta_test_logdensities.rds"))
