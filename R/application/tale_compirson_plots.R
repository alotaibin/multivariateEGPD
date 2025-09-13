################################################################################
# Tail Comparison QQ Plots: Stations 11 & 23
################################################################################

# Required packages
library(mev)        # fit.gpd(), fbvpot()
library(evd)        # pgev(), rbvevd()
library(ggplot2)
library(patchwork)  # shared-legend layouts
library(RColorBrewer)

################################################################################
# Helper: Inverse GPD Transformation (used for mapping simulated data)
################################################################################
Fgpd_inv <- function(par, data, simdata, threshold_data, threshold_simdata) {
  sigma <- par[1]
  xi    <- par[2]
  L     <- threshold_simdata 
  L_data<- threshold_data
  n     <- length(simdata)

  # Proportion of simulated exceedances
  delta_n <- sum(simdata >= L) / n
  trans_gpd <- numeric(n)

  for (i in seq_len(n)) {
    if (simdata[i] >= L) {
      p_value <- 1 - (1 - simdata[i]) / delta_n
      p_value <- max(min(p_value, 0.9999), 1e-6)  # clip to (0,1)
      trans_gpd[i] <- L_data + qgpd(p_value, loc = 0, scale = sigma, shape = xi)
    } else {
      # Below threshold: empirical quantile from observed data
      trans_gpd[i] <- as.numeric(quantile(data, probs = simdata[i], na.rm = TRUE))
    }
  }
  return(trans_gpd)
}

################################################################################
# Load rainfall data (Stations 11 & 23)
################################################################################
df_11_23_scaled <- read.csv("R/data/pair_11_23_matrix_scaled.csv",
                            stringsAsFactors = FALSE)
mat_11_23_scaled <- as.matrix(df_11_23_scaled[,-1])  # drop date column if present
dim(mat_11_23_scaled)  # (2 × 2166)
real1 <- mat_11_23_scaled[1, ]
real2 <- mat_11_23_scaled[2, ]
bivariate <- data.frame(real1, real2)

################################################################################
# Fit marginal GPDs (95% thresholds for upper tail, 90% for lower tail)
################################################################################
thresh_up <- c(quantile(real1, 0.95, na.rm = TRUE),
               quantile(real2, 0.95, na.rm = TRUE))
thresh_low <- c(quantile(-real1, 0.90, na.rm = TRUE),
                quantile(-real2, 0.90, na.rm = TRUE))

fit_up1  <- fit.gpd(real1, thresh_up[1])
fit_up2  <- fit.gpd(real2, thresh_up[2])
fit_low1 <- fit.gpd(-real1, thresh_low[1])
fit_low2 <- fit.gpd(-real2, thresh_low[2])

################################################################################
# Fit bivariate POT models (log, bilog, HR, CT) for upper and lower tails
################################################################################
# Upper
fit_log_up   <- fbvpot(bivariate, threshold = thresh_up, model = "log",
                       scale1 = fit_up1$param[1], shape1 = fit_up1$param[2],
                       scale2 = fit_up2$param[1], shape2 = fit_up2$param[2], std.err = FALSE)
fit_bilog_up <- fbvpot(bivariate, threshold = thresh_up, model = "bilog",
                       scale1 = fit_up1$param[1], shape1 = fit_up1$param[2],
                       scale2 = fit_up2$param[1], shape2 = fit_up2$param[2], std.err = FALSE)
fit_hr_up    <- fbvpot(bivariate, threshold = thresh_up, model = "hr",
                       scale1 = fit_up1$param[1], shape1 = fit_up1$param[2],
                       scale2 = fit_up2$param[1], shape2 = fit_up2$param[2], std.err = FALSE)
fit_ct_up    <- fbvpot(bivariate, threshold = thresh_up, model = "ct",
                       scale1 = fit_up1$param[1], shape1 = fit_up1$param[2],
                       scale2 = fit_up2$param[1], shape2 = fit_up2$param[2], std.err = FALSE)

# Lower (negated data for lower tail)
biv_low <- data.frame(-real1, -real2)
fit_log_low   <- fbvpot(biv_low, threshold = thresh_low, model = "log",
                        scale1 = fit_low1$param[1], shape1 = fit_low1$param[2] + 0.2,
                        scale2 = fit_low2$param[1], shape2 = fit_low2$param[2] + 0.2, std.err = FALSE)
fit_bilog_low <- fbvpot(biv_low, threshold = thresh_low, model = "bilog",
                        scale1 = fit_low1$param[1], shape1 = fit_low1$param[2] + 0.2,
                        scale2 = fit_low2$param[1], shape2 = fit_low2$param[2] + 0.2, std.err = FALSE)
fit_hr_low    <- fbvpot(biv_low, threshold = thresh_low, model = "hr",
                        scale1 = fit_low1$param[1], shape1 = fit_low1$param[2] + 0.2,
                        scale2 = fit_low2$param[1], shape2 = fit_low2$param[2] + 0.2, std.err = FALSE)
fit_ct_low    <- fbvpot(biv_low, threshold = thresh_low, model = "ct",
                        scale1 = fit_low1$param[1], shape1 = fit_low1$param[2] + 0.2,
                        scale2 = fit_low2$param[1], shape2 = fit_low2$param[2] + 0.2, std.err = FALSE)

################################################################################
# Simulate from fitted models and transform back to GPD margins
################################################################################
m_data <- length(real1)

# Upper simulations
sim_log_up   <- rbvevd(m_data, dep   = fit_log_up$estimate, model = "log")
sim_bilog_up <- rbvevd(m_data, alpha = fit_bilog_up$estimate[1], beta = fit_bilog_up$estimate[2], model = "bilog")
sim_hr_up    <- rbvevd(m_data, dep   = fit_hr_up$estimate, model = "hr")
sim_ct_up    <- rbvevd(m_data, alpha = fit_ct_up$estimate[1], beta = fit_ct_up$estimate[2], model = "ct")

# Transform via Uniform → GPD margins
cdf_log_up1   <- pgev(sim_log_up[,1]);   cdf_log_up2   <- pgev(sim_log_up[,2])
cdf_bilog_up1 <- pgev(sim_bilog_up[,1]); cdf_bilog_up2 <- pgev(sim_bilog_up[,2])
cdf_hr_up1    <- pgev(sim_hr_up[,1]);    cdf_hr_up2    <- pgev(sim_hr_up[,2])
cdf_ct_up1    <- pgev(sim_ct_up[,1]);    cdf_ct_up2    <- pgev(sim_ct_up[,2])

gpd_log_up   <- cbind(
  Fgpd_inv(fit_up1$param, real1, cdf_log_up1,   thresh_up[1], 0.95),
  Fgpd_inv(fit_up2$param, real2, cdf_log_up2,   thresh_up[2], 0.95)
)
gpd_bilog_up <- cbind(
  Fgpd_inv(fit_up1$param, real1, cdf_bilog_up1, thresh_up[1], 0.95),
  Fgpd_inv(fit_up2$param, real2, cdf_bilog_up2, thresh_up[2], 0.95)
)
gpd_hr_up    <- cbind(
  Fgpd_inv(fit_up1$param, real1, cdf_hr_up1,    thresh_up[1], 0.95),
  Fgpd_inv(fit_up2$param, real2, cdf_hr_up2,    thresh_up[2], 0.95)
)
gpd_ct_up    <- cbind(
  Fgpd_inv(fit_up1$param, real1, cdf_ct_up1,    thresh_up[1], 0.95),
  Fgpd_inv(fit_up2$param, real2, cdf_ct_up2,    thresh_up[2], 0.95)
)

################################################################################
# Assemble observed + simulated datasets (upper & lower tails, incl. MEGPD)
################################################################################
observed_data <- cbind(real1, real2)

## Simulate data under estimated parameters
n2 <- ncol(mat_11_23_scaled) # 2166 obs
sample22 <- rXY(n2,
                estimates[1:3],
                c(estimates[4],estimates[4]),
                c(estimates[5],estimates[5]),
                estimates[6],
                PLOT=FALSE, onlyXY=FALSE)

ssX2 <- sample22$XY[,1]
ssY2 <- sample22$XY[,2]

sim_data_list_up <- list(
  Log   = gpd_log_up,
  Bilog = gpd_bilog_up,
  HR    = gpd_hr_up,
  CT    = gpd_ct_up,
  MEGPD = cbind(ssX2, ssY2)   # neural model samples
)

sim_data_list_low <- list(
  Log   = -cbind(Fgpd_inv(fit_low1$param, -real1, pgev(rbvevd(m_data, dep = fit_log_low$estimate, model = "log")[,1]), thresh_low[1], 0.90),
                 Fgpd_inv(fit_low2$param, -real2, pgev(rbvevd(m_data, dep = fit_log_low$estimate, model = "log")[,2]), thresh_low[2], 0.90)),
  Bilog = -cbind(Fgpd_inv(fit_low1$param, -real1, pgev(rbvevd(m_data, alpha = fit_bilog_low$estimate[1], beta = fit_bilog_low$estimate[2], model = "bilog")[,1]), thresh_low[1], 0.90),
                 Fgpd_inv(fit_low2$param, -real2, pgev(rbvevd(m_data, alpha = fit_bilog_low$estimate[1], beta = fit_bilog_low$estimate[2], model = "bilog")[,2]), thresh_low[2], 0.90)),
  HR    = -cbind(Fgpd_inv(fit_low1$param, -real1, pgev(rbvevd(m_data, dep = fit_hr_low$estimate, model = "hr")[,1]), thresh_low[1], 0.90),
                 Fgpd_inv(fit_low2$param, -real2, pgev(rbvevd(m_data, dep = fit_hr_low$estimate, model = "hr")[,2]), thresh_low[2], 0.90)),
  CT    = -cbind(Fgpd_inv(fit_low1$param, -real1, pgev(rbvevd(m_data, alpha = fit_ct_low$estimate[1], beta = fit_ct_low$estimate[2], model = "ct")[,1]), thresh_low[1], 0.90),
                 Fgpd_inv(fit_low2$param, -real2, pgev(rbvevd(m_data, alpha = fit_ct_low$estimate[1], beta = fit_ct_low$estimate[2], model = "ct")[,2]), thresh_low[2], 0.90)),
  MEGPD = cbind(ssX2, ssY2)
)

################################################################################
# Panel plotting function
################################################################################
draw_panel <- function(th, ob, list_data, title, colors) {
  plot(th, ob, col = colors[1], pch = 20,
       xlab = "", ylab = "", main = title,
       xlim = range(th, na.rm=TRUE),
       ylim = range(ob, na.rm=TRUE),
       asp  = 1, panel.first = grid(nx = NULL, ny = NULL, col="lightgrey", lty="dotted"))
  abline(0,1, lty=2, col="gray50")
  for (i in 2:length(list_data)) {
    points(list_data[[i]]$Theoretical,
           list_data[[i]]$Observed,
           col = colors[i], pch = 20)
  }
}

################################################################################
# Final Plot: 2×3 layout (upper = first row, lower = second row)
################################################################################
colors <- c("purple2", "goldenrod", "steelblue", "seagreen3", "tomato")
set_labels <- names(sim_data_list_up)

par(mfrow = c(2,3), oma = c(5,5,4,8), mar = c(2,2,1,1), xpd = FALSE)

# Upper
draw_panel(quantile(sim_data_list_up[[1]][,1], seq(0,1,0.01)),
           quantile(observed_data[,1], seq(0,1,0.01)),
           lapply(sim_data_list_up, function(sim) data.frame(Theoretical = quantile(sim[,1], seq(0,1,0.01)), Observed = quantile(observed_data[,1], seq(0,1,0.01)))),
           "Station 11", colors)

draw_panel(quantile(sim_data_list_up[[1]][,2], seq(0,1,0.01)),
           quantile(observed_data[,2], seq(0,1,0.01)),
           lapply(sim_data_list_up, function(sim) data.frame(Theoretical = quantile(sim[,2], seq(0,1,0.01)), Observed = quantile(observed_data[,2], seq(0,1,0.01)))),
           "Station 23", colors)

draw_panel(quantile(rowSums(sim_data_list_up[[1]]), seq(0,1,0.01)),
           quantile(rowSums(observed_data), seq(0,1,0.01)),
           lapply(sim_data_list_up, function(sim) data.frame(Theoretical = quantile(rowSums(sim), seq(0,1,0.01)), Observed = quantile(rowSums(observed_data), seq(0,1,0.01)))),
           "Joint (11+23)", colors)

# Lower
draw_panel(quantile(sim_data_list_low[[1]][,1], seq(0,1,0.01)),
           quantile(real1, seq(0,1,0.01)),
           lapply(sim_data_list_low, function(sim) data.frame(Theoretical = quantile(sim[,1], seq(0,1,0.01)), Observed = quantile(real1, seq(0,1,0.01)))),
           "", colors)

draw_panel(quantile(sim_data_list_low[[1]][,2], seq(0,1,0.01)),
           quantile(real2, seq(0,1,0.01)),
           lapply(sim_data_list_low, function(sim) data.frame(Theoretical = quantile(sim[,2], seq(0,1,0.01)), Observed = quantile(real2, seq(0,1,0.01)))),
           "", colors)

draw_panel(quantile(rowSums(sim_data_list_low[[1]]), seq(0,1,0.01)),
           quantile(rowSums(cbind(real1,real2)), seq(0,1,0.01)),
           lapply(sim_data_list_low, function(sim) data.frame(Theoretical = quantile(rowSums(sim), seq(0,1,0.01)), Observed = quantile(rowSums(cbind(real1,real2)), seq(0,1,0.01)))),
           "", colors)

# Axis labels
mtext("Theoretical Quantiles", side=1, outer=TRUE, line=2.5, cex=.8)
mtext("Sample Quantiles", side=2, outer=TRUE, line=2.5, cex=.8)

# Legend
par(xpd = NA)
legend("right", legend = set_labels, col = colors, pch = 20, bty = "n", cex = 0.8)

# Save plot
dir.create("Figures/Application", showWarnings = FALSE, recursive = TRUE)
png("Figures/Application/qq_panels_11_23.png", width=8, height=6, units="in", res=300)
dev.off()
