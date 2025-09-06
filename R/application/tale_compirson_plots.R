library(mev)        # for fit.gpd() and fbvpot()
library(evd)        # for pgev() and rbvevd()
library(ggplot2)
library(patchwork)  # for shared-legend layouts
library(RColorBrewer)

################################################################################
#  Upper‐tail QQ comparison for stations 11 & 23
################################################################################


# 1. Define Fgpd_inv (same as in old code)

Fgpd_inv <- function(par, data, simdata, threshold_data, threshold_simdata) {
  # par: c(sigma, xi) for GPD on data margin
  sigma <- par[1]
  xi    <- par[2]
  L     <- threshold_simdata 
  L_data<- threshold_data
  n     <- length(simdata)
  # proportion of simulated exceeding threshold L
  delta_n <- sum(simdata >= L)/n
  trans_gpd <- numeric(n)
  
  for (i in seq_len(n)) {
    if (simdata[i] >= L) {
      p_value <- 1 - (1 - simdata[i])/delta_n  # for qgpd
      # Adjust for edge cases
      if (p_value <= 0) {
        p_value <- abs(p_value)
      } else if (p_value >= 1) {
        p_value <- 0.9999
      }
      # threshold_data + GPD excess
      trans_gpd[i] <- L_data + qgpd(p_value, loc = 0, scale = sigma, shape = xi)
    } else {
      # below threshold: transform via empirical quantile of the data margin
      # simdata[i] ∈ [0,1), treat as probability
      trans_gpd[i] <- as.numeric(quantile(data, probs = simdata[i], na.rm = TRUE))
    }
  }
  return(trans_gpd)
}
# 2. Load raw rainfall data and extract station 11 & 23
################################################################################

# call real rainfall data:
# real1 <- mat_11_23_scaled[1,]
# real2 <- mat_11_23_scaled[2,]

# real1 <- mat_11_15_scaled[1,]
# real2 <- mat_11_15_scaled[2,]

real1 <- mat_15_23_scaled[1,]
real2 <- mat_15_23_scaled[2,]

# If you prefer threshold-fitting including zeros: adjust accordingly.
# But here we focus on positive margins.

# Empirical chi-plot for upper tail:
par(mfrow = c(1,2))
chiplot(cbind(real2, real1))

################################################################################
# 3. Define thresholds for tail-fitting (upper joint tail)
################################################################################
# Use e.g. 95th percentile of each margin (on positive-only series real1, real2):
thresh <- c(
  as.numeric(quantile(real1, 0.95, na.rm = TRUE)),
  as.numeric(quantile(real2, 0.95, na.rm = TRUE))
)
cat("Thresholds for station 15 & 23 (95%):", thresh, "\n")

# On the scatter, plot threshold lines:
par(mfrow = c(1,1))
plot(real1, real2, pch = 20, cex = 0.6,
     xlab = "Station 11", ylab = "Station 23",
     main = "Threshold lines for tail fitting")
abline(v = thresh[1], h = thresh[2], col = "red", lty = 2)

################################################################################
# 4. Fit marginal GPD (univariate) for each station upper tail
################################################################################
# fit.gpd from mev: syntax fit.gpd(data, threshold)
fit_up1 <- fit.gpd(real1, thresh[1])
fit_up2 <- fit.gpd(real2, thresh[2])
cat("Station 15 GPD fit params (scale, shape):", fit_up1$param, "\n")
cat("Station 23 GPD fit params (scale, shape):", fit_up2$param, "\n")

################################################################################
# 5. Fit bivariate tail models (classical) via fbvpot
################################################################################
# fbvpot: x = data frame or matrix with two columns, threshold vector, plus initial scale/shape
bivariate <- data.frame(real1, real2)
# Fit logistic, bilogistic, HR, CT models:
fit1.test.up <- fbvpot(x = bivariate, threshold = thresh,
                       model = "log",
                       scale1 = fit_up1$param[1], scale2 = fit_up2$param[1],
                       shape1 = fit_up1$param[2], shape2 = fit_up2$param[2],
                       std.err = FALSE)
fit2.test.up <- fbvpot(x = bivariate, threshold = thresh,
                       model = "bilog",
                       scale1 = fit_up1$param[1], scale2 = fit_up2$param[1],
                       shape1 = fit_up1$param[2], shape2 = fit_up2$param[2],
                       std.err = FALSE)
fit3.test.up <- fbvpot(x = bivariate, threshold = thresh,
                       model = "hr",
                       scale1 = fit_up1$param[1], scale2 = fit_up2$param[1],
                       shape1 = fit_up1$param[2], shape2 = fit_up2$param[2],
                       std.err = FALSE)
fit4.test.up <- fbvpot(x = bivariate, threshold = thresh,
                       model = "ct",
                       scale1 = fit_up1$param[1], scale2 = fit_up2$param[1],
                       shape1 = fit_up1$param[2], shape2 = fit_up2$param[2],
                       std.err = FALSE)
# Print summaries:
print(fit1.test.up)
print(fit2.test.up)
print(fit3.test.up)
print(fit4.test.up)

################################################################################
# 6. Simulate from each fitted classical model, with GEV margins
################################################################################
# m = number of replicates for simulation: choose m = length(real1) to match data size,
# or m = large (1e7) if desired. Here we set m_data = length(real1) for direct comparison:
m_data <- nrow(bivariate)

# Simulate from each model:
sim.fit1.up <- rbvevd(m_data, dep = fit1.test.up$estimate, model = "log")
sim.fit2.up <- rbvevd(m_data, alpha = fit2.test.up$estimate[1],
                      beta = fit2.test.up$estimate[2], model = "bilog")
sim.fit3.up <- rbvevd(m_data, dep = fit3.test.up$estimate, model = "hr")
sim.fit4.up <- rbvevd(m_data, alpha = fit4.test.up$estimate[1],
                      beta = fit4.test.up$estimate[2], model = "ct")
# Each sim.fit*_up is an m_data × 2 numeric matrix, with GEV margins (unit Fréchet or GEV?).
# We will transform these to uniform via pgev, then to GPD margins matching real data via Fgpd_inv.

################################################################################
# 7. Transform simulated GEV to uniform [0,1], then to GPD margins on original data
################################################################################
# First, transform each simulated margin to uniform via GEV CDF:
# Note: pgev in evd expects GEV parameters? By default pgev(x) uses default location=0, scale=1, shape=1? 
# But in river code they used pgev(sim.fit1_up[,1]) directly. This assumes sim.fit1.up margins follow the GEV fitted? 
# In classical EV, rbvevd produces unit GEV margins with some parameterization. 
# We follow the same pattern:
cdf1.up_loc1_ <- pgev(sim.fit1.up[,1])
cdf1.up_loc2_ <- pgev(sim.fit1.up[,2])
cdf2.up_loc1_ <- pgev(sim.fit2.up[,1])
cdf2.up_loc2_ <- pgev(sim.fit2.up[,2])
cdf3.up_loc1_ <- pgev(sim.fit3.up[,1])
cdf3.up_loc2_ <- pgev(sim.fit3.up[,2])
cdf4.up_loc1_ <- pgev(sim.fit4.up[,1])
cdf4.up_loc2_ <- pgev(sim.fit4.up[,2])

# Now transform to GPD margins matching real data via Fgpd_inv:
thresh1 <- c(0.95, 0.95)  # threshold on uniform scale: same 0.95 used in river code
# But ensure these align: threshold_simdata = 0.95 means upper 5%. 
# threshold_data should be thresh from real data in original units.

gpd1.up_loc1 <- Fgpd_inv(par = c(fit_up1$param[1], fit_up1$param[2]),
                         data = real1,
                         simdata = cdf1.up_loc1_,
                         threshold_data = thresh[1],
                         threshold_simdata = thresh1[1])
gpd1.up_loc2 <- Fgpd_inv(par = c(fit_up2$param[1], fit_up2$param[2]),
                         data = real2,
                         simdata = cdf1.up_loc2_,
                         threshold_data = thresh[2],
                         threshold_simdata = thresh1[2])

gpd2.up_loc1 <- Fgpd_inv(par = c(fit_up1$param[1], fit_up1$param[2]),
                         data = real1,
                         simdata = cdf2.up_loc1_,
                         threshold_data = thresh[1],
                         threshold_simdata = thresh1[1])
gpd2.up_loc2 <- Fgpd_inv(par = c(fit_up2$param[1], fit_up2$param[2]),
                         data = real2,
                         simdata = cdf2.up_loc2_,
                         threshold_data = thresh[2],
                         threshold_simdata = thresh1[2])

gpd3.up_loc1 <- Fgpd_inv(par = c(fit_up1$param[1], fit_up1$param[2]),
                         data = real1,
                         simdata = cdf3.up_loc1_,
                         threshold_data = thresh[1],
                         threshold_simdata = thresh1[1])
gpd3.up_loc2 <- Fgpd_inv(par = c(fit_up2$param[1], fit_up2$param[2]),
                         data = real2,
                         simdata = cdf3.up_loc2_,
                         threshold_data = thresh[2],
                         threshold_simdata = thresh1[2])

gpd4.up_loc1 <- Fgpd_inv(par = c(fit_up1$param[1], fit_up1$param[2]),
                         data = real1,
                         simdata = cdf4.up_loc1_,
                         threshold_data = thresh[1],
                         threshold_simdata = thresh1[1])
gpd4.up_loc2 <- Fgpd_inv(par = c(fit_up2$param[1], fit_up2$param[2]),
                         data = real2,
                         simdata = cdf4.up_loc2_,
                         threshold_data = thresh[2],
                         threshold_simdata = thresh1[2])

################################################################################
# 8. Prepare observed exceedances matrix and MEGPD simulated data
################################################################################
# Observed exceedances: here we take the raw positive pair (real1, real2) above threshold?
# In your old code: riverdata13 <- cbind(WyeSummerFlow_data[,1], WyeSummerFlow_data[,3])
# They then applied Fgpd_inv only for simulated; for observed_data they passed full bivariate.
observed_data <- cbind(real1, real2)

#MEGPDsample <- cbind(ssX1, ssY1)
#MEGPDsample <- cbind(ssX2, ssY2)
MEGPDsample <- cbind(ssX3, ssY3)


################################################################################
# 9. Assemble list of simulated data sets for comparison
################################################################################
sim_data_list_up <- list(
  gpd_log   = cbind(gpd1.up_loc1, gpd1.up_loc2),
  gpd_bilog = cbind(gpd2.up_loc1, gpd2.up_loc2),
  gpd_hr    = cbind(gpd3.up_loc1, gpd3.up_loc2),
  gpd_ct    = cbind(gpd4.up_loc1, gpd4.up_loc2),
  megpd     = MEGPDsample
)
set_labels <- names(sim_data_list_up)  # c("gpd_log","gpd_bilog","gpd_hr","gpd_ct","megpd")

################################################################################
# 10. Diagnostic QQ plots (margin 1, margin 2, joint) comparing observed vs each sim
################################################################################
# Use ggplot2 for nicer output. We compute quantiles at seq(0,1,by=0.01).
probs <- seq(0, 1, by = 0.01)

# 10a. Prepare data frames for margin 1:
sim_data_df_list1_up <- lapply(seq_along(sim_data_list_up), function(i) {
  simdat <- sim_data_list_up[[i]][,1]
  data.frame(
    Theoretical = quantile(simdat, probs = probs, na.rm = TRUE),
    Observed    = quantile(observed_data[,1], probs = probs, na.rm = TRUE),
    Model       = set_labels[i]
  )
})
# Combine into one data.frame
df1 <- do.call(rbind, lapply(seq_along(sim_data_df_list1_up), function(i) {
  cbind(sim_data_df_list1_up[[i]], Model = set_labels[i])
}))

# 10b. Margin 2:
sim_data_df_list2_up <- lapply(seq_along(sim_data_list_up), function(i) {
  simdat <- sim_data_list_up[[i]][,2]
  data.frame(
    Theoretical = quantile(simdat, probs = probs, na.rm = TRUE),
    Observed    = quantile(observed_data[,2], probs = probs, na.rm = TRUE),
    Model       = set_labels[i]
  )
})
df2 <- do.call(rbind, lapply(seq_along(sim_data_df_list2_up), function(i) {
  cbind(sim_data_df_list2_up[[i]], Model = set_labels[i])
}))

# 10c. Joint sum margin (X+Y):
sim_data_df_list3_up <- lapply(seq_along(sim_data_list_up), function(i) {
  simdat <- rowSums(sim_data_list_up[[i]])
  data.frame(
    Theoretical = quantile(simdat, probs = probs, na.rm = TRUE),
    Observed    = quantile(observed_data[,1] + observed_data[,2], probs = probs, na.rm = TRUE),
    Model       = set_labels[i]
  )
})
df3 <- do.call(rbind, lapply(seq_along(sim_data_df_list3_up), function(i) {
  cbind(sim_data_df_list3_up[[i]], Model = set_labels[i])
}))

# 10d. Plotting with ggplot2: three side-by-side QQ plots
# Choose distinct colors:
palette_colors <- RColorBrewer::brewer.pal(max(3, length(set_labels)), "Set1")
names(palette_colors) <- set_labels

# Function to produce one QQ ggplot
make_qq_plot <- function(df, xlab_title, ylab_title) {
  ggplot(df, aes(x = Theoretical, y = Observed, color = Model)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    scale_color_manual(values = palette_colors) +
    labs(x = xlab_title, y = ylab_title) +
    coord_equal() +
    theme_minimal(base_size = 12) +
    theme(
      panel.border = element_rect(color = "black", fill = NA),
      panel.background = element_blank(),
      strip.background = element_blank(),
      legend.position = "bottom"
    )
}

p_margin1 <- make_qq_plot(df1, xlab_title = "Theoretical Quantiles", ylab_title = "Observed Quantiles (Station 11)")
p_margin2 <- make_qq_plot(df2, xlab_title = "Theoretical Quantiles", ylab_title = "Observed Quantiles (Station 23)")
p_joint   <- make_qq_plot(df3, xlab_title = "Theoretical Quantiles", ylab_title = "Observed Quantiles (Joint sum)")

# Arrange side by side with a shared legend
library(patchwork)
combined_qq <- (p_margin1 + p_margin2 + p_joint) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
print(combined_qq)

################################################################################
# End of script
################################################################################

######################################################
library(ggplot2)
# Create a list of labels for simulation sets (customize as needed)
set_labels <- c("Log", "Bilog", "HR", "CT", "MEGPD")


# Create a list of data frames for each simulation set and assign labels

## For margin 1
sim_data_df_list1_up <-  lapply(1:length(sim_data_list_up), function(i) {
  data.frame(Theoretical = quantile(sim_data_list_up[[i]][,1], probs = seq(0, 1, by = 0.01)),
             Observed = quantile(observed_data[,1], probs = seq(0, 1, by = 0.01)),
             Model = set_labels[i])})

## For margin 2
sim_data_df_list2_up <-  lapply(1:length(sim_data_list_up), function(i) {
  data.frame(Theoretical = quantile(sim_data_list_up[[i]][,2], probs = seq(0, 1, by = 0.01)),
             Observed = quantile(observed_data[,2], probs = seq(0, 1, by = 0.01)),
             Model = set_labels[i])})

## For the joint (margin 1+2)
sim_data_df_list3_up <-  lapply(1:length(sim_data_list_up), function(i) {
  data.frame(Theoretical = quantile(sim_data_list_up[[i]][,1]+sim_data_list_up[[i]][,2], probs = seq(0, 1, by = 0.01)),
             Observed = quantile(observed_data[,1]+observed_data[,2], probs = seq(0, 1, by = 0.01)),
             Model = set_labels[i])})

# Combine the data frames into one
sim_data_df1_up <- do.call(rbind, sim_data_df_list1_up)
sim_data_df2_up <- do.call(rbind, sim_data_df_list2_up)
sim_data_df3_up <- do.call(rbind, sim_data_df_list3_up)
######################################################


# Set the desired font size for axis labels and tick marks
axis_label_size <- 1.2  # Adjust this value as needed
axis_tick_size <- 1.0   # Adjust this value as needed

# QQ plot for margin 1
par(mfrow = c(1, 3))
# Define a vector of colors matching the "Accent" palette
colors <- c("#FFFF99", "#386CB0", "#7FC97F", "#BEAED4", "#FDC086")

plot(sim_data_df_list1_up[[1]]$Theoretical, sim_data_df_list1_up[[1]]$Observed,
     col = colors[1], pch = 20, xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles", main = "",
     xlim = c(0, max(sim_data_df1_up$Theoretical)), 
     ylim = c(0, max(sim_data_df1_up$Observed)),
     cex.axis = axis_tick_size, cex.lab = axis_label_size)  # Adjust font sizes

abline(0, 1, col = "black", lty = 2)
for (i in 2:length(sim_data_df_list1_up)) {
  points(sim_data_df_list1_up[[i]]$Theoretical, sim_data_df_list1_up[[i]]$Observed, pch = 20, col = colors[i])
  legend("bottomright", legend = set_labels, col = colors, pch = 20, cex = 1)
}
#abline(h=thresh[2],v=thresh[1], col = "red", lty = 2)


# QQ plot for margin 2
plot(sim_data_df_list2_up[[1]]$Theoretical, sim_data_df_list2_up[[1]]$Observed,
     col = colors[1], pch = 20, xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles", main = "",
     xlim = c(0, max(sim_data_df2_up$Theoretical)), 
     ylim = c(0, max(sim_data_df2_up$Observed)),
     cex.axis = axis_tick_size, cex.lab = axis_label_size)  # Adjust font sizes

abline(0, 1, col = "black", lty = 2)
for (i in 2:length(sim_data_df_list2_up)) {
  points(sim_data_df_list2_up[[i]]$Theoretical, sim_data_df_list2_up[[i]]$Observed, pch = 20, col = colors[i])
  legend("bottomright", legend = set_labels, col = colors, pch = 20, cex = 1)
}
#abline(h=thresh[2],v=thresh[1], col = "red", lty = 2)


# QQ plot for joint margins
plot(sim_data_df_list3_up[[1]]$Theoretical, sim_data_df_list3_up[[1]]$Observed,
     col = colors[1], pch = 20, xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles", main = "",
     xlim = c(0, max(sim_data_df3_up$Theoretical)), 
     ylim = c(0, max(sim_data_df3_up$Observed)),
     cex.axis = axis_tick_size, cex.lab = axis_label_size)  # Adjust font sizes

abline(0, 1, col = "black", lty = 2)
for (i in 2:length(sim_data_df_list3_up)) {
  points(sim_data_df_list3_up[[i]]$Theoretical, sim_data_df_list3_up[[i]]$Observed, pch = 20, col = colors[i])
  legend("bottomright", legend = set_labels, col = colors, pch = 20, cex = 1)
}
#abline(h=thresh[2],v=thresh[1], col = "red", lty = 2)

# Create a list of labels for simulation sets (customize as needed)
set_labels <- c("Log", "Bilog", "HR", "CT", "MEGPD")


# Create a list of data frames for each simulation set and assign labels

## For margin 1
sim_data_df_list1_up <-  lapply(1:length(sim_data_list_up), function(i) {
  data.frame(Theoretical = quantile(sim_data_list_up[[i]][,1], probs = seq(0, 1, by = 0.01)),
             Observed = quantile(observed_data[,1], probs = seq(0, 1, by = 0.01)),
             Model = set_labels[i])})

## For margin 2
sim_data_df_list2_up <-  lapply(1:length(sim_data_list_up), function(i) {
  data.frame(Theoretical = quantile(sim_data_list_up[[i]][,2], probs = seq(0, 1, by = 0.01)),
             Observed = quantile(observed_data[,2], probs = seq(0, 1, by = 0.01)),
             Model = set_labels[i])})

## For the joint (margin 1+2)
sim_data_df_list3_up <-  lapply(1:length(sim_data_list_up), function(i) {
  data.frame(Theoretical = quantile(sim_data_list_up[[i]][,1]+sim_data_list_up[[i]][,2], probs = seq(0, 1, by = 0.01)),
             Observed = quantile(observed_data[,1]+observed_data[,2], probs = seq(0, 1, by = 0.01)),
             Model = set_labels[i])})

# Combine the data frames into one
sim_data_df1_up <- do.call(rbind, sim_data_df_list1_up)
sim_data_df2_up <- do.call(rbind, sim_data_df_list2_up)
sim_data_df3_up <- do.call(rbind, sim_data_df_list3_up)
######################################################


# Set the desired font size for axis labels and tick marks
axis_label_size <- 1.2  # Adjust this value as needed
axis_tick_size <- 1.0   # Adjust this value as needed

# QQ plot for margin 1
par(mfrow = c(1, 3))
# Define a vector of colors matching the "Accent" palette
colors <- c("#FFFF99", "#386CB0", "#7FC97F", "#BEAED4", "#FDC086")

upper_tale1 <- plot(sim_data_df_list1_up[[1]]$Theoretical, sim_data_df_list1_up[[1]]$Observed,
     col = colors[1], pch = 20, xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles", main = "",
     xlim = c(0, max(sim_data_df1_up$Theoretical)), 
     ylim = c(0, max(sim_data_df1_up$Observed)),
     cex.axis = axis_tick_size, cex.lab = axis_label_size)  # Adjust font sizes

abline(0, 1, col = "black", lty = 2)
for (i in 2:length(sim_data_df_list1_up)) {
  points(sim_data_df_list1_up[[i]]$Theoretical, sim_data_df_list1_up[[i]]$Observed, pch = 20, col = colors[i])
  legend("bottomright", legend = set_labels, col = colors, pch = 20, cex = 1)
}
#abline(h=thresh[2],v=thresh[1], col = "red", lty = 2)


# QQ plot for margin 2
plot(sim_data_df_list2_up[[1]]$Theoretical, sim_data_df_list2_up[[1]]$Observed,
     col = colors[1], pch = 20, xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles", main = "",
     xlim = c(0, max(sim_data_df2_up$Theoretical)), 
     ylim = c(0, max(sim_data_df2_up$Observed)),
     cex.axis = axis_tick_size, cex.lab = axis_label_size)  # Adjust font sizes

abline(0, 1, col = "black", lty = 2)
for (i in 2:length(sim_data_df_list2_up)) {
  points(sim_data_df_list2_up[[i]]$Theoretical, sim_data_df_list2_up[[i]]$Observed, pch = 20, col = colors[i])
  legend("bottomright", legend = set_labels, col = colors, pch = 20, cex = 1)
}
#abline(h=thresh[2],v=thresh[1], col = "red", lty = 2)


# QQ plot for joint margins
upper_joint <- plot(sim_data_df_list3_up[[1]]$Theoretical, sim_data_df_list3_up[[1]]$Observed,
     col = colors[1], pch = 20, xlab = "Theoretical Quantiles",
     ylab = "Sample Quantiles", main = "",
     xlim = c(0, max(sim_data_df3_up$Theoretical)), 
     ylim = c(0, max(sim_data_df3_up$Observed)),
     cex.axis = axis_tick_size, cex.lab = axis_label_size)  # Adjust font sizes

abline(0, 1, col = "black", lty = 2)
for (i in 2:length(sim_data_df_list3_up)) {
  points(sim_data_df_list3_up[[i]]$Theoretical, sim_data_df_list3_up[[i]]$Observed, pch = 20, col = colors[i])
  legend("bottomright", legend = set_labels, col = colors, pch = 20, cex = 1)
}
#abline(h=thresh[2],v=thresh[1], col = "red", lty = 2)

################################################################################
##  Lower‐tail QQ comparison for stations 11 & 23
################################################################################


## 1. Real rainfall → lower tail (negate so we look at “exceedances” below)
# real11 <- -1 * mat_11_23_scaled[1, ]
# real23 <- -1 * mat_11_23_scaled[2, ]

# real11 <- -1 * mat_11_15_scaled[1, ]
# real23 <- -1 * mat_11_15_scaled[2, ]

real11 <- -1 * mat_15_23_scaled[1, ]
real23 <- -1 * mat_15_23_scaled[2, ]

## 2. 90% thresholds on the negated series
threshL <- c(
  as.numeric(quantile(real11, 0.90, na.rm=TRUE)),
  as.numeric(quantile(real23, 0.90, na.rm=TRUE))
)
cat("Lower-tail 90% thresholds:", threshL, "\n")

## 3. Univariate GPD fits at those thresholds
fit_low11 <- fit.gpd(real11, threshL[1])
fit_low23 <- fit.gpd(real23, threshL[2])

## 4. Bivariate censored GPD fits for the lower tail
biv_low <- data.frame(real11, real23)

fit1.test.low <- fbvpot(
  x         = biv_low,
  threshold = threshL,
  model     = "log",
  scale1    = fit_low11$param[1],
  shape1    = fit_low11$param[2]+0.2,
  scale2    = fit_low23$param[1],
  shape2    = fit_low23$param[2]+0.2,
  std.err   = FALSE
)
fit2.test.low <- fbvpot(
  x         = biv_low, threshold = threshL, model = "bilog",
  scale1    = fit_low11$param[1], shape1 = fit_low11$param[2]+0.2,
  scale2    = fit_low23$param[1], shape2 = fit_low23$param[2]+0.2,
  std.err   = FALSE
)
fit3.test.low <- fbvpot(
  x         = biv_low, threshold = threshL, model = "hr",
  scale1    = fit_low11$param[1], shape1 = fit_low11$param[2]+0.2,
  scale2    = fit_low23$param[1], shape2 = fit_low23$param[2]+0.2,
  std.err   = FALSE
)
fit4.test.low <- fbvpot(
  x         = biv_low, threshold = threshL, model = "ct",
  scale1    = fit_low11$param[1], shape1 = fit_low11$param[2]+0.2,
  scale2    = fit_low23$param[1], shape2 = fit_low23$param[2]+0.2,
  std.err   = FALSE
)

## 5. Simulate from each model (match real data length)
m_data     <- length(real11)
sim1.low   <- rbvevd(m_data, dep   = fit1.test.low$estimate, model = "log")
sim2.low   <- rbvevd(m_data, alpha = fit2.test.low$estimate[1],
                     beta  = fit2.test.low$estimate[2], model = "bilog")
sim3.low   <- rbvevd(m_data, dep   = fit3.test.low$estimate, model = "hr")
sim4.low   <- rbvevd(m_data, alpha = fit4.test.low$estimate[1],
                     beta  = fit4.test.low$estimate[2], model = "ct")

## 6. GEV → Uniform
u1_1 <- pgev(sim1.low[,1]); u1_2 <- pgev(sim1.low[,2])
u2_1 <- pgev(sim2.low[,1]); u2_2 <- pgev(sim2.low[,2])
u3_1 <- pgev(sim3.low[,1]); u3_2 <- pgev(sim3.low[,2])
u4_1 <- pgev(sim4.low[,1]); u4_2 <- pgev(sim4.low[,2])

## 7. Uniform → your data’s GPD margins via Fgpd_inv()
#    threshold_simdata = 0.90  since we used the 90% quantile above
tsim <- 0.90
gpd1.low_11 <- Fgpd_inv(
  par             = fit_low11$param, data = real11,
  simdata         = u1_1, threshold_data = threshL[1],
  threshold_simdata = tsim
)
gpd1.low_23 <- Fgpd_inv(
  par             = fit_low23$param, data = real23,
  simdata         = u1_2, threshold_data = threshL[2],
  threshold_simdata = tsim
)
gpd2.low_11 <- Fgpd_inv(fit_low11$param, real11, u2_1, threshL[1], tsim)
gpd2.low_23 <- Fgpd_inv(fit_low23$param, real23, u2_2, threshL[2], tsim)
gpd3.low_11 <- Fgpd_inv(fit_low11$param, real11, u3_1, threshL[1], tsim)
gpd3.low_23 <- Fgpd_inv(fit_low23$param, real23, u3_2, threshL[2], tsim)
gpd4.low_11 <- Fgpd_inv(fit_low11$param, real11, u4_1, threshL[1], tsim)
gpd4.low_23 <- Fgpd_inv(fit_low23$param, real23, u4_2, threshL[2], tsim)

## 8. Assemble MEGPD samples (negate to go to “lower” tail)
MEGPDsample_low <- cbind(
  `MEGPD_11` = ssX3,
  `MEGPD_23` = ssY3
)

## 9. Put all 5 simulated matrices in a named list
sim_data_list_low <- list(
  gpd_log   = -cbind(gpd1.low_11, gpd1.low_23),
  gpd_bilog = -cbind(gpd2.low_11, gpd2.low_23),
  gpd_hr    = -cbind(gpd3.low_11, gpd3.low_23),
  gpd_ct    = -cbind(gpd4.low_11, gpd4.low_23),
  megpd     = MEGPDsample_low
)
set_labels_low <- names(sim_data_list_low)

## 10. Build three QQ data-frames (margins 1 & 2 and joint sum)
probs <- seq(0,1,by=0.01)
observed_low <- cbind(`Obs_11` = real1, `Obs_23` = real2)

qq_list1_low <- lapply(seq_along(sim_data_list_low), function(i) {
  simv <- sim_data_list_low[[i]][,1]
  data.frame(
    Theoretical = quantile(simv,   probs, na.rm=TRUE),
    Observed    = quantile(observed_low[,1], probs, na.rm=TRUE),
    Model       = set_labels_low[i]
  )
})
df1_low <- do.call(rbind, qq_list1_low)

qq_list2_low <- lapply(seq_along(sim_data_list_low), function(i) {
  simv <- sim_data_list_low[[i]][,2]
  data.frame(
    Theoretical = quantile(simv,   probs, na.rm=TRUE),
    Observed    = quantile(observed_low[,2], probs, na.rm=TRUE),
    Model       = set_labels_low[i]
  )
})
df2_low <- do.call(rbind, qq_list2_low)

qq_list3_low <- lapply(seq_along(sim_data_list_low), function(i) {
  simv <- rowSums(sim_data_list_low[[i]])
  data.frame(
    Theoretical = quantile(simv,   probs, na.rm=TRUE),
    Observed    = quantile(rowSums(observed_low), probs, na.rm=TRUE),
    Model       = set_labels_low[i]
  )
})
df3_low <- do.call(rbind, qq_list3_low)

## 11. Plot with ggplot2 / patchwork, matching your earlier style
palette_low <- RColorBrewer::brewer.pal(length(set_labels_low), "Set1")
names(palette_low) <- set_labels_low

make_qq <- function(df, xlab, ylab) {
  ggplot(df, aes(x=Theoretical, y=Observed, color=Model)) +
    geom_ribbon(aes(ymin=Theoretical, ymax=Theoretical), fill=NA) +  # no ribbon here
    geom_point(size=0.8, alpha=0.6) +
    geom_abline(slope=1, intercept=0, linetype="dashed", color="gray50") +
    scale_color_manual(values=palette_low) +
    labs(x=xlab, y=ylab) +
    coord_equal() +
    theme_minimal(base_size=11) +
    theme(
      panel.border      = element_rect(color="black", fill=NA),
      panel.background  = element_blank(),
      legend.position   = "bottom",
      strip.background  = element_blank()
    )
}

p1_low <- make_qq(df1_low, xlab="Theoretical Quantiles", ylab="Observed Quantiles (11)")
p2_low <- make_qq(df2_low, xlab="Theoretical Quantiles", ylab="Observed Quantiles (23)")
p3_low <- make_qq(df3_low, xlab="Theoretical Quantiles", ylab="Observed Quantiles (Joint sum)")

(p1_low + p2_low + p3_low) +
  plot_layout(guides="collect") & theme(legend.position="bottom")
#############################################
#better aesthitic plot:

set_labels <- c("Log","Bilog","HR","CT","MEGPD")
colors    <- c("#FFFF99","#386CB0","#7FC97F","#BEAED4","#FDC086")

# build the three lists of data.frames exactly as before
## margin 1 (station 11)
sim_data_df_list1_low <- lapply(seq_along(sim_data_list_low), function(i) {
  data.frame(
    Theoretical = quantile(sim_data_list_low[[i]][,1],
                           probs = seq(0,1,by=0.01), na.rm=TRUE),
    Observed    = quantile(observed_low[,1],
                           probs = seq(0,1,by=0.01), na.rm=TRUE),
    Model       = set_labels[i]
  )
})

## margin 2 (station 23)
sim_data_df_list2_low <- lapply(seq_along(sim_data_list_low), function(i) {
  data.frame(
    Theoretical = quantile(sim_data_list_low[[i]][,2],
                           probs = seq(0,1,by=0.01), na.rm=TRUE),
    Observed    = quantile(observed_low[,2],
                           probs = seq(0,1,by=0.01), na.rm=TRUE),
    Model       = set_labels[i]
  )
})

## joint sum (11+23)
sim_data_df_list3_low <- lapply(seq_along(sim_data_list_low), function(i) {
  data.frame(
    Theoretical = quantile(rowSums(sim_data_list_low[[i]]),
                           probs = seq(0,1,by=0.01), na.rm=TRUE),
    Observed    = quantile(rowSums(observed_low),
                           probs = seq(0,1,by=0.01), na.rm=TRUE),
    Model       = set_labels[i]
  )
})

# 3. glue each list into one big data.frame
sim_data_df1_low <- do.call(rbind, sim_data_df_list1_low)
sim_data_df2_low <- do.call(rbind, sim_data_df_list2_low)
sim_data_df3_low <- do.call(rbind, sim_data_df_list3_low)

# 4. axis/font sizes
axis_label_size <- 1.2
axis_tick_size  <- 1.0

# 5. now the three side‐by‐side base‐R QQ panels

par(mfrow=c(1,3))

## margin 1
low_tale1 <-plot(
  sim_data_df_list1_low[[1]]$Theoretical,
  sim_data_df_list1_low[[1]]$Observed,
  col = colors[1], pch = 20,
  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
  main = "",
  xlim = c(0, max(sim_data_df1_low$Theoretical)),
  ylim = c(0, max(sim_data_df1_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, col="black", lty=2)
for(i in 2:length(sim_data_df_list1_low)){
  points(
    sim_data_df_list1_low[[i]]$Theoretical,
    sim_data_df_list1_low[[i]]$Observed,
    pch=20, col=colors[i]
  )
}
legend(
  "bottomright", legend=set_labels,
  col=colors, pch=20, cex=0.8
)

## margin 2
low_tale2 <- plot(
  sim_data_df_list2_low[[1]]$Theoretical,
  sim_data_df_list2_low[[1]]$Observed,
  col = colors[1], pch = 20,
  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
  main = "",
  xlim = c(0, max(sim_data_df2_low$Theoretical)),
  ylim = c(0, max(sim_data_df2_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, col="black", lty=2)
for(i in 2:length(sim_data_df_list2_low)){
  points(
    sim_data_df_list2_low[[i]]$Theoretical,
    sim_data_df_list2_low[[i]]$Observed,
    pch=20, col=colors[i]
  )
}
legend(
  "bottomright", legend=set_labels,
  col=colors, pch=20, cex=0.8
)

## joint sum
joint_low <- plot(
  sim_data_df_list3_low[[1]]$Theoretical,
  sim_data_df_list3_low[[1]]$Observed,
  col = colors[1], pch = 20,
  xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
  main = "",
  xlim = c(0, max(sim_data_df3_low$Theoretical)),
  ylim = c(0, max(sim_data_df3_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, col="black", lty=2)
for(i in 2:length(sim_data_df_list3_low)){
  points(
    sim_data_df_list3_low[[i]]$Theoretical,
    sim_data_df_list3_low[[i]]$Observed,
    pch=20, col=colors[i]
  )
}
legend(
  "bottomright", legend=set_labels,
  col=colors, pch=20, cex=0.8
)

####################################
### combine all 6 figures:

library(patchwork)



# And that both sets use the same `set_labels` and `palette_colors`
#  upper_tale1 + upper_tale2 + upper_joint + lower_tale1   + lower_tale2 + lower_joint 

combined_qq_all <- (
  # top row: upper‐tail
  upper_tale1 + upper_tale2 + upper_joint +
    plot_layout(ncol = 3)
) / (
  # bottom row: lower‐tail
  lower_tale1   + lower_tale2 + lower_joint   +
    plot_layout(ncol = 3)
) + 
  plot_layout(
    heights = c(1, 1),    # equal heights for both rows
    guides  = "collect"   # collect a single legend
  ) & 
  theme(legend.position = "bottom")

# draw it
print(combined_qq_all)


###############################################################
### combine all plots in one:
#— set up the 2×3 layout with outer margins for shared labels
par(
  mfrow = c(2,3),
  oma   = c(5, 5, 4, 8),   # bottom, left, top, right outer margins
  mar   = c(2, 2, 1, 1),    # small panel margins
  xpd   = FALSE 
)
# define your colors & set_labels once:
colors     <- c("goldenrod", "steelblue", "seagreen3", "purple2", "tomato")
set_labels <- c("Log","Bilog","HR","CT","MEGPD")

## 1) upper_tail1
plot(
  sim_data_df_list1_up[[1]]$Theoretical,
  sim_data_df_list1_up[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab = "upper tale",        # no per-panel axes
  main  = "Ammerzoden",
  xlim  = c(0, max(sim_data_df1_up$Theoretical)),
  ylim  = c(0, max(sim_data_df1_up$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size,
  abline(0,1, lty=2, col="gray50")
)

for(i in 2:length(sim_data_df_list1_up)){
  points(
    sim_data_df_list1_up[[i]]$Theoretical,
    sim_data_df_list1_up[[i]]$Observed,
    col = colors[i], pch=20
  )
}


## 2) upper_tail2
plot(
  sim_data_df_list2_up[[1]]$Theoretical,
  sim_data_df_list2_up[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab  = "",
  main  = "Zaltbommel",
  xlim  = c(0, max(sim_data_df2_up$Theoretical)),
  ylim  = c(0, max(sim_data_df2_up$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list2_up)){
  points(
    sim_data_df_list2_up[[i]]$Theoretical,
    sim_data_df_list2_up[[i]]$Observed,
    col = colors[i], pch=20
  )
}

## 3) upper_joint
plot(
  sim_data_df_list3_up[[1]]$Theoretical,
  sim_data_df_list3_up[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "",  ylab  = "",
  main  = "Ammerzoden + Zaltbommel",
  xlim  = c(0, max(sim_data_df3_up$Theoretical)),
  ylim  = c(0, max(sim_data_df3_up$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list3_up)){
  points(
    sim_data_df_list3_up[[i]]$Theoretical,
    sim_data_df_list3_up[[i]]$Observed,
    col = colors[i], pch=20
  )
}


## 4) lower_tail1
plot(
  sim_data_df_list1_low[[1]]$Theoretical,
  sim_data_df_list1_low[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab  = "lower tale",
  xlim  = c(0, max(sim_data_df1_low$Theoretical)),
  ylim  = c(0, max(sim_data_df1_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list1_low)){
  points(
    sim_data_df_list1_low[[i]]$Theoretical,
    sim_data_df_list1_low[[i]]$Observed,
    col = colors[i], pch=20
  )
}

## 5) lower_tail2
plot(
  sim_data_df_list2_low[[1]]$Theoretical,
  sim_data_df_list2_low[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab  = "",
  xlim  = c(0, max(sim_data_df2_low$Theoretical)),
  ylim  = c(0, max(sim_data_df2_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list2_low)){
  points(
    sim_data_df_list2_low[[i]]$Theoretical,
    sim_data_df_list2_low[[i]]$Observed,
    col = colors[i], pch=20
  )
}

## 6) lower_joint
plot(
  sim_data_df_list3_low[[1]]$Theoretical,
  sim_data_df_list3_low[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab  = "",
  xlim  = c(0, max(sim_data_df3_low$Theoretical)),
  ylim  = c(0, max(sim_data_df3_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list3_low)){
  points(
    sim_data_df_list3_low[[i]]$Theoretical,
    sim_data_df_list3_low[[i]]$Observed,
    col = colors[i], pch=20
  )
}

#— now add the *one* shared x‐ and y‐axis labels
mtext("Theoretical Quantiles",
      side = 1,        # bottom
      outer = TRUE,
      line  = 2.5,
      cex   = .85)

mtext("Sample Quantiles",
      side = 2,        # left
      outer = TRUE,
      line  = 2.5,
      cex   = .85)
# allow drawing in outer margin
# allow drawing in outer margin:
par(xpd = NA)

# pick your new colors:
legend_cols <- c("goldenrod", "steelblue", "seagreen3", "purple2", "tomato")

# place legend at mid‐height of the figure:
legend(
  x      = grconvertX(1.1, from = "npc", to = "user") + 0.05 * diff(par("usr")[1:2]),
  y      = grconvertY(1.2, from = "npc", to = "user"),
  legend = set_labels,
  col    = legend_cols,
  pch    = 20,
  bty    = "n",
  cex    = 0.8
)


par(
  mfrow = c(2,3),
  oma   = c(5, 5, 4, 8),   # bottom, left, top, right outer margins
  mar   = c(2, 2, 1, 1),    # small panel margins
  xpd   = FALSE 
)
# define your colors & set_labels once:
colors     <- c("goldenrod", "steelblue", "seagreen3", "purple2", "tomato")
set_labels <- c("Log","Bilog","HR","CT","MEGPD")

## 1) upper_tail1
plot(
  sim_data_df_list1_up[[1]]$Theoretical,
  sim_data_df_list1_up[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab = "upper tale",        # no per-panel axes
  main  = "Ammerzoden",
  xlim  = c(0, max(sim_data_df1_up$Theoretical)),
  ylim  = c(0, max(sim_data_df1_up$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size,
  abline(0,1, lty=2, col="gray50")
)

for(i in 2:length(sim_data_df_list1_up)){
  points(
    sim_data_df_list1_up[[i]]$Theoretical,
    sim_data_df_list1_up[[i]]$Observed,
    col = colors[i], pch=20
  )
}


## 2) upper_tail2
plot(
  sim_data_df_list2_up[[1]]$Theoretical,
  sim_data_df_list2_up[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab  = "",
  main  = "Zaltbommel",
  xlim  = c(0, max(sim_data_df2_up$Theoretical)),
  ylim  = c(0, max(sim_data_df2_up$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list2_up)){
  points(
    sim_data_df_list2_up[[i]]$Theoretical,
    sim_data_df_list2_up[[i]]$Observed,
    col = colors[i], pch=20
  )
}

## 3) upper_joint
plot(
  sim_data_df_list3_up[[1]]$Theoretical,
  sim_data_df_list3_up[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "",  ylab  = "",
  main  = "Ammerzoden + Zaltbommel",
  xlim  = c(0, max(sim_data_df3_up$Theoretical)),
  ylim  = c(0, max(sim_data_df3_up$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list3_up)){
  points(
    sim_data_df_list3_up[[i]]$Theoretical,
    sim_data_df_list3_up[[i]]$Observed,
    col = colors[i], pch=20
  )
}


## 4) lower_tail1
plot(
  sim_data_df_list1_low[[1]]$Theoretical,
  sim_data_df_list1_low[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab  = "lower tale",
  xlim  = c(0, max(sim_data_df1_low$Theoretical)),
  ylim  = c(0, max(sim_data_df1_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list1_low)){
  points(
    sim_data_df_list1_low[[i]]$Theoretical,
    sim_data_df_list1_low[[i]]$Observed,
    col = colors[i], pch=20
  )
}

## 5) lower_tail2
plot(
  sim_data_df_list2_low[[1]]$Theoretical,
  sim_data_df_list2_low[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab  = "",
  xlim  = c(0, max(sim_data_df2_low$Theoretical)),
  ylim  = c(0, max(sim_data_df2_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list2_low)){
  points(
    sim_data_df_list2_low[[i]]$Theoretical,
    sim_data_df_list2_low[[i]]$Observed,
    col = colors[i], pch=20
  )
}

## 6) lower_joint
plot(
  sim_data_df_list3_low[[1]]$Theoretical,
  sim_data_df_list3_low[[1]]$Observed,
  col   = colors[1], pch = 20,
  xlab  = "", ylab  = "",
  xlim  = c(0, max(sim_data_df3_low$Theoretical)),
  ylim  = c(0, max(sim_data_df3_low$Observed)),
  cex.axis = axis_tick_size, cex.lab = axis_label_size
)
abline(0,1, lty=2, col="gray50")
for(i in 2:length(sim_data_df_list3_low)){
  points(
    sim_data_df_list3_low[[i]]$Theoretical,
    sim_data_df_list3_low[[i]]$Observed,
    col = colors[i], pch=20
  )
}

#— now add the *one* shared x‐ and y‐axis labels
mtext("Theoretical Quantiles",
      side = 1,        # bottom
      outer = TRUE,
      line  = 2.5,
      cex   = .85)

mtext("Sample Quantiles",
      side = 2,        # left
      outer = TRUE,
      line  = 2.5,
      cex   = .85)
# allow drawing in outer margin
# allow drawing in outer margin:
par(xpd = NA)

# pick your new colors:
legend_cols <- c("goldenrod", "steelblue", "seagreen3", "purple2", "tomato")

# place legend at mid‐height of the figure:
legend(
  x      = grconvertX(1.1, from = "npc", to = "user") + 0.05 * diff(par("usr")[1:2]),
  y      = grconvertY(1.2, from = "npc", to = "user"),
  legend = set_labels,
  col    = legend_cols,
  pch    = 20,
  bty    = "n",
  cex    = 0.8
)

#################################

# 1) set up the 2×3 layout with shared outer margins
par(
  mfrow = c(2,3),
  oma   = c(5, 5, 4, 8),   # bottom, left, top, right outer margins
  mar   = c(2, 2, 1, 1),   # small panel margins
  xpd   = FALSE            # keep clipping on for panel drawing
)

# your color / label definitions
colors     <- c("purple2", "goldenrod", "steelblue", "seagreen3", "tomato")
set_labels <- c("Log","Bilog","HR","CT","MEGPD")
legend_cols<- colors

# helper to draw one panel:
draw_panel <- function(th, ob, list_data, main_title) {
  plot(
    th, ob, 
    col   = colors[1], pch = 20,
    xlab  = "", ylab = "",
    main  = main_title,
    xlim  = range(th, na.rm=TRUE),
    ylim  = range(ob, na.rm=TRUE),
    asp   = 1,                       # enforce 1:1 ratio
    panel.first = {                  # draw grid behind data
      grid(nx = NULL, ny = NULL, col="lightgrey", lty="dotted")
    },
    cex.axis = axis_tick_size, 
    cex.lab  = axis_label_size
  )
  abline(0,1, lty=2, col="gray50")
  for(i in 2:length(list_data)){
    points(
      list_data[[i]]$Theoretical,
      list_data[[i]]$Observed,
      col = colors[i], pch=20
    )
  }
}

# 2) draw the six QQ panels
draw_panel(
  sim_data_df_list1_up[[1]]$Theoretical,
  sim_data_df_list1_up[[1]]$Observed,
  sim_data_df_list1_up,
  "Giersbergen"
)
draw_panel(
  sim_data_df_list2_up[[1]]$Theoretical,
  sim_data_df_list2_up[[1]]$Observed,
  sim_data_df_list2_up,
  "Zaltbommel"
)
draw_panel(
  sim_data_df_list3_up[[1]]$Theoretical,
  sim_data_df_list3_up[[1]]$Observed,
  sim_data_df_list3_up,
  "Giersbergen + Zaltbommel"
)
draw_panel(
  sim_data_df_list1_low[[1]]$Theoretical,
  sim_data_df_list1_low[[1]]$Observed,
  sim_data_df_list1_low,
  ""
)
draw_panel(
  sim_data_df_list2_low[[1]]$Theoretical,
  sim_data_df_list2_low[[1]]$Observed,
  sim_data_df_list2_low,
  ""
)
draw_panel(
  sim_data_df_list3_low[[1]]$Theoretical,
  sim_data_df_list3_low[[1]]$Observed,
  sim_data_df_list3_low,
  ""
)

# 3) shared axis labels
mtext("Theoretical Quantiles",
      side  = 1,        # bottom
      outer = TRUE,
      line  = 2.5,
      cex = .8)
mtext("Sample Quantiles",
      side  = 2,        # left
      outer = TRUE,
      line  = 2.5,
      cex= .8)

# 4) single legend in outer right margin
par(xpd = NA)  # allow drawing in outer margin
# place legend at mid‐height of the figure:
legend(
  x      = grconvertX(1.1, from = "npc", to = "user") + 0.05 * diff(par("usr")[1:2]),
  y      = grconvertY(1.2, from = "npc", to = "user"),
  legend = set_labels,
  col    = legend_cols,
  pch    = 20,
  bty    = "n",
  cex    = 0.8
)
#####
#####save the plots:

# open a PNG device
png(
  filename = "qq_panels_15_23.png",
  width    = 8,      # width in inches
  height   = 6,      # height in inches
  units    = "in",
  res      = 300     # resolution in dpi
)

# --- your plotting code goes here ---
par(
  mfrow = c(2,3),
  oma   = c(5, 5, 4, 8),
  mar   = c(2, 2, 1, 1),
  xpd   = FALSE
)

colors     <- c("purple2","goldenrod","steelblue","seagreen3","tomato")
set_labels <- c("Log","Bilog","HR","CT","MEGPD")
legend_cols<- colors

draw_panel <- function(th, ob, list_data, main_title) {
  plot(
    th, ob,
    col   = colors[1], pch = 20,
    xlab  = "", ylab = "",
    main  = main_title,
    xlim  = range(th, na.rm=TRUE),
    ylim  = range(ob, na.rm=TRUE),
    asp   = 1,
    panel.first = {
      grid(nx = NULL, ny = NULL, col="lightgrey", lty="dotted")
    },
    cex.axis = axis_tick_size,
    cex.lab  = axis_label_size
  )
  abline(0,1, lty=2, col="gray50")
  for(i in 2:length(list_data)){
    points(
      list_data[[i]]$Theoretical,
      list_data[[i]]$Observed,
      col = colors[i], pch = 20
    )
  }
}

# draw the six panels
draw_panel(sim_data_df_list1_up[[1]]$Theoretical, sim_data_df_list1_up[[1]]$Observed, sim_data_df_list1_up, "Ammerzoden")
draw_panel(sim_data_df_list2_up[[1]]$Theoretical, sim_data_df_list2_up[[1]]$Observed, sim_data_df_list2_up, "Zaltbommel")
draw_panel(sim_data_df_list3_up[[1]]$Theoretical, sim_data_df_list3_up[[1]]$Observed, sim_data_df_list3_up, "Ammerzoden + Zaltbommel")
draw_panel(sim_data_df_list1_low[[1]]$Theoretical, sim_data_df_list1_low[[1]]$Observed, sim_data_df_list1_low, "")
draw_panel(sim_data_df_list2_low[[1]]$Theoretical, sim_data_df_list2_low[[1]]$Observed, sim_data_df_list2_low, "")
draw_panel(sim_data_df_list3_low[[1]]$Theoretical, sim_data_df_list3_low[[1]]$Observed, sim_data_df_list3_low, "")

# shared labels
mtext("Theoretical Quantiles", side=1, outer=TRUE, line=2.5, cex=.8)
mtext("Sample Quantiles",      side=2, outer=TRUE, line=2.5, cex=.8)

# legend
par(xpd=NA)
legend(
  x      = grconvertX(1.1, "npc", "user") + .05 * diff(par("usr")[1:2]),
  y      = grconvertY(1.2, "npc", "user"),
  legend = set_labels,
  col    = legend_cols,
  pch    = 20,
  bty    = "n",
  cex    = 0.8
)
# --- end plotting code ---

dev.off()
