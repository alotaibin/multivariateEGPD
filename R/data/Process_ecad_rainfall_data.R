##############################################################################
## Downloading and Processing ECAD Rainfall Data:
##############################################################################

options(timeout = 2000)

download.file(
  "https://knmi-ecad-assets-prd.s3.amazonaws.com/download/ECA_blend_rr.zip",
  destfile = "ECA_blend_rr.zip",
  mode     = "wb",
  method   = "libcurl"
)

# Unzip into “ECA_blend_rr/”
unlink("ECA_blend_rr", recursive = TRUE)               # remove any old folder
unzip("ECA_blend_rr.zip", exdir = "ECA_blend_rr")      # extract all station files

# Verifying the station files
files <- list.files("ECA_blend_rr_data/ECA_blend_rr", recursive = TRUE)
length(grep("^RR_STAID0.*\\.txt$", files))            # 16320 should be ~17125
#head(grep("^RR_STAID0.*\\.txt$", files, value = TRUE), 10)

##############################################################################
## Build & filter ECAD NL rainfall data, then extract 3 stations
##############################################################################

# 0. Dependencies
library(imputeTS)
library(xts)

# 1. Paths and parameters
out_dir  <- "ECA_blend_rr_data/ECA_blend_rr"
meta_f   <- file.path(out_dir, "stations.txt") #Carlo: ECA_blend_station_rr.txt

# Season & missing-data cutoff
SY       <- 1999
EY       <- 2024
sel_months <- c(10,11,12,1,2)
PP         <- 0.90

# 2. Read & filter station metadata
info_all <- read.csv(meta_f, skip=17, stringsAsFactors=FALSE)
info_all$CN <- gsub(" ", "", info_all$CN, fixed=TRUE)
info_nl <- subset(info_all, CN == "NL")

# 3. Build file list for each station
make_fname <- function(id) sprintf("RR_STAID0%05d.txt", id)
fn_raw <- file.path(out_dir, sapply(info_nl$STAID, make_fname))

# 4. Determine global date range
min_year <- integer(length(fn_raw))
max_year <- integer(length(fn_raw))
for (i in seq_along(fn_raw)) {
  d <- read.table(fn_raw[i], header=TRUE, sep=",", skip=20,
                  colClasses=c("character","numeric","integer"))
  dates <- as.character(d$DATE)
  min_year[i] <- as.integer(substr(dates[1],     1,4))
  max_year[i] <- as.integer(substr(dates[length(dates)],1,4))
}
mindate <- as.Date(paste0(min(min_year),"0101"),"%Y%m%d")
maxdate <- as.Date(paste0(max(max_year),"1231"),"%Y%m%d")
all_dates <- seq(mindate, maxdate, by="day")
n_days    <- length(all_dates) #65744
n_sites   <- length(fn_raw) #549

# 5. Assemble the raw daily matrix
raw_mat <- matrix(NA_real_, n_days, n_sites)
for (i in seq_along(fn_raw)) {
  d <- read.table(fn_raw[i], header=TRUE, sep=",", skip=20,
                  colClasses=c("character","numeric","integer"))
  y <- d$RR
  y[d$Q_RR != 0] <- NA
  day0 <- as.Date(as.character(d$DATE[1]), "%Y%m%d")
  off  <- as.integer(day0 - mindate)
  raw_mat[(off+1):(off+length(y)), i] <- y
}
colnames(raw_mat) <- paste0("S", info_nl$STAID)
rownames(raw_mat) <- as.character(all_dates)

# 6. (Optional) write intermediate CSVs
write.csv(cbind(info_nl, minyear=min_year, maxyear=max_year),
          "info-stations-NL-RR.txt", row.names=FALSE)
write.csv(data.frame(day=all_dates, raw_mat),
          "data-stations-NL-RR.txt", row.names=FALSE)

# 7. Read back & temporal filter
dat <- read.csv("ECA_blend_rr_data/data-stations-NL-RR.txt", stringsAsFactors=FALSE)
dat$day <- as.Date(dat$day)
dat <- subset(dat, format(day,"%Y") %in% as.character(SY:EY))  #550 stations

# 8. Spatial filter: Noord-Brabant box
info2 <- read.csv("ECA_blend_rr_data/info-stations-NL-RR.txt", stringsAsFactors=FALSE)
print(names(info2)) #"STAID"   "STANAME" "CN"  "LAT"  "LON" "HGHT" "minyear" "maxyear"

#Convert the DMS LAT/LON into decimal degrees
convert.lat <- function(coord) {
  d1 <- as.numeric(substr(coord,  1, 3))
  d2 <- as.numeric(substr(coord,  5, 6))
  d3 <- as.numeric(substr(coord,  8, 9))
  dec <- abs(d1) + d2/60 + d3/3600
  if (d1 < 0) dec <- -dec
  dec
}
convert.long <- function(coord) {
  d1 <- as.numeric(substr(coord,  1, 4))
  d2 <- as.numeric(substr(coord,  6, 7))
  d3 <- as.numeric(substr(coord,  9,10))
  dec <- abs(d1) + d2/60 + d3/3600
  if (d1 < 0) dec <- -dec
  dec
}

info2$dec.lat  <- sapply(as.character(info2$LAT),  convert.lat)
info2$dec.lon  <- sapply(as.character(info2$LON),  convert.long)

#define our box
lat_min <- 51.25; lat_max <- 51.85
lon_min <- 4.43;  lon_max <- 5.48

sel_geo <- with(info2,
                dec.lat > lat_min & dec.lat < lat_max &
                  dec.lon > lon_min & dec.lon < lon_max
)

keep_ids <- info2$STAID[sel_geo]
# check
#print(keep_ids) #41 stations

# build the column names in dat
cols_geo  <- paste0("S", keep_ids)
cols_exist <- intersect(c("day", cols_geo), names(dat))
dat2       <- dat[, cols_exist]  # 42 stations
#View(dat2)

#Subset your dat to those station‐columns
PP    <- 0.90
frac_ok <- colSums(!is.na(dat2[,-1])) / nrow(dat2)
keep2   <- names(frac_ok)[frac_ok > PP]
dat3    <- dat2[, c("day", keep2)] #32 stations
#View(dat3)

#Interpolate remaining NAs
for (nm in keep2) {
  dat3[[nm]] <- na_interpolation(dat3[[nm]])
}

#Subset to fall–winter months (Oct–Feb)
sel_months <- c(10,11,12,1,2)

xt  <- xts(dat3[,-1], order.by = dat3$day)
mon <- as.integer(format(index(xt), "%m"))

xt2 <- xt[ mon %in% sel_months, ]
mat_season <- coredata(xt2)
colnames(mat_season) <- colnames(xt2)


# Build the `stations` metadata to match mat_season’s columns

# Start from the spatially‐filtered metadata
info_spatial <- info2[ sel_geo, ]

# Then apply your missingness filter: keep2 are names like "S2357","S2417",…
#    Strip "S" to match the integer STAID column
sta_keep_ids <- as.integer(sub("^S", "", keep2))

info_filt <- info_spatial[ info_spatial$STAID %in% sta_keep_ids, ]

# 3. Now reorder info_filt so that its rows are in the same order as mat_season’s columns
#    (colnames(mat_season) are exactly keep2, in the same order)
stations <- info_filt[ match(sta_keep_ids, info_filt$STAID), ]

# 4. Quick sanity check: column names vs STAID
stopifnot( all( paste0("S", stations$STAID) == colnames(mat_season) ) )

# 5. Now stations has one row per column of mat_season, with all metadata (incl. dec.lat, dec.lon)
View(stations)

# Extract stations 11, 15 & 27
sel_cols   <- c(11, 15, 27)
final_mat  <- mat_season[, sel_cols]
final_info <- info2[ sel_geo, ][ sel_cols, ]

print(final_info)
dim(final_mat)   # days × 3

View(final_mat)

#################################################################################
#################################################################################
# --- 1) Rebuild `stations` to match mat_season --------------------------

# info2   : raw NL metadata (with dec.lat, dec.lon) from step 8
# sel_geo : logical vector of spatially‐filtered stations
# keep2   : character vector of columns “Sxxxx” that survived missing‐data filter
# mat_season: the fall–winter matrix (cols = keep2, rows = days)

# a) Spatially filtered metadata
info_spatial <- info2[ sel_geo, ]

# b) Metadata after missingness filter
sta_ids      <- as.integer(sub("^S", "", keep2))
info_filt    <- info_spatial[ info_spatial$STAID %in% sta_ids, ]

# c) Reorder to line up with mat_season columns
stations     <- info_filt[ match(sta_ids, info_filt$STAID), ]

# Sanity check
stopifnot(all(paste0("S", stations$STAID) == colnames(mat_season)))

# --- 2) Define planar‐projection function -------------------------------

lonlat.to.planar <- function(lon.lat, miles = FALSE) {
  rdistearth <- function(loc1, loc2 = loc1, miles = FALSE) {
    R   <- if (miles) 3963.34 else 6378.388
    rad <- pi/180
    lat1 <- loc1[,2]*rad; lon1 <- loc1[,1]*rad
    lat2 <- loc2[,2]*rad; lon2 <- loc2[,1]*rad
    PP1  <- cbind(cos(lat1)*cos(lon1), cos(lat1)*sin(lon1), sin(lat1))
    PP2  <- cbind(cos(lat2)*cos(lon2), cos(lat2)*sin(lon2), sin(lat2))
    cosd <- PP1 %*% t(PP2)
    pmax <- pmin(cosd, 1)
    R * acos(pmax)
  }
  x  <- lon.lat[,1]; y <- lon.lat[,2]
  mx <- mean(x); my <- mean(y)
  sy <- rdistearth(cbind(rep(mx,2), range(y)))[2,1]
  sx <- rdistearth(cbind(range(x), rep(my,2)))[2,1]
  facx <- sx/(max(x)-min(x)); facy <- sy/(max(y)-min(y))
  cbind((x-mx)*facx, (y-my)*facy)
}

# --- 3) Compute planar coords & mask chosen stations ---------------------

coords_raw <- cbind(stations$dec.lon, stations$dec.lat)
coords_p   <- lonlat.to.planar(coords_raw)

is_chosen  <- seq_len(nrow(coords_p)) %in% sel_cols

# --- 4) Plot all stations + labels, highlight chosen ---------------------

par(mfrow = c(1,1), mar = c(4,4,2,1))

plot(
  coords_p, type = "n",
  xlab = "Projected X", ylab = "Projected Y",
  main = "Filtered Stations (NL fall–winter)\nChosen in Red"
)

# a) all points in light gray
points(coords_p, pch = 16, cex = 0.6, col = "lightgray")

# b) label **every** station in dark gray
text(
  coords_p,
  labels = seq_len(nrow(coords_p)),
  col    = "darkgray",
  cex    = 0.7,
  pos    = 3
)

# c) overplot chosen points in red
points(coords_p[is_chosen, ], pch = 16, cex = 0.8, col = "red")

# d) re‐label chosen stations in bold red
text(
  coords_p[is_chosen, ],
  labels = sel_cols,
  col    = "red",
  font   = 2,
  cex    = 0.7,
  pos    = 3
)

# e) legend
legend(
  "topright",
  legend = c("Other stations", "Chosen stations"),
  pch    = 16,
  col    = c("lightgray", "red")
)
#################################################################################
#################################################################################
## save:
# Include the dates as a column (if rownames are dates)
out_df <- data.frame(
  final_mat,
  row.names = NULL
)

write.csv(
  out_df,
  file = "RR_three_stations_NL_fall-winter.csv",
  row.names = FALSE
)


#################################################################################
#################################################################################

# Your candidate station *positions*:
potential <- c(6, 25, 23, 28, 22, 11, 15, 27)
alldata <- mat_season
# Get the *names* of those columns in mat_season
candidate_names <- colnames(alldata)[ potential ]

#1) Plot all their time series
plot.ts(
  mat_season[, potential],
  ylab = "Rainfall (mm)",
  main = "Daily Rainfall for Candidate Stations (6, 25, 23, 28, 22, 11, 15, 27)"
)


# 2) Count zeros in each station
ts_mat <- mat_season[, potential]
zero_counts <- colSums(ts_mat == 0, na.rm = TRUE)
zero_df <- data.frame(
  station = candidate_names,
  zeros   = zero_counts,
  row.names = NULL
)
print(zero_df)


# 3) POT fit at the 95% threshold for each candidate
cat("### POT fits (shape & scale) at 95% threshold ###\n\n")
for (k in potential) {
  y <- as.numeric(alldata[, k])
  q <- quantile(y, probs = 0.95, na.rm = TRUE)
  fit <- evd::fpot(y, threshold = q)
  cat("Station", k, "-", stations[k, "STANAME"], "\n")
  print(fit$estimate)      # shape & scale
  cat("\n")
}

# 4) Chi-plots for each pair of candidates
cat("### Chi-plots for each pair ###\n")
pairs <- combn(potential, 2, simplify = FALSE)
par(mfrow = c(2, 2))
for (pr in pairs) {
  title <- paste0("Stations ", pr[1], " vs ", pr[2])
  evd::chiplot(alldata[, pr], which = 1, main1 = title)
}

# 5) Autocorrelation of the nonzero series
cat("### ACF of nonzero series ###\n\n")
# First build a matrix where zeros → NA
data_nz <- alldata
data_nz[data_nz == 0] <- NA

par(mfrow = c(2, 2))
for (k in potential) {
  acf(data_nz[, k], na.action = na.pass, main = paste("ACF Station", k))
}

#################################################################################
#################################################################################

### Choose three stations:

# 1. Extract the series for station 11 (first pick) and each candidate
base_series <- mat_season[, sel_cols[1]]          # e.g. column 11
cand_mat    <- mat_season[, potential]   # matrix of candidates

# 2. Compute Pearson correlations (use complete.obs to drop any NA pairs)
cors <- sapply(seq_along(potential), function(i) {
  cor(base_series, cand_mat[,i], use = "complete.obs")
})

# 3. Summarize
cand_ids <- colnames(mat_season)[potential]
cor_df   <- data.frame(
  idx    = potential,
  STAID  = cand_ids,
  cor_to_first = cors
)
print(cor_df)

# 4. Pick the candidate with *lowest* absolute correlation
best     <- cor_df$idx[ which.min(abs(cor_df$cor_to_first)) ]
best_id  <- cor_df$STAID[ cor_df$idx == best ]

cat("Replacing station", sel_cols[3], "with station", best, "(", best_id, ")\n")

# 5. Build your new sel vector
new_sel <- c(sel_cols[1], sel_cols[2], best)

#################################################################################
# Chosen stations: station_indx <- c(11, 23, 15)
# Helper to build, rescale, and save a pair matrix
make_and_save_scaled_pair <- function(i, j, name) {
  # 1) Extract the two‐column submatrix and drop zeros
  m <- mat_season[, c(i, j), drop=FALSE]
  keep <- (m[,1] != 0) & (m[,2] != 0)
  m <- m[keep, , drop=FALSE]
  # 2) Compute empirical SD of each column
  sds <- apply(m, 2, sd, na.rm = TRUE)
  # 3) Rescale
  m_scaled <- sweep(m, 2, sds, `/`)
  # 4) Save scaling factors for later (so you can “un‐scale” if needed)
  scale_df <- data.frame(
    station = paste0("S", stations$STAID[c(i,j)]),
    sd      = sds
  )
  write.csv(
    scale_df,
    file = sprintf("pair_%s_scales.csv", name),
    row.names = FALSE
  )
  # 5) Save the scaled matrix
  colnames(m_scaled) <- paste0("S", stations$STAID[c(i,j)])
  write.csv(
    as.data.frame(m_scaled),
    file = sprintf("pair_%s_matrix_scaled.csv", name),
    row.names = FALSE
  )
}
# Loop over your three pairs
for (nm in names(pairs)) {
  idx <- pairs[[nm]]
  make_and_save_scaled_pair(idx[1], idx[2], nm)
}
df_11_15_scaled <- read.csv("pair_11_15_matrix_scaled.csv", stringsAsFactors = FALSE)
df_11_23_scaled <- read.csv("pair_11_23_matrix_scaled.csv", stringsAsFactors = FALSE)
df_15_23_scaled <- read.csv("pair_15_23_matrix_scaled.csv", stringsAsFactors = FALSE)
dim(df_11_15_scaled)
dim(df_11_23_scaled)
dim(df_15_23_scaled)
mat_11_15_scaled <- as.matrix(df_11_15_scaled)
mat_11_23_scaled <- as.matrix(df_11_23_scaled)
mat_15_23_scaled <- as.matrix(df_15_23_scaled)
