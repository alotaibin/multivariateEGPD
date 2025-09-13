##############################################################################
## Downloading and Processing ECAD Rainfall Data
##############################################################################

options(timeout = 2000)

# Download raw data
download.file(
  "https://knmi-ecad-assets-prd.s3.amazonaws.com/download/ECA_blend_rr.zip",
  destfile = "ECA_blend_rr.zip",
  mode     = "wb",
  method   = "libcurl"
)

# Unzip into intermediates/
unlink("intermediates/ECA_blend_rr", recursive = TRUE)
unzip("ECA_blend_rr.zip", exdir = "intermediates/ECA_blend_rr")

##############################################################################
## Build & filter ECAD NL rainfall data, then extract 3 stations
##############################################################################

# Dependencies
library(imputeTS)
library(xts)

# Paths and parameters
out_dir    <- "intermediates/ECA_blend_rr"
meta_f     <- file.path(out_dir, "stations.txt")
SY         <- 1999
EY         <- 2024
sel_months <- c(10,11,12,1,2)
PP         <- 0.90

# Read & filter metadata
info_all <- read.csv(meta_f, skip=17, stringsAsFactors=FALSE)
info_all$CN <- gsub(" ", "", info_all$CN, fixed=TRUE)
info_nl <- subset(info_all, CN == "NL")

# Assemble filenames
make_fname <- function(id) sprintf("RR_STAID0%05d.txt", id)
fn_raw <- file.path(out_dir, sapply(info_nl$STAID, make_fname))

# Determine global date range
min_year <- integer(length(fn_raw))
max_year <- integer(length(fn_raw))
for (i in seq_along(fn_raw)) {
  d <- read.table(fn_raw[i], header=TRUE, sep=",", skip=20,
                  colClasses=c("character","numeric","integer"))
  dates <- as.character(d$DATE)
  min_year[i] <- as.integer(substr(dates[1],     1,4))
  max_year[i] <- as.integer(substr(dates[length(dates)],1,4))
}
mindate   <- as.Date(paste0(min(min_year),"0101"),"%Y%m%d")
maxdate   <- as.Date(paste0(max(max_year),"1231"),"%Y%m%d")
all_dates <- seq(mindate, maxdate, by="day")

# Raw daily matrix
raw_mat <- matrix(NA_real_, length(all_dates), length(fn_raw))
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

# Write intermediates
write.csv(cbind(info_nl, minyear=min_year, maxyear=max_year),
          "intermediates/info-stations-NL-RR.txt", row.names=FALSE)
write.csv(data.frame(day=all_dates, raw_mat),
          "intermediates/data-stations-NL-RR.txt", row.names=FALSE)

##############################################################################
## Temporal + spatial filtering
##############################################################################

dat   <- read.csv("intermediates/data-stations-NL-RR.txt", stringsAsFactors=FALSE)
dat$day <- as.Date(dat$day)
dat   <- subset(dat, format(day,"%Y") %in% as.character(SY:EY))

info2 <- read.csv("intermediates/info-stations-NL-RR.txt", stringsAsFactors=FALSE)

# Convert coords to decimal degrees
convert.lat <- function(coord) {
  d1 <- as.numeric(substr(coord,1,3))
  d2 <- as.numeric(substr(coord,5,6))
  d3 <- as.numeric(substr(coord,8,9))
  dec <- abs(d1) + d2/60 + d3/3600
  if (d1 < 0) dec <- -dec
  dec
}
convert.long <- function(coord) {
  d1 <- as.numeric(substr(coord,1,4))
  d2 <- as.numeric(substr(coord,6,7))
  d3 <- as.numeric(substr(coord,9,10))
  dec <- abs(d1) + d2/60 + d3/3600
  if (d1 < 0) dec <- -dec
  dec
}
info2$dec.lat <- sapply(as.character(info2$LAT), convert.lat)
info2$dec.lon <- sapply(as.character(info2$LON), convert.long)

# Spatial box: Noord-Brabant
lat_min <- 51.25; lat_max <- 51.85
lon_min <- 4.43;  lon_max <- 5.48
sel_geo <- with(info2,
                dec.lat > lat_min & dec.lat < lat_max &
                dec.lon > lon_min & dec.lon < lon_max)

keep_ids <- info2$STAID[sel_geo]
cols_geo <- paste0("S", keep_ids)
dat2     <- dat[, c("day", intersect(cols_geo, names(dat)))]

# Filter missingness
frac_ok <- colSums(!is.na(dat2[,-1])) / nrow(dat2)
keep2   <- names(frac_ok)[frac_ok > PP]
dat3    <- dat2[, c("day", keep2)]

# Interpolate missing
for (nm in keep2) dat3[[nm]] <- na_interpolation(dat3[[nm]])

# Restrict to fall–winter
xt       <- xts(dat3[,-1], order.by=dat3$day)
mon      <- as.integer(format(index(xt), "%m"))
xt2      <- xt[mon %in% sel_months, ]
mat_season <- coredata(xt2)
colnames(mat_season) <- colnames(xt2)

##############################################################################
## Select 3 stations and build pair matrices
##############################################################################

stations <- info2[ sel_geo, ]
sel_cols <- c(11, 23, 15)  # chosen station indices
final_mat  <- mat_season[, sel_cols]

# Save full processed matrix
write.csv(final_mat, "intermediates/RR_three_stations_NL_fall-winter.csv", row.names=FALSE)

# Save pairwise scaled versions
make_and_save_scaled_pair <- function(i, j, name) {
  m <- mat_season[, c(i, j), drop=FALSE]
  keep <- (m[,1] != 0) & (m[,2] != 0)
  m <- m[keep, , drop=FALSE]
  sds <- apply(m, 2, sd, na.rm=TRUE)
  m_scaled <- sweep(m, 2, sds, `/`)
  colnames(m_scaled) <- paste0("S", stations$STAID[c(i,j)])
  write.csv(as.data.frame(m_scaled),
            file = sprintf("intermediates/pair_%s_matrix_scaled.csv", name),
            row.names=FALSE)
}
make_and_save_scaled_pair(sel_cols[1], sel_cols[2], "11_23")
make_and_save_scaled_pair(sel_cols[1], sel_cols[3], "11_15")
make_and_save_scaled_pair(sel_cols[2], sel_cols[3], "15_23")

##############################################################################
## Save station map (for reproducibility)
##############################################################################
coords_p <- cbind(stations$dec.lon, stations$dec.lat)
png("Figures/Application/stations_map.png", width=6, height=5, units="in", res=300)
plot(coords_p, pch=16, col="gray", xlab="Longitude", ylab="Latitude",
     main="Filtered NL Stations\nChosen in Red")
points(coords_p[sel_cols, ], pch=16, col="red")
dev.off()
