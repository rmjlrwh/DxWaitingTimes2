# parameters for running main models

# Briefly load data to get latest months
# Load data
load(file.path(data_out, "icb_combined_mergers_alldx.RData"))

# Latest month vs 12 months before
latest_month    <- floor_date(max(icb_combined_mergers_alldx$Date, na.rm = TRUE), "month")
prev_year_month <- latest_month %m-% years(1)

# Remove data
rm(icb_combined_mergers_alldx)

# File path for saving
run_label <- "main"
