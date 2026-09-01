# # parameters for running sensitivity models

# Briefly load data to get latest months
# Load data
load(file.path(data_out, "icb_combined_mergers_alldx.RData"))

# Select 4 months ago vs 12 months before that
latest_month    <- floor_date(max(icb_combined_mergers_alldx$Date, na.rm = TRUE), "month") %m-% months(5)
prev_year_month <- latest_month %m-% years(1)

# Remove data
rm(icb_combined_mergers_alldx)

# File path for saving
run_label <- paste0("sens_", format(latest_month, "%Y-%m"))