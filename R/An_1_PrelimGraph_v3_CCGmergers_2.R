# AN1 prelimggraph

rm(list=ls())
gc()
# Load packages
library(tidyverse)
library(patchwork)
library(openxlsx)

source("R/config.R")

# Load data
load(file.path(data_dir, "Metadata_Indicators_ICB.RData"))
load(file.path(data_out, "icb_combined_mergers.RData"))


## ########################################################################
# Dx waits
## All dx waits - Create plots (6+ weeks, natural and log scales)

# Filter for both 6+ and 13+ week waiting times
icb_combined_mergers_alldx <- icb_combined_mergers %>%
  mutate(
    NHSCode_PostMerge = factor(NHSCode_PostMerge),
    IndicatorName = factor(IndicatorName),
    Value = Count / Denominator,
    Value_pct = Value * 100,
    # Extract waiting period (6 or 13 weeks)
    WaitingPeriod = case_when(
      str_detect(IndicatorName, "Number waiting 6")  ~ "6+ weeks",
      str_detect(IndicatorName, "Number waiting 13") ~ "13+ weeks",
      TRUE ~ NA_character_
    )
   
  )

# Create cleaned test names for facet labels with specific ordering
icb_combined_mergers_alldx <- icb_combined_mergers_alldx %>%
  mutate(TestName = str_extract(IndicatorName, "(?<=\\().*?(?=\\))") %>%
           str_trim() %>%
           str_to_title() %>%
           str_replace("^Ct$", "CT") %>%
           str_replace("^Mri$", "MRI") %>%
           str_replace("Flexi.*", "Flexible Sigmoidoscopy")
  ) %>%
  mutate(TestName = factor(TestName, levels = c(
    "MRI",
    "Gastroscopy",
    "CT",
    "Flexible Sigmoidoscopy",
    "Non Obstetric Ultrasound",
    "Colonoscopy",
    "Echocardiography",
    "Cystoscopy"
  )))

# Filter just to 6 weeks for this analysis
icb_combined_mergers_alldx <- icb_combined_mergers_alldx |>
  filter(str_detect(IndicatorName, "Number waiting 6"))


## ########################################################################
## Create plot

# Function to create summary stats and plot
create_waiting_plot <- function(data, waiting_period, use_log_scale = FALSE) {
  
  # Filter for specific waiting period
  plot_data <- data %>% filter(WaitingPeriod == waiting_period)
  
  # For log scale, filter out zero/negative values
  if (use_log_scale) {
    plot_data <- plot_data %>% filter(Value_pct > 0)
  }
  
  # Calculate sub-ICB distribution summary stats (median + IQR)
  summary_stats <- plot_data %>%
    group_by(Date, TestName) %>%
    summarise(
      median_val = median(Value_pct, na.rm = TRUE),
      q25        = quantile(Value_pct, 0.25, na.rm = TRUE),
      q75        = quantile(Value_pct, 0.75, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Calculate national total proportion (sum counts / sum denominators across all sub-ICBs)
  national_stats <- plot_data %>%
    group_by(Date, TestName) %>%
    summarise(
      national_count       = sum(Count, na.rm = TRUE),
      national_denominator = sum(Denominator, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(national_prob = national_count / national_denominator * 100)
  
  # Join for plotting
  combined_stats <- summary_stats %>%
    left_join(national_stats, by = c("Date", "TestName"))
  
  # Create base plot
  plot <- ggplot() +
    geom_line(data = plot_data,
              aes(x = Date, y = Value_pct, group = interaction(NHSCode_PostMerge, IndicatorName)),
              size = 0.3, alpha = 0.1, color = "gray70") +
    # IQR ribbon (sub-ICB spread)
    geom_ribbon(data = combined_stats,
                aes(x = Date, ymin = q25, ymax = q75),
                fill = "steelblue", alpha = 0.3) +
    # Median sub-ICB line
    geom_line(data = combined_stats,
              aes(x = Date, y = median_val, colour = "Median local area"),
              size = 0.8) +
    # National total proportion line
    geom_line(data = combined_stats,
              aes(x = Date, y = national_prob, colour = "National total"),
              size = 0.8) +
    scale_colour_manual(
      name   = NULL,
      values = c("Median local area" = "steelblue", "National total" = "black")
    ) +
    facet_wrap(~TestName, scales = "free_y", ncol = 2) +
    scale_x_date(
      date_labels = "%b %Y",
      date_breaks = "1 year"
    ) +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1),
      plot.title      = element_text(face = "bold", size = 14),
      strip.text      = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
  
  # Add scale and labels
  if (use_log_scale) {
    plot <- plot +
      scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))
      ) +
      labs(
        title    = paste("Percentage of patients waiting", waiting_period, "for a test, by local area (Log Scale)"),
        subtitle = "Individual local areas (grey), median local area (blue), national total (black), local area IQR (blue ribbon)",
        x        = "Time Period",
        y        = "% (log scale)"
      )
  } else {
    plot <- plot +
      scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0, 100)) +
      labs(
        title    = paste("Percentage of patients waiting", waiting_period, "for a test, by local area"),
        subtitle = "Individual local areas (grey), median local area (blue), national total (black), local area IQR (blue ribbon)",
        x        = "Time Period",
        y        = "%"
      )
  }
  
  return(plot)
}

# Create plots
trend_plot_6weeks_natural <- create_waiting_plot(icb_combined_mergers_alldx, "6+ weeks", use_log_scale = FALSE)
trend_plot_6weeks_log     <- create_waiting_plot(icb_combined_mergers_alldx, "6+ weeks", use_log_scale = TRUE)

# Display
trend_plot_6weeks_natural
trend_plot_6weeks_log

# Save plots
ggsave(file.path(output,"trend_plot_6weeks_natural.png"), plot = trend_plot_6weeks_natural, width = 3000, height = 4000, units = "px")
ggsave(file.path(output,"trend_plot_6weeks_log.png"), plot = trend_plot_6weeks_log, width = 3000, height = 4000, units = "px")

# Save dataset for reuse in later code
save(icb_combined_mergers_alldx, file = file.path(data_out,"icb_combined_mergers_alldx.RData"))


## #####################################################################
# APPENDIX TABLE — all months from January 2018
## #####################################################################

# Prepare data with YearMonth column
icb_combined_mergers_dx_extract <- icb_combined_mergers_alldx %>%
  mutate(YearMonth = format(Date, "%Y-%m"))

# Compute national totals and sub-ICB distribution per test and month
appendix_table <- icb_combined_mergers_dx_extract %>%
  group_by(TestName, YearMonth) %>%
  summarise(
    national_count       = as.integer(round(sum(Count, na.rm = TRUE))),
    national_denominator = as.integer(round(sum(Denominator, na.rm = TRUE))),
    n_icbs               = sum(!is.na(Value_pct) & Denominator > 0),
    national_prob        = round(national_count / national_denominator * 100, 2),
    median_pct           = round(median(Value_pct, na.rm = TRUE), 2),
    q25_pct              = round(quantile(Value_pct, 0.25, na.rm = TRUE), 2),
    q75_pct              = round(quantile(Value_pct, 0.75, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  # Order tests
  mutate(TestName_order = case_when(
    TestName == "MRI"                      ~ 1,
    TestName == "CT"                       ~ 2,
    TestName == "Non Obstetric Ultrasound" ~ 3,
    TestName == "Echocardiography"         ~ 4,
    TestName == "Gastroscopy"              ~ 5,
    TestName == "Flexible Sigmoidoscopy"   ~ 6,
    TestName == "Colonoscopy"              ~ 7,
    TestName == "Cystoscopy"               ~ 8,
    TRUE ~ 99
  )) %>%
  arrange(TestName_order, YearMonth) %>%
  select(
    TestName, YearMonth,
    national_count, national_denominator, national_prob,
    median_pct, q25_pct, q75_pct, n_icbs
  ) %>%
  rename(
    "Diagnostic test name"              = TestName,
    "Year-Month"                        = YearMonth,
    "National count"                    = national_count,
    "National denominator"              = national_denominator,
    "National % (waiting ≥6 weeks)"     = national_prob,
    "Median local area % (waiting ≥6 wks)" = median_pct,
    "25th percentile local area %"         = q25_pct,
    "75th percentile local area %"         = q75_pct,
    "N local areas"                        = n_icbs
  )

save(appendix_table, file = file.path(output,"appendix_table.RData"))


## #####################################################################
# MAIN RESULTS TABLE — 6-monthly snapshots, long format
## #####################################################################

latest_month <- max(icb_combined_mergers_dx_extract$Date, na.rm = TRUE)
start_date   <- ymd("2018-01-01")

# Generate 6-monthly dates back from latest month
six_month_dates  <- seq(from = latest_month, to = start_date, by = "-6 months")
six_month_labels <- format(six_month_dates, "%Y-%m")

# Filter appendix to 6-monthly snapshots, keep long format
main_results_table <- appendix_table %>%
  filter(`Year-Month` %in% six_month_labels)

save(main_results_table, file = file.path(output,"main_results_table.RData"))


## #####################################################################
# MAIN RESULTS TABLE 2 — latest month vs same month previous year, wide
## #####################################################################

prev_year_month   <- latest_month %m-% years(1)
comparison_labels <- format(c(latest_month, prev_year_month), "%Y-%m")

main_results_table2 <- appendix_table %>%
  filter(`Year-Month` %in% comparison_labels) %>%
  pivot_wider(
    id_cols     = `Diagnostic test name`,
    names_from  = `Year-Month`,
    values_from = c(`National count`, `National denominator`,
                    `National % (waiting ≥6 weeks)`,
                    `Median local area % (waiting ≥6 wks)`,
                    `25th percentile local area %`,
                    `75th percentile local area %`),
    names_glue  = "{.value} ({`Year-Month`})"
  )

# Order columns: previous year first, then latest month
col_order2 <- c(
  "Diagnostic test name",
  grep(paste0("\\(", comparison_labels[2], "\\)"), names(main_results_table2), value = TRUE),
  grep(paste0("\\(", comparison_labels[1], "\\)"), names(main_results_table2), value = TRUE)
)
main_results_table2 <- main_results_table2 %>%
  select(any_of(col_order2))

save(main_results_table2, file = file.path(output,"main_results_table2.RData"))




