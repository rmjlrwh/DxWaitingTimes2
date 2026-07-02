#AN2 models - FILE B: Format results, build plots, export to Excel

rm(list=ls())
gc()
# Load packages
library(tidyverse)
library(lme4)
library(patchwork)
library(openxlsx)

source("R/config.R")


## #####################################################################

# Define comparison months- Current Month vs Same month previous year
# This is done by referencing the params files
# params = main analsyis, params_sensitivity for the sensitivity analysis

# source("R/params.R")  # swap to params_sensitivity.R for the sensitivity run
 source("R/params_sensitivity.R")  # swap to params_sensitivity.R for the sensitivity run


## #####################################################################


# Load sub icb name file
load(file.path(data_out, "subicb.RData"))

# Load modelling outputs from file 1
load(file.path(data_out, paste0(run_label,"_AN2_model_outputs.RData")))

# Load saved models (needed for prediction plots)
all_models <- readRDS(file.path(data_out, paste0(run_label,"_all_models_auto.RData")))

# Load data (needed for histogram plots)
load(file.path(data_out,  "icb_combined_mergers_alldx.RData"))

## ####################################################################

# Define comparison months-  Month x vs Same month previous year
comparison_months <- c(latest_month, prev_year_month)
comparison_labels <- format(comparison_months, "%Y-%m")

## #####################################################################

# Rebuild data_2month for plots
data_2month <- icb_combined_mergers_alldx |>
  mutate(YearMonth = format(Date, "%Y-%m")) |>
  filter(YearMonth %in% comparison_labels) |>
  mutate(
    NHSCode_PostMerge = factor(NHSCode_PostMerge),
    Timeperiod        = factor(Timeperiod),
    YearMonth         = factor(YearMonth)
  )

all_tests <- c(imaging_tests, endoscopy_tests)


## #################################################################
## Format modelled prev vs latest comparison

modelled_results_table <- modelled_results |>
  pivot_wider(
    id_cols     = c(test),
    names_from  = YearMonth,
    values_from = c(
      mean_logit,
      variance_logit,
      sd_logit,
      sd_diff,
      sd_diff_lower_95,
      sd_diff_upper_95,
      sd_diff_pvalue,
      pct25,
      pct75
    ),
    names_glue  = "{.value} ({YearMonth})"
  )

# Order columns so previous year comes first
modelled_results_table <- modelled_results_table |>
  select(
    test,
    contains(paste0("(", comparison_labels[2], ")")),
    contains(paste0("(", comparison_labels[1], ")")),
    contains("sd_diff")
  )

# Keep just the sd diff for one time period and rename
modelled_results_table <- modelled_results_table |>
  select(-which(
    str_detect(names(modelled_results_table), paste0("(", comparison_labels[2], ")")) 
    & str_detect(names(modelled_results_table), "sd_diff")))


## #################################################################
## Format obs vs modelled tables

# OR difference
obs_vs_modelled <- obs_vs_modelled |>
  mutate(
    obs_pct75 = obs_pct75 / 100,
    obs_pct25 = obs_pct25 / 100,
    mod_pct75 = mod_pct75 / 100,
    mod_pct25 = mod_pct25 / 100,
  ) |>
  mutate(
    obs_or_diff = (obs_pct75/(1-obs_pct75)) / (obs_pct25/(1-obs_pct25)),
    mod_or_diff = (mod_pct75/(1-mod_pct75)) / (mod_pct25/(1-mod_pct25))
  ) |>
  mutate(
    obs_pct75 = obs_pct75 * 100,
    obs_pct25 = obs_pct25 * 100,
    mod_pct75 = mod_pct75 * 100,
    mod_pct25 = mod_pct25 * 100,
  ) |>
  mutate(
    across(ends_with("or_diff"), ~ round(.x, 2)),
    across(c(obs_pct25, obs_pct75, mod_pct25, mod_pct75), ~ round(.x, 1))
  )

# Shortned main results table
main_results_table <- obs_vs_modelled |>
  select(
    test,
    YearMonth,
    obs_pct25, obs_pct75, obs_or_diff,
    mod_pct25, mod_pct75, mod_or_diff
  )


## #####################################################################
# Build histogram plots for all tests, stored in named lists
# One list for prev year, one for latest month

plots_hist_prev <- list()
plots_hist_latest <- list()

for (test in all_tests) {
  
  dat <- data_2month %>%
    filter(TestName == test) %>%
    mutate(YearMonth = factor(as.character(YearMonth), levels = c(ref_prev, ref_latest))) %>%
    droplevels()
  
  # Create observed density data frame
  observed_test <- dat %>%
    mutate(
      probability = Count / Denominator,
      Timeperiod = case_when(
        YearMonth == ref_prev   ~ label_prev,
        YearMonth == ref_latest ~ label_latest)
    )
  
  # Filter modelled results to test and month of interest
  mod_row_prev <- modelled_results |> filter(test == !!test, YearMonth == ref_prev)
  mod_row_latest <- modelled_results |> filter(test == !!test, YearMonth == ref_latest)
  
  # Create modelled density for histogram, by test and previous vs latest month
  make_density_df <- function(mean_logit, sd) {
    
    # set up variable containg equally space points 1-1000 in proportion space 
    tibble(density_x = (1:1000) / 1000) |>
      mutate(
        
        # Convert density points onto log odds scale
        density_x_log = qlogis(density_x),
        
        # *calculate PDF for these points on log-odds scale using the intercept (mean_logit) and 
        # standard deviation of the random effect from the mixed model
        density_log   = dnorm(density_x_log, mean = mean_logit, sd = sd),
        
        # scale these values to account for the fact that points are 
        # differently spaced on the proportion scale
        density       = density_log *
          (exp(-2 * density_x_log) + (2 * exp(-density_x_log)) + 1) /
          exp(-density_x_log),
        
        # scale x-values and density by 100 to move from proportion to percentage
        density_x     = density_x * 100,
        density       = density / 100
      )
  }
  
  density_prev <- make_density_df(mod_row_prev$mean_logit, mod_row_prev$sd_logit)
  density_latest <- make_density_df(mod_row_latest$mean_logit, mod_row_latest$sd_logit)
  
  # Histogram for previous month
  plots_hist_prev[[test]] <- ggplot() +
    
    # Observed density
    geom_histogram(
      data     = observed_test |> filter(Timeperiod == label_prev),
      aes(
        x = probability * 100,
        y = after_stat(density),
        fill = "Observed"
      ),
      boundary = 0,
      binwidth = 4,
      colour   = "white"
    ) +
    
    # Modelled density
    geom_line(
      data = density_prev |> filter(density_x < 80),
      aes(
        x = density_x,
        y = density,
        colour = "Modelled"
      ),
      linewidth = 1
    ) +
    scale_fill_manual(values = c("Observed" = "lightgrey")) +
    scale_colour_manual(values = c("Modelled" = "blue")) +
    labs(
      title = test,
      x     = "% waiting >6 weeks",
      y     = "Density",
      fill  = "",
      colour = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
  # Histogram for latest month
  plots_hist_latest[[test]] <- ggplot() +
    geom_histogram(
      data     = observed_test |> filter(Timeperiod == label_latest),
      aes(
        x = probability * 100,
        y = after_stat(density),
        fill = "Observed"
      ),
      boundary = 0,
      binwidth = 4,
      colour   = "white"
    ) +
    geom_line(
      data = density_latest |> filter(density_x < 80),
      aes(
        x = density_x,
        y = density,
        colour = "Modelled"
      ),
      linewidth = 1
    ) +
    scale_fill_manual(values = c("Observed" = "lightgrey")) +
    scale_colour_manual(values = c("Modelled" = "blue")) +
    labs(
      title = test,
      x     = "% waiting >6 weeks",
      y     = "Density",
      fill  = "",
      colour = ""
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")
  
}


## #####################################################################
## Patchworks — PREV

patchwork_hist_imaging_prev <- wrap_plots(plots_hist_prev[imaging_tests], ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Observed vs modelled distribution — Imaging (", label_prev_nice, ")"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "bottom")


ggsave(file = file.path(output,paste0(run_label,"_fig_hist_imaging_", ref_prev, ".png")), plot = patchwork_hist_imaging_prev,
       width = 16, height = 4 * ceiling(length(imaging_tests) / 2), dpi = 300)

patchwork_hist_endoscopy_prev <- wrap_plots(plots_hist_prev[endoscopy_tests], ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Observed vs modelled distribution — Endoscopy (", label_prev_nice, ")"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "bottom")

ggsave(file = file.path(output,paste0(run_label,"_fig_hist_endoscopy_", ref_prev, ".png")), plot = patchwork_hist_endoscopy_prev,
       width = 16, height = 4 * ceiling(length(endoscopy_tests) / 2), dpi = 300)


## Patchworks — LATEST

patchwork_hist_imaging_latest <- wrap_plots(plots_hist_latest[imaging_tests], ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Observed vs modelled distribution — Imaging (", label_latest_nice, ")"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "bottom")

ggsave(file = file.path(output,paste0(run_label,"_fig_hist_imaging_", ref_latest, ".png")), plot = patchwork_hist_imaging_latest,
       width = 16, height = 4 * ceiling(length(imaging_tests) / 2), dpi = 300)

patchwork_hist_endoscopy_latest <- wrap_plots(plots_hist_latest[endoscopy_tests], ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Observed vs modelled distribution — Endoscopy (", label_latest_nice, ")"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "bottom")

ggsave(file = file.path(output,paste0(run_label,"_fig_hist_endoscopy_", ref_latest, ".png")), plot = patchwork_hist_endoscopy_latest,
       width = 16, height = 4 * ceiling(length(endoscopy_tests) / 2), dpi = 300)

patchwork_hist_imaging_prev
patchwork_hist_endoscopy_prev
patchwork_hist_imaging_latest
patchwork_hist_endoscopy_latest


## #####################################################################
## Build predicted vs observed plots by sub-ICB rank

plots_p <- list()

for (test in all_tests) {
  
  newdf <- all_newdf_df |> filter(TestName == test)
  
  ## ------------------------------------------------------------------
  ## Plot p: predicted + observed by sub-ICB rank
  plots_p[[test]] <- ggplot(
    data = newdf,
    aes(
      y      = Predicted_prob,
      ymin   = Predicted_lcl,
      ymax   = Predicted_ucl,
      x      = ploty,
      colour = Timeperiod,
      fill   = Timeperiod
    )
  ) +
    geom_point(shape = 15) +
    geom_errorbar(width = 0) +
    geom_point(shape = 4, aes(y = Observed_prob, x = ploty - 0.2)) +
    theme_minimal() +
    scale_x_continuous(name = "", breaks = NULL, minor_breaks = NULL, trans = "reverse") +
    scale_y_continuous(name = paste("% waiting 6+ weeks")) +
    labs(title = test)
  
}


## #####################################################################

## Patchwork of previous vs current month predictions plots

## Patchwork of plot p — imaging
patchwork_p_imaging <- wrap_plots(plots_p[imaging_tests], ncol = 2) +
  plot_layout(guides = "collect") +  
  plot_annotation(
    title = "Predicted and observed % waiting 6+ weeks by local area — Imaging",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file = file.path(output,paste0(run_label,"_fig_predictions_imaging.png")), plot = patchwork_p_imaging,
       width = 16, height = 4 * ceiling(length(imaging_tests) / 2), dpi = 300)


## Patchwork of plot p — endoscopy
patchwork_p_endoscopy <- wrap_plots(plots_p[endoscopy_tests], ncol = 2) +
  plot_layout(guides = "collect") +  
  plot_annotation(
    title = "Predicted and observed % waiting 6+ weeks by local area — Endoscopy",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file = file.path(output,paste0(run_label,"_fig_predictions_endoscopy.png")), plot = patchwork_p_endoscopy,
       width = 16, height = 4 * ceiling(length(endoscopy_tests) / 2), dpi = 300)


patchwork_p_imaging
patchwork_p_endoscopy


## #####################################################################

## #####################################################################
## Build sub-ICB appendix table (predicted vs observed + p-values)

appendix_subicb <- all_compare2_df |>
  
  # Keep only the two comparison months
  filter(Timeperiod %in% c(ref_prev, ref_latest)) |>
  
  mutate(
    Timeperiod = factor(Timeperiod, levels = c(ref_prev, ref_latest))
  ) |>
  
  select(
    TestName,
    NHSCode_PostMerge,
    Timeperiod,
    
    Predicted_prob,
    Observed_prob,
    
    Predicted_lcl,
    Predicted_ucl,
    Observed_lcl,
    Observed_ucl,
    
    Observed_Count,
    Observed_Denominator,
    
    Predicted_diff,
    Predicted_diff_lcl,
    Predicted_diff_ucl,
    p.value
  ) |>
  
  # Pivot to wide so prev vs latest are side-by-side
  pivot_wider(
    id_cols = c(TestName, NHSCode_PostMerge, p.value,
                Predicted_diff, Predicted_diff_lcl, Predicted_diff_ucl),
    names_from = Timeperiod,
    values_from = c(
      Predicted_prob,
      Observed_prob,
      Predicted_lcl,
      Predicted_ucl,
      Observed_lcl,
      Observed_ucl,
      Observed_Count,
      Observed_Denominator
    ),
    names_glue = "{.value} ({Timeperiod})"
  ) |> 
  
  rename(
    test = TestName
  ) |>
  
  arrange(test, desc(paste0("Predicted_prob (", ref_prev,")")))

# Add sub ICB names
appendix_subicb <- appendix_subicb |>
  left_join(subicb, by = "NHSCode_PostMerge")

## Clean column ordering
appendix_subicb <- appendix_subicb |>
  select(
    test,
    NHSCode_PostMerge,
    AreaName_PostMerge,
    
    # Previous year first, then 2025 before 2026
    (contains("Observed") & contains("2025")),
    (contains("Observed") & contains("2026")),
    (contains("Predicted") & contains("2025")),
    (contains("Predicted") & contains("2026")),
    
    # Differences + p-value at end
    Predicted_diff,
    Predicted_diff_lcl,
    Predicted_diff_ucl,
    p.value
  )

## Round appendix values (2dp except counts/denominators)

appendix_subicb <- appendix_subicb |>
  mutate(
    across(
      where(is.numeric) & 
        !contains("Count") & 
        !contains("Denominator") &
        !contains("value"),
      ~ round(.x, 2)
    ),
    across(contains("value"),
           ~ round(.x,3)),
    across(
      contains("Count") | contains("Denominator"),
      ~ as.integer(round(.x))
    )
  )


## #####################################################################
## Export results to excel workbook
## #####################################################################

## #####################################################################
# Helper: write data frame to xlsx with a styled label row above the data

write_list_to_xlsx <- function(dfs, path, labels_dict = var_labels) {
  
  wb <- createWorkbook()
  
  label_style <- createStyle(
    fontColour     = "#000000",
    fgFill         = "#DCE6F1",
    textDecoration = "bold",
    wrapText       = TRUE,
    valign         = "top",
    border         = "Bottom",
    borderColour   = "#4472C4"
  )
  
  for (sheet_name in names(dfs)) {
    
    df <- dfs[[sheet_name]]
    rownames(df) <- NULL
    
    addWorksheet(wb, sheet_name)
    
    labels <- setNames(
      ifelse(names(df) %in% names(labels_dict),
             labels_dict[names(df)],
             names(df)),
      names(df)
    )
    
    labels[is.na(labels)] <- names(df)
    
    label_row <- as.data.frame(as.list(labels), stringsAsFactors = FALSE)
    names(label_row) <- names(df)
    
    writeData(wb, sheet_name, label_row, startRow = 1, colNames = FALSE, rowNames = FALSE)
    writeData(wb, sheet_name, df,        startRow = 2, colNames = FALSE, rowNames = FALSE)
    
    addStyle(wb, sheet_name, label_style,
             rows = 1, cols = seq_along(df), gridExpand = TRUE)
    
    setColWidths(wb, sheet_name, cols = seq_along(df), widths = "auto")
    freezePane(wb, sheet_name, firstActiveRow = 2)
  }
  
  saveWorkbook(wb, path, overwrite = TRUE)
}


## #####################################################################
# Define variable labels
# Global definition that can be applied depending on what variables are there


var_labels <- c(
  test           = "Diagnostic test",
  YearMonth      = "Time period",
  
  # ---- Observed vs modelled ----
  obs_mean_logit     = "Observed mean (log-odds scale)",
  obs_variance_logit = "Observed variance (log-odds scale)",
  obs_sd_logit       = "Observed SD (log-odds scale)",
  obs_pct25          = "Observed 25th percentile (%)",
  obs_pct75          = "Observed 75th percentile (%)",
  
  mod_mean_logit     = "Modelled mean (log-odds scale)",
  mod_variance_logit = "Modelled variance (log-odds scale)",
  mod_sd_logit       = "Modelled SD (log-odds scale)",
  mod_pct25          = "Modelled 25th percentile (%)",
  mod_pct75          = "Modelled 75th percentile (%)",
  
  # ---- or differences ----
  obs_or_diff = "Observed OR difference (75th / 25th percentile)",
  mod_or_diff = "Modelled OR difference (75th / 25th percentile)",
  
  # ---- Model parameters ----
  fixef_timeperiod_logit = "Fixed effect for year (log-odds)",
  fixef_timeperiod_OR    = "Fixed effect for year (odds ratio)",
  fixef_timeperiod_pval  = "P-value for fixed effect",
  
  re_intercept_variance  = "Random effect variance: intercept",
  re_slope_variance      = "Random effect variance: time slope",
  re_intercept_slope_cov = "Covariance between random intercept and time slope",
  re_intercept_slope_cor = "Correlation between random intercept and time slope", 
  
  # ---- Labels for predicted sub icb results -------
  NHSCode_PostMerge = "Local area",
  AreaName_PostMerge = "Local area name",
  
  Observed_Count = "Observed count (waiting ≥6 weeks)",
  Observed_Denominator = "Observed denominator (total waiting)",
  
  Predicted_diff     = "Absolute change (modelled)",
  Predicted_diff_lcl = "Lower CI (change)",
  Predicted_diff_ucl = "Upper CI (change)",
  p.value            = "P-value for time period comparison"
)


# ---- Dynamically add observed & predicted labels by month ----
for (m in comparison_labels) {
  
  m_nice <- format(ymd(paste0(m, "-01")), "%b %Y")
  
  # Observed
  var_labels[paste0("Observed_prob (", m, ")")]        <- paste0("Observed % (waiting ≥6 weeks)", m_nice)
  var_labels[paste0("Observed_lcl (", m, ")")]         <- paste0("Observed lower CI — ", m_nice)
  var_labels[paste0("Observed_ucl (", m, ")")]         <- paste0("Observed upper CI — ", m_nice)
  var_labels[paste0("Observed_Count (", m, ")")]       <- paste0("Observed count (waiting ≥6 weeks) — ", m_nice)
  var_labels[paste0("Observed_Denominator (", m, ")")] <- paste0("Observed denominator (total waiting) — ", m_nice)
  
  # Modelled (Predicted)
  var_labels[paste0("Predicted_prob (", m, ")")] <- paste0("Modelled % (waiting ≥6 weeks) — ", m_nice)
  var_labels[paste0("Predicted_lcl (", m, ")")]  <- paste0("Modelled lower CI — ", m_nice)
  var_labels[paste0("Predicted_ucl (", m, ")")]  <- paste0("Modelled upper CI — ", m_nice)
}

# ---- Dynamically add modelled results labels ----
for (m in comparison_labels) {
  
  m_nice <- format(ymd(paste0(m, "-01")), "%b %Y")
  
  var_labels[paste0("mean_logit (", m, ")")]     <- paste("Mean (log-odds scale) —", m_nice)
  var_labels[paste0("variance_logit (", m, ")")] <- paste("Variance (log-odds scale) —", m_nice)
  var_labels[paste0("sd_logit (", m, ")")]       <- paste("SD (log-odds scale) —", m_nice)
  var_labels[paste0("pct25 (", m, ")")]          <- paste("25th percentile (%) —", m_nice)
  var_labels[paste0("pct75 (", m, ")")]          <- paste("75th percentile (%) —", m_nice)
  
  # SD difference columns — only meaningful for latest month but will exist
  # as NA for prev month after pivot, so label both for completeness
  var_labels[paste0("sd_diff (", m, ")")]          <- paste("Difference in SD (latest minus previous, log-odds scale) —", m_nice)
  var_labels[paste0("sd_diff_lower_95 (", m, ")")] <- paste("SD difference lower bound (95% CI) —", m_nice)
  var_labels[paste0("sd_diff_upper_95 (", m, ")")] <- paste("SD difference upper bound (95% CI) —", m_nice)
  var_labels[paste0("sd_diff_pvalue (", m, ")")]   <- paste("P-value for change in SD —", m_nice)
  
}


## #####################################################################
# Load results_table2 from code file 1 and join model OR and p-value,
# then export as xlsx

load(file.path(output, "main_results_table2.RData"))

# main_results_table2 <- main_results_table2 |>
  # rename(test = "Diagnostic test")

model_params_slim <- model_params |>
  select(test, fixef_timeperiod_OR, fixef_timeperiod_pval)

results_table2_with_model <- main_results_table2 |>
  left_join(model_params_slim, by = "test")


## #####################################################

# Load appendix table 
load(file.path(output,"appendix_table.RData"))


## #####################################################################
# Helper: order rows by test (using test_order from config), then by
# any additional grouping variables supplied in `by`

order_by_test <- function(df, test_col = "test", by = NULL) {
  df[[test_col]] <- factor(df[[test_col]], levels = test_order_tables)
  sort_vars <- c(test_col, by)
  df[do.call(order, lapply(sort_vars, function(v) df[[v]])), ]
}

## #####################################################################
# Export

write_list_to_xlsx(
  list(
    "Observed vs modelled"           = order_by_test(main_results_table),
    "Modelled results"               = order_by_test(modelled_results_table),
    "Model parameters"               = order_by_test(model_params),
    "National results"               = order_by_test(results_table2_with_model),
    "Appendix- Observed vs modelled" = order_by_test(obs_vs_modelled),
    "National results - appendix"    = order_by_test(appendix_table),
    "Appendix - local predictions"   = order_by_test(appendix_subicb, by = "NHSCode_PostMerge")
  ),
  path = file.path(output, paste0(run_label, "_all_results_combined3.xlsx"))
)



## ##########################################################
# Useful code

# View available models
names(all_models)

# View specific model
model_gastroscopy <- all_models[["Gastroscopy"]]
summary(model_gastroscopy)