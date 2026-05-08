#AN2 models 

rm(list=ls())
gc()
# Load packages
library(tidyverse)
library(binom)
library(lme4)
library(patchwork)
library(marginaleffects)
library(openxlsx)

source("R/config.R")

# Load data
load(file.path(data_out, "icb_combined_mergers_alldx.RData"))

# Load sub icb name file
load(file.path(data_out, "subicb.RData"))


## #####################################################################

# Define comparison months- Current Month vs Same month previous year
latest_month <- floor_date(max(icb_combined_mergers_alldx$Date, na.rm = TRUE), "month")
prev_year_month <- latest_month %m-% years(1)

comparison_months <- c(latest_month, prev_year_month)
comparison_labels <- format(comparison_months, "%Y-%m")

## #####################################################################

# Run functions file, which depends on definition of comparison months above
# Includes functions run across multiple R code files e.g. labels and excel exports
source("R/functions.R")



## #####################################################################

# Filter table
data_2month <- icb_combined_mergers_alldx |>
  mutate(YearMonth = format(Date, "%Y-%m")) |>
  filter(YearMonth %in% comparison_labels)


## #####################################################################
# # Format factor vars
data_2month <- data_2month %>%
  mutate(
    NHSCode_PostMerge = factor(NHSCode_PostMerge),
    Timeperiod        = factor(Timeperiod),
    YearMonth = factor(YearMonth)
  )


## #####################################################################
# Define test groups
all_tests       <- unique(data_2month$TestName)

# Define reference levels for modelling
ref_prev   <- comparison_labels[2]  # e.g. "2025-01"
ref_latest <- comparison_labels[1]  # e.g. "2026-01"

# Define labels for outputs

# Names for data columns
label_prev   <- format(ymd(paste0(ref_prev, "-01")), "%b_%Y")
label_latest <- format(ymd(paste0(ref_latest, "-01")), "%b_%Y")

# Names for graphs
label_prev_nice   <- format(ymd(paste0(ref_prev, "-01")), "%b %Y")
label_latest_nice <- format(ymd(paste0(ref_latest, "-01")), "%b %Y")


## #####################################################################
# Loop over all tests to fit models and collect results

modelled_results_list  <- list()
obs_vs_mod_list        <- list()
model_params_list      <- list()
all_models             <- list()

n_boot <- 1000

for (test in all_tests) {
  
  dat <- data_2month %>%
    filter(TestName == test) %>%
    mutate(YearMonth = factor(YearMonth, levels = c(ref_prev, ref_latest)))
  
  ## ── Observed ────────────────────────────────────────────────────────
  
  # Format variables
  observed <- dat %>%
    mutate(
      probability = Count / Denominator, 
      Timeperiod = format(ymd(paste0(YearMonth, "-01")), "%b_%Y")
    ) %>%
    select(TestName, NHSCode_PostMerge, Timeperiod, YearMonth, Count, Denominator, probability)
  
  observed_log <- observed |>
    
    #transform indicator to log odds scale
    mutate(probability_log = qlogis(probability)) |> 
    
    # continuity correction to produce approximate values on the log-odds scale
    mutate(probability_log = case_when(
      Count == 0   ~ qlogis(0.5 / Denominator),
      TRUE         ~ probability_log
    )) |>
    mutate(probability_log = case_when(
      probability == 1 ~ qlogis((Denominator - 0.5) / Denominator),
      TRUE             ~ probability_log
    ))
  
  # Extract summary statistics for observed indicator on log-odds scale
  obs_log_summary <- observed_log %>%
    group_by(YearMonth, Timeperiod) %>%
    summarise(
      obs_mean_logit = mean(probability_log, na.rm = TRUE),
      obs_sd_logit   = sd(probability_log,   na.rm = TRUE),
      obs_var_logit  = obs_sd_logit^2,
      .groups = "drop"
    )
  
  # Use summary statistics to extract 25th and 75th percentiles
  # Also convert back to natural scale
  obs_pctiles <- observed_log %>%
    group_by(YearMonth, Timeperiod) %>%
    summarise(
      obs_pct25 = plogis(quantile(probability_log, 0.25, na.rm = TRUE)),
      obs_pct75 = plogis(quantile(probability_log, 0.75, na.rm = TRUE)),
      .groups = "drop"
    )
  
  ## ── Models ──────────────────────────────────────────────────────────
  
  model <- glmer(
    cbind(Count, Denominator - Count) ~ YearMonth + (YearMonth | NHSCode_PostMerge),
    family = binomial,
    data   = dat
  )
  
  all_models[[test]] <- model
  
  ## ── Variance-covariance & fixed effects ─────────────────────────────
  
  vc          <- as.matrix(VarCorr(model)$NHSCode_PostMerge)
  fixef_      <- fixef(model)
  time_effect <- names(fixef_)[2]
  
  # Extract intercept (mean logit) for previous month
  mean_logit_prev   <- fixef_["(Intercept)"]
  
  # Extract intercept (mean logit) for latest month
  mean_logit_latest <- fixef_["(Intercept)"] + fixef_[time_effect]
  
  # Variance for previous month
  var_prev <- vc["(Intercept)", "(Intercept)"]
  
  # Variance for later month =
  # variance in previous month 
  #+ random slope variance (how much variance changes by latest month) 
  # + 2 x random intercept-slope covariance
  var_latest <- vc["(Intercept)", "(Intercept)"] +
    vc[time_effect, time_effect] +
    2 * vc["(Intercept)", time_effect]
  
  # SD for previous month
  sd_prev   <- sqrt(var_prev)
  
  # SD for latest month
  sd_latest <- sqrt(var_latest)
  
  # Difference in SD
  sd_diff   <- sd_latest - sd_prev
  
  
  ## ── Bootstrap CIs for SD difference ─────────────────────────────────
  
  ids   <- unique(dat$NHSCode_PostMerge)
  n_ids <- length(ids)
  
  boot_sd_diffs <- numeric(n_boot)
  
  for (b in seq_len(n_boot)) {
    
    # Sample IDs with replacement. Sampling total number of subicbs in dataset (106) using n_ids
    sampled_ids <- sample(ids, size = n_ids, replace = TRUE)
  
    # Build bootstrapped dataset, using sampled IDs
    # relabelling so duplicated IDs are treated as distinct sub-ICBs
    boot_dat <- map_dfr(seq_along(sampled_ids), function(i) {
      dat |>
        filter(NHSCode_PostMerge == sampled_ids[i]) |>
        mutate(NHSCode_PostMerge = factor(paste0("boot_", i)))
    })
  
    # Refit model on bootstrapped sample
    boot_model <- glmer(
      cbind(Count, Denominator - Count) ~ YearMonth + (YearMonth | NHSCode_PostMerge),
      family = binomial,
      data   = boot_dat
    )
    
    # Extract SDs and their difference from the bootstrapped model
    boot_vc          <- as.matrix(VarCorr(boot_model)$NHSCode_PostMerge)
    boot_time_effect <- names(fixef(boot_model))[2]
    
    boot_var_prev   <- boot_vc["(Intercept)", "(Intercept)"]
    boot_var_latest <- boot_vc["(Intercept)", "(Intercept)"] +
      boot_vc[boot_time_effect, boot_time_effect] +
      2 * boot_vc["(Intercept)", boot_time_effect]
    
    boot_sd_prev   <- sqrt(boot_var_prev)
    boot_sd_latest <- sqrt(boot_var_latest)
    
    # Creates list of differences in SDs from each boostrapping round
    # nb [b] references the var name in the forloop header
    boot_sd_diffs[b] <- boot_sd_latest - boot_sd_prev
    
  }
  
  # SD of bootstrap differences = standard error of bootstrapped differences between SDs in prev month vs latest
  se_sd_diff <- sd(boot_sd_diffs)
  
  # 95% CI using normal distribution
  ci_lower <- sd_diff - qnorm(0.975) * se_sd_diff
  ci_upper <- sd_diff + qnorm(0.975) * se_sd_diff
  
  # z_stat = difference in SDs in actual dataset / standard error for the SDs from the bootstrapping
  # Converts difference in SDs into standard errors away from zero
  # If z_stat = 0, no change
  z_stat  <- sd_diff / se_sd_diff
  
  # Draw normal distribution of z_stat to create p value
  # probability of getting a z score that extreme by chance
  # 2 sided (hence -abs to make all values positive)
  p_value <- 2 * pnorm(-abs(z_stat))
  

  ## ── Percentile range results ─────────────────────────────────────────
  
  range_results <- tibble(
    pctile     = c(0.25, 0.75),
    
    logit_prev   = (qnorm(pctile) * sd_prev)   + mean_logit_prev,
    logit_latest = (qnorm(pctile) * sd_latest)  + mean_logit_latest,
    
    perc_prev    = plogis(logit_prev),
    perc_latest  = plogis(logit_latest)
  )
  
  p25 <- range_results |> filter(pctile == 0.25)
  p75 <- range_results |> filter(pctile == 0.75)
  
  ## ── Model summary parameters ─────────────────────────────────────────
  model_params_list[[test]] <- tibble(
    test                   = test,
    fixef_timeperiod_logit = round(fixef_[time_effect],                                                                3),
    fixef_timeperiod_OR    = round(exp(fixef_[time_effect]),                                                           3),
    fixef_timeperiod_pval  = round(summary(model)$coefficients[time_effect, "Pr(>|z|)"],                              3),
    re_intercept_variance  = round(vc["(Intercept)", "(Intercept)"],                                                   3),
    re_slope_variance      = round(vc[time_effect, time_effect],                                                       3),
    re_intercept_slope_cov = round(vc["(Intercept)", time_effect],                                                     3),
    re_intercept_slope_cor = round(vc["(Intercept)", time_effect] /
                                     sqrt(vc["(Intercept)", "(Intercept)"] * vc[time_effect, time_effect]),            3)
  )
  
  ## ── Modelled results rows ────────────────────────────────────────────
  modelled_results_list[[test]] <- tibble(
    test        = test,
    YearMonth   = c(ref_prev,      ref_latest),
    mean_logit  = round(c(mean_logit_prev,  mean_logit_latest), 2),
    variance_logit = round(c(var_prev,      var_latest),        2),
    sd_logit    = round(c(sd_prev,          sd_latest),         2),
    sd_diff     = round(c(sd_diff,         sd_diff),           4),
    sd_diff_lower_95 = round(c(ci_lower,    ci_lower),          4),
    sd_diff_upper_95 = round(c(ci_upper,    ci_upper),          4),
    sd_diff_pvalue   = round(c(p_value,    p_value),           4),
    pct25       = round(c(p25$perc_prev,    p25$perc_latest) * 100, 2),
    pct75       = round(c(p75$perc_prev,    p75$perc_latest) * 100, 2),
  )
  
  ## ── Observed vs modelled rows ────────────────────────────────────────
  obs_perc       <- obs_pctiles     |> filter(YearMonth == ref_prev)
  obs_latest     <- obs_pctiles     |> filter(YearMonth == ref_latest)
  obs_sum_prev   <- obs_log_summary |> filter(YearMonth == ref_prev)
  obs_sum_latest <- obs_log_summary |> filter(YearMonth == ref_latest)
  
  obs_vs_mod_list[[test]] <- tibble(
    test        = test,
    YearMonth = c(ref_prev, ref_latest),
    
    obs_mean_logit     = round(c(obs_sum_prev$obs_mean_logit, obs_sum_latest$obs_mean_logit), 2),
    obs_variance_logit = round(c(obs_sum_prev$obs_var_logit,  obs_sum_latest$obs_var_logit),  2),
    obs_sd_logit       = round(c(obs_sum_prev$obs_sd_logit,   obs_sum_latest$obs_sd_logit),   2),
    obs_pct25          = round(c(obs_perc$obs_pct25, obs_latest$obs_pct25) * 100, 2),
    obs_pct75          = round(c(obs_perc$obs_pct75, obs_latest$obs_pct75) * 100, 2),
    
    mod_mean_logit     = round(c(mean_logit_prev, mean_logit_latest), 2),
    mod_variance_logit = round(c(var_prev,        var_latest),        2),
    mod_sd_logit       = round(c(sd_prev,         sd_latest),         2),
    mod_pct25          = round(c(p25$perc_prev,   p25$perc_latest)  * 100, 2),
    mod_pct75          = round(c(p75$perc_prev,   p75$perc_latest)  * 100, 2)
  )
  
}

saveRDS(all_models, file = file.path(data_out,"all_models_auto.RData"))


## #####################################################################
# Bind results

modelled_results <- bind_rows(modelled_results_list)
obs_vs_modelled  <- bind_rows(obs_vs_mod_list)
model_params     <- bind_rows(model_params_list)


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

# Fold difference
obs_vs_modelled <- obs_vs_modelled |>
  mutate(
    obs_fold_diff = obs_pct75 / obs_pct25,
    mod_fold_diff = mod_pct75 / mod_pct25
  ) |>
  mutate(
    across(ends_with("fold_diff"), ~ round(.x, 2)),
    across(c(obs_pct25, obs_pct75, mod_pct25, mod_pct75), ~ round(.x, 1))
  )


# Shortned main results table
main_results_table <- obs_vs_modelled |>
  select(
    test,
    YearMonth,
    obs_pct25, obs_pct75, obs_fold_diff,
    mod_pct25, mod_pct75, mod_fold_diff
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

ggsave(file = file.path(output,paste0("fig_hist_imaging_", ref_prev, ".png")), plot = patchwork_hist_imaging_prev,
       width = 16, height = 4 * ceiling(length(imaging_tests) / 2), dpi = 300)

patchwork_hist_endoscopy_prev <- wrap_plots(plots_hist_prev[endoscopy_tests], ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Observed vs modelled distribution — Endoscopy (", label_prev_nice, ")"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "bottom")

ggsave(file = file.path(output,paste0("fig_hist_endoscopy_", ref_prev, ".png")), plot = patchwork_hist_endoscopy_prev,
       width = 16, height = 4 * ceiling(length(endoscopy_tests) / 2), dpi = 300)


## Patchworks — LATEST

patchwork_hist_imaging_latest <- wrap_plots(plots_hist_latest[imaging_tests], ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Observed vs modelled distribution — Imaging (", label_latest_nice, ")"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "bottom")

ggsave(file = file.path(output,paste0("fig_hist_imaging_", ref_latest, ".png")), plot = patchwork_hist_imaging_latest,
       width = 16, height = 4 * ceiling(length(imaging_tests) / 2), dpi = 300)

patchwork_hist_endoscopy_latest <- wrap_plots(plots_hist_latest[endoscopy_tests], ncol = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = paste0("Observed vs modelled distribution — Endoscopy (", label_latest_nice, ")"),
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  ) &
  theme(legend.position = "bottom")

ggsave(file = file.path(output,paste0("fig_hist_endoscopy_", ref_latest, ".png")), plot = patchwork_hist_endoscopy_latest,
       width = 16, height = 4 * ceiling(length(endoscopy_tests) / 2), dpi = 300)

patchwork_hist_imaging_prev
patchwork_hist_endoscopy_prev
patchwork_hist_imaging_latest
patchwork_hist_endoscopy_latest




## ###########################################
## Predict for all tests
## ###########################################

all_tests <- c(imaging_tests, endoscopy_tests)

all_newdf    <- list()
all_compare2 <- list()
plots_p      <- list()

for (test in all_tests) {
  
  model <- all_models[[test]]
  
  ## ------------------------------------------------------------------
  ## Get test-specific data
  dat <- data_2month %>%
    filter(TestName == test) %>%
    mutate(YearMonth = factor(as.character(YearMonth), levels = c(ref_prev, ref_latest))) %>%
    droplevels()
  
  ## ------------------------------------------------------------------
  ## Predictions ONLY on that test's data
  # Using marginaleffects::predictions
  newdf <- predictions(
    model,
    newdata = dat,     #
    re.form = NULL
  ) |>
    mutate(
      Observed_prob = Count / Denominator,
      Observed_lcl  = binom.confint(Count, Denominator, method = "wilson")$lower,
      Observed_ucl  = binom.confint(Count, Denominator, method = "wilson")$upper
    ) |>
    select(
      TestName, NHSCode_PostMerge, Timeperiod, YearMonth, 
      Predicted_prob = estimate,
      Predicted_lcl  = conf.low,
      Predicted_ucl  = conf.high,
      Observed_prob, Observed_lcl, Observed_ucl,
      Observed_Count = Count,
      Observed_Denominator = Denominator
    ) 
  
  ## ------------------------------------------------------------------
  ## Create ranking based ONLY on previous (same month a year ago)
  # This will be used to rank observations in a then vs now comparison graph later
  rank_prev <- newdf |>
    filter(Timeperiod == ref_prev) |>
    arrange(desc(Predicted_prob)) |>
    mutate(ploty = row_number()) |>
    select(NHSCode_PostMerge, ploty)
  
  ## ------------------------------------------------------------------
  ## Join ranking back to ALL data (both years)
  newdf <- newdf |>
    left_join(rank_prev, by = "NHSCode_PostMerge")
  
  
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
  
  ## ------------------------------------------------------------------
  ## Comparisons
  # Compare previous vs latest month using marginaleffects::comparisons
  # Although not explicitly spcified, this generates P Values for the comparison
  compare <- comparisons(
    model,
    variables = list(YearMonth = c(ref_prev, ref_latest)),
    by        = "NHSCode_PostMerge",
    type      = "response"
  ) |>
    rename(
      Predicted_diff     = estimate,
      Predicted_diff_lcl = conf.low,
      Predicted_diff_ucl = conf.high
    )
  

  ## ------------------------------------------------------------------
  ##
  # Join comparisons into actual predictions data frame
  compare2 <- newdf |>
    select(TestName, NHSCode_PostMerge, Timeperiod,
           Predicted_prob, Predicted_lcl, Predicted_ucl,
           Observed_prob, Observed_lcl, Observed_ucl,
           Observed_Count, Observed_Denominator
    ) |>
    left_join(compare, by = "NHSCode_PostMerge")
  
  all_newdf[[test]]    <- newdf
  all_compare2[[test]] <- compare2
}

## #####################################################################
## Combine and save
all_newdf_df    <- bind_rows(all_newdf)
all_compare2_df <- bind_rows(all_compare2)

save(all_compare2_df, file = file.path(data_out,"compare2.RData"))


## #####################################################################

## Patchwork of previous vs current month predictions plots

## Patchwork of plot p — imaging
patchwork_p_imaging <- wrap_plots(plots_p[imaging_tests], ncol = 2) +
  plot_layout(guides = "collect") +  
  plot_annotation(
    title = "Predicted and observed % waiting 6+ weeks by local area — Imaging",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file = file.path(output,"fig_predictions_imaging.png"), plot = patchwork_p_imaging,
       width = 16, height = 4 * ceiling(length(imaging_tests) / 2), dpi = 300)


## Patchwork of plot p — endoscopy
patchwork_p_endoscopy <- wrap_plots(plots_p[endoscopy_tests], ncol = 2) +
  plot_layout(guides = "collect") +  
  plot_annotation(
    title = "Predicted and observed % waiting 6+ weeks by local area — Endoscopy",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

ggsave(file = file.path(output,"fig_predictions_endoscopy.png"), plot = patchwork_p_endoscopy,
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
# Load results_table2 from code file 1 and join model OR and p-value,
# then export as xlsx

load(file.path(output, "main_results_table2.RData"))

main_results_table2 <- main_results_table2 |>
  rename(test = "Diagnostic test name")

model_params_slim <- model_params |>
  select(test, fixef_timeperiod_OR, fixef_timeperiod_pval)

results_table2_with_model <- main_results_table2 |>
  left_join(model_params_slim, by = "test")


## #####################################################

# Load appendix table 
load(file.path(output, "appendix_table.RData"))





## #####################################################################
# Export

write_list_to_xlsx(
  list(
    "Observed vs modelled"                         = main_results_table,
    "Modelled results"                             = modelled_results_table,
    "Model parameters"                             = model_params,
    "National results"                             = results_table2_with_model,
    "Appendix- Observed vs modelled" = obs_vs_modelled,
    "National results - appendix" = appendix_table, 
    "Appendix - local predictions"     = appendix_subicb
  ),
  path = file.path(output, "all_results_combined2.xlsx")
)


## ##########################################################
# Useful code

# View available models
names(all_models)

# View specific model
model_gastroscopy <- all_models[["Gastroscopy"]]
summary(model_gastroscopy)

