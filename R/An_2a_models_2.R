#AN2 models - FILE A: Run models, predictions, comparisons, save outputs

rm(list=ls())
gc()
# Load packages
library(tidyverse)
library(binom)
library(lme4)
library(patchwork)
library(marginaleffects)
library(openxlsx)

## #####################################################################
# Load file paths
source("R/config.R")

# Define comparison months- Current Month vs Same month previous year
# This is done by referencing the params files
# params = main analsyis, params_sensitivity for the sensitivity analysis

#source("R/params.R")  # swap to params_sensitivity.R for the sensitivity run
 source("R/params_sensitivity.R")  # swap to params_sensitivity.R for the sensitivity run

## #####################################################################

# Load data
load(file.path(data_out, "icb_combined_mergers_alldx.RData"))

# Load sub icb name file
load(file.path(data_out, "subicb.RData"))


## #################################################################

# # Filter table for test runs just to MRI
# icb_combined_mergers_alldx <- icb_combined_mergers_alldx |>
#   filter(TestName == "MRI")

## ####################################################################

# Define comparison months-  Month x vs Same month previous year
comparison_months <- c(latest_month, prev_year_month)
comparison_labels <- format(comparison_months, "%Y-%m")


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
    mutate(YearMonth = factor(YearMonth, levels = c(ref_prev, ref_latest))) |>
    select(TestName, NHSCode_PostMerge, YearMonth, Count, Denominator)  # Remove unnecessary vars to run faster
  
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
  #   variance in previous month (baseline spread across Sub-ICBs)
  # + variance of the random slope (spread of Sub-ICB deviations from the national trend)
  # + 2 x covariance between random intercept and slope (whether starting position predicts divergence)
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

saveRDS(all_models, file = file.path(data_out, paste0(run_label, "_all_models_auto.RData")))

## #####################################################################
# Bind results

modelled_results <- bind_rows(modelled_results_list)
obs_vs_modelled  <- bind_rows(obs_vs_mod_list)
model_params     <- bind_rows(model_params_list)


## ###########################################
## Predict for all tests
## ###########################################

 all_tests <- c(imaging_tests, endoscopy_tests)

all_newdf    <- list()
all_compare2 <- list()

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
    newdata = dat,
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

save(all_compare2_df, file = file.path(data_out, paste0(run_label,"_compare2.RData")))




## #####################################################################
## Save all modelling outputs needed for file 2

save(
  modelled_results,
  obs_vs_modelled,
  model_params,
  all_newdf_df,
  all_compare2_df,
  comparison_labels,
  ref_prev,
  ref_latest,
  label_prev,
  label_latest,
  label_prev_nice,
  label_latest_nice,
  file = file.path(data_out, paste0(run_label,"_AN2_model_outputs.RData"))
)



