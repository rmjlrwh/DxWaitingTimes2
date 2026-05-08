# ── Single-test practice script: Gastroscopy ─────────────────────────────────
# Deconstructed from the full loop; bootstrap CIs replace the reversed-model approach

library(tidyverse)
library(lme4)
library(boot)
library(openxlsx)

source("R/config.R")

load(file.path(data_out, "icb_combined_mergers_alldx.RData"))
load(file.path(data_out, "subicb.RData"))


## TO DO: Remove columns not needed



# ── Time period setup ─────────────────────────────────────────────────────────

latest_month    <- floor_date(max(icb_combined_mergers_alldx$Date, na.rm = TRUE), "month")
prev_year_month <- latest_month %m-% years(1)
comparison_months <- c(latest_month, prev_year_month)
comparison_labels <- format(comparison_months, "%Y-%m")

ref_prev   <- comparison_labels[2]
ref_latest <- comparison_labels[1]

label_prev        <- format(ymd(paste0(ref_prev,   "-01")), "%b_%Y")
label_latest      <- format(ymd(paste0(ref_latest, "-01")), "%b_%Y")
label_prev_nice   <- format(ymd(paste0(ref_prev,   "-01")), "%b %Y")
label_latest_nice <- format(ymd(paste0(ref_latest, "-01")), "%b %Y")

# ── Filter to one test specified below, two months ────────────────────────────────────────

data_2month <- icb_combined_mergers_alldx |>
  mutate(YearMonth = format(Date, "%Y-%m")) |>
  filter(YearMonth %in% comparison_labels) |>
  mutate(
    NHSCode_PostMerge = factor(NHSCode_PostMerge),
    Timeperiod        = factor(Timeperiod),
    YearMonth         = factor(YearMonth)
  )

test <- "MRI"

dat <- data_2month |>
  filter(TestName == test) |>
  mutate(YearMonth = factor(YearMonth, levels = c(ref_prev, ref_latest)))



# ── Fit model ───────────────────────────────────────────────────────

model <- glmer(
  cbind(Count, Denominator - Count) ~ YearMonth + (YearMonth | NHSCode_PostMerge),
  family = binomial,
  data   = dat
)

# ── Extract real-data SDs and difference in SDs ───────────────────────────────

vc          <- as.matrix(VarCorr(model)$NHSCode_PostMerge)
fixef_      <- fixef(model)
time_effect <- names(fixef_)[2]

var_prev   <- vc["(Intercept)", "(Intercept)"]
var_latest <- vc["(Intercept)", "(Intercept)"] +
  vc[time_effect,   time_effect] +
  2 * vc["(Intercept)", time_effect]

sd_prev   <- sqrt(var_prev)
sd_latest <- sqrt(var_latest)
sd_diff   <- sd_latest - sd_prev          # create CI for this SD difference



##################################################################################

# ##########################################
## Boostrap CIs for SD difference
# ##########################################

n_boot <- 1000

# Unique sub-ICB IDs present in the data
ids <- unique(dat$NHSCode_PostMerge)
n_ids <- length(ids)

boot_sd_diffs <- numeric(n_boot)

for (b in seq_len(n_boot)) # Repeat for number of boostraps
  
  {
  
  # Sample IDs with replacement. Sampling total number of subicbs in dataset (106) using n_ids
  sampled_ids <- sample(ids, size = n_ids, replace = TRUE)
  
  # Build bootstrapped dataset, using sampled IDs
  # relabelling so duplicated IDs are treated as distinct sub-ICBs
  
  boot_dat <- map_dfr(seq_along(sampled_ids), function(i) # NB gives map_dfr no. rows to iterate over that matches number of sampled_ids (106)
    {
    dat |>
      filter(NHSCode_PostMerge == sampled_ids[i]) |> # Filter the main data to row positions that match those in the randomly drawn sampled_ids
      mutate(NHSCode_PostMerge = factor(paste0("boot_", i))) # Give subICBs a new unique ID if the same one has been drawn twice
  })
  
  # Check each boot ID has both months
  check <- boot_dat |>
    distinct(NHSCode_PostMerge, YearMonth) |>
    count(NHSCode_PostMerge)
  
  if (any(check$n != 2)) {
    warning(paste("Bootstrap sample", b, "has incomplete time pairs"))
  }
  
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


# ── Derive CI and p-value for the difference in SDs ──────────────────────────

# SD of bootstrap differences = standard error of bootstrapped differences between SDs in prev month vs latest
se_sd_diffs <- sd(boot_sd_diffs)

# z_stat = difference in SDs in actual dataset / standard error for the SDs from the bootstrapping
# Converts difference in SDs into standard errors away from zero
# If z_stat = 0, no change
z_stat  <- sd_diff / se_sd_diffs

# Draw normal distribution of z_stat to create p value
# probability of getting a z score that extreme by chance
# 2 sided (hence -abs to make all values positive)
p_value <- 2 * pnorm(-abs(z_stat))

# 95% CI using normal distribution
ci_lower <- sd_diff - qnorm(0.975) * se_sd_diffs
ci_upper <- sd_diff + qnorm(0.975) * se_sd_diffs



# ── Collect results ───────────────────────────────────────────────────────────

bootstrap_results_manual <- tibble(
  test          = test,
  sd_prev       = round(sd_prev,   4),
  sd_latest     = round(sd_latest, 4),
  sd_diff       = round(sd_diff,   4),
  ci_lower_95   = round(ci_lower,  4),
  ci_upper_95   = round(ci_upper,  4),
  se_sd_diffs    = round(se_sd_diffs, 4),
  z_stat        = round(z_stat,    3),
  p_value       = round(p_value,   4),
  n_boot        = n_boot
)

print(bootstrap_results_manual)





## #####################################################################
## Compare to results from Boot package, and different distribution (e.g. norm vs bca)
## #####################################################################

# Prep data
data_bs <- dat |>
  select(
    NHSCode_PostMerge,
    YearMonth,
    Count,
    Denominator
  ) |>
  mutate(NHSCode_PostMerge = as.numeric(NHSCode_PostMerge))

# Bootstrapping... make some one-row-per-id data to simplify the sampling
data_bs <- data_bs |>
  pivot_wider(
    id_cols = NHSCode_PostMerge,
    names_from = YearMonth,
    values_from = c(Count, Denominator)
  )


# Function for modelling --------------------------------------------------

# Function to extract model components
calc_stats <- function(model) {
  
  # Extract model components
  vc          <- as.matrix(VarCorr(model)$NHSCode_PostMerge)
  fixef_      <- fixef(model)
  time_effect <- names(fixef_)[2]
  
  # Calc variances and intercepts
  var_prev   <- vc["(Intercept)", "(Intercept)"]
  var_latest <- vc["(Intercept)", "(Intercept)"] +
    vc[time_effect,   time_effect] +
    2 * vc["(Intercept)", time_effect]
  
  sd_prev   <- sqrt(var_prev)
  sd_latest <- sqrt(var_latest)
  diff_sds   <- sd_latest - sd_prev          # create CI for this SD difference
  
  return(diff_sds)

}


# Function to calculate different SD
diff_sds <- function(data_bs, NHSCode_PostMerge = NHSCode_PostMerge) {
  
  # Ready data frame for modelling
  
  bs_df <- data_bs[NHSCode_PostMerge,] |>
    mutate(NHSCode_PostMerge = row_number()) |>
    pivot_longer(cols = c(
      paste0("Count_",ref_prev),
      paste0("Count_",ref_latest),
      paste0("Denominator_",ref_prev),
      paste0("Denominator_",ref_latest)
      )) |>
    mutate(
      YearMonth = str_extract(name, "\\d{4}-\\d{2}") |> factor(),
      name = str_remove(name, "_\\d{4}-\\d{2}") |> factor()
    ) |>
    pivot_wider(
      id_cols = c(NHSCode_PostMerge, YearMonth),
      names_from = name,
      values_from = value
    ) 
  
  # Fit the model
  bs_model <- glmer(
    cbind(Count, Denominator - Count) ~ YearMonth + (YearMonth | NHSCode_PostMerge),
    family = binomial,
    data   = bs_df
  )
  
  # Extract bootstrap results
  res_bs <- calc_stats(model = bs_model)
  return(res_bs)
  
}

## #####################################################################
# Fit model for Gastroscopy only
bs_results <- boot::boot(data = data_bs, statistic = diff_sds, R = 1000, sim = "ordinary", stype = "i")
bs_results |> boot::boot.ci(type = c("norm","basic", "perc", "bca"))


## Ready for outputs

# Convert boot output to tidy table
bootstrap_results_package <- tibble(
  statistic_original = bs_results$t0,
  mean_boot          = mean(bs_results$t),
  sd_boot            = sd(bs_results$t),
  ci_norm_lower      = boot::boot.ci(bs_results, type = "norm")$normal[2],
  ci_norm_upper      = boot::boot.ci(bs_results, type = "norm")$normal[3],
  ci_perc_lower      = boot::boot.ci(bs_results, type = "perc")$percent[4],
  ci_perc_upper      = boot::boot.ci(bs_results, type = "perc")$percent[5],
  ci_bca_lower       = boot::boot.ci(bs_results, type = "bca")$bca[4],
  ci_bca_upper       = boot::boot.ci(bs_results, type = "bca")$bca[5]
)

## Save to workbook

write_list_to_xlsx(
  list(
    "Bootstrap summary (manual)"   = bootstrap_results_manual,
    "Bootstrap results (boot pkg)" = bootstrap_results_package
    ),
  path = file.path(output, "bootstrap_results.xlsx")
)



