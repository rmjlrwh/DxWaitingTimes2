# ── Single-test practice script: Gastroscopy ─────────────────────────────────
# Deconstructed from the full loop; bootstrap CIs replace the reversed-model approach

library(tidyverse)
library(lme4)

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

# ── Filter to Gastroscopy, two months ────────────────────────────────────────

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

n_boot <- 500

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
ci_lower <- sd_diff - 1.96 * se_sd_diffs
ci_upper <- sd_diff + 1.96 * se_sd_diffs



# ── Collect results ───────────────────────────────────────────────────────────

results <- tibble(
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

print(results)





