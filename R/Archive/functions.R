#Functions file

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
  test           = "Diagnostic test name",
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
  
  # ---- Fold differences ----
  obs_fold_diff = "Observed fold difference (75th / 25th percentile)",
  mod_fold_diff = "Modelled fold difference (75th / 25th percentile)",
  
  # ---- Model parameters ----
  fixef_timeperiod_logit = "Fixed effect for current month (log-odds)",
  fixef_timeperiod_OR    = "Fixed effect for current month (odds ratio)",
  fixef_timeperiod_pval  = "P-value for fixed effect",
  
  re_intercept_variance  = "Random effect variance: intercept",
  re_slope_variance      = "Random effect variance: time slope",
  re_intercept_slope_cov = "Covariance between random intercept and time slope",
  re_intercept_slope_cor = "Correlation between random intercept and time slope", 
  
  # ---- Labels for predicted sub icb results -------
  NHSCode_PostMerge = "Local area",
  AreaName_PostMerge = "Local area name",
  
  Observed_Count = "Observed count",
  Observed_Denominator = "Observed denominator",
  
  Predicted_diff     = "Absolute change (modelled)",
  Predicted_diff_lcl = "Lower CI (change)",
  Predicted_diff_ucl = "Upper CI (change)",
  p.value            = "P-value for time period comparison"
)


# ---- Dynamically add observed & predicted labels by month ----
for (m in comparison_labels) {
  
  m_nice <- format(ymd(paste0(m, "-01")), "%b %Y")
  
  # Observed
  var_labels[paste0("Observed_prob (", m, ")")]        <- paste0("Observed % — ", m_nice)
  var_labels[paste0("Observed_lcl (", m, ")")]         <- paste0("Observed lower CI — ", m_nice)
  var_labels[paste0("Observed_ucl (", m, ")")]         <- paste0("Observed upper CI — ", m_nice)
  var_labels[paste0("Observed_Count (", m, ")")]       <- paste0("Observed count — ", m_nice)
  var_labels[paste0("Observed_Denominator (", m, ")")] <- paste0("Observed denominator — ", m_nice)
  
  # Modelled (Predicted)
  var_labels[paste0("Predicted_prob (", m, ")")] <- paste0("Modelled % — ", m_nice)
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

