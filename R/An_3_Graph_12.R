# AN 3 - Graphs

rm(list=ls())
gc()

library(tidyverse)
library(openxlsx)
library(sf)
library(patchwork)

source("R/config.R")

# Load data
load(file.path(data_out, "compare2.RData"))

# Load sub icb boundaries
sub_icb_boundaries <- st_read(file.path(data_dir,"Sub_Integrated_Care_Board_Locations_April_2023_EN_BSC_-5867181813749012508 (1)/SICBL_APR_2023_EN_BSC.shp"))


## #####################################################################
# Define comparison months
latest_month    <- floor_date(max(as.Date(paste0(
  str_replace(unique(all_compare2_df$Timeperiod), "_", "-"), "-01")),
  na.rm = TRUE), "month")
prev_year_month <- latest_month %m-% years(1)

comparison_months  <- c(latest_month, prev_year_month)
comparison_labels  <- format(comparison_months, "%Y-%m")

ref_latest <- comparison_labels[1]  # e.g. "2026-01"
ref_prev   <- comparison_labels[2]  # e.g. "2025-01"

# Underscore versions for column name matching
ref_latest_col <- str_replace(ref_latest, "-", "_")
ref_prev_col   <- str_replace(ref_prev,   "-", "_")

# Nice labels for titles
label_prev_nice   <- format(latest_month %m-% years(1), "%B %Y")
label_latest_nice <- format(latest_month,                "%B %Y")

## ########################################################################
## Widen and compute changes
## ########################################################################

all_compare2_df <- all_compare2_df |>
  mutate(Timeperiod = str_replace(Timeperiod, "-", "_")) |>
  pivot_wider(
    id_cols     = c(NHSCode_PostMerge, TestName),
    names_from  = Timeperiod,
    values_from = c(Observed_prob, Observed_lcl, Observed_ucl,
                    Predicted_prob, Predicted_lcl, Predicted_ucl,
                    Observed_Count, Observed_Denominator,
                    Predicted_diff, Predicted_diff_lcl, Predicted_diff_ucl,
                    p.value)
  )

all_compare2_df <- all_compare2_df |>
  select(-paste0("p.value_", ref_prev_col)) |>
  mutate(Predicted_diff_rel = (
    (.data[[paste0("Predicted_prob_", ref_latest_col)]] - .data[[paste0("Predicted_prob_", ref_prev_col)]]) /
      .data[[paste0("Predicted_prob_", ref_prev_col)]]
  ))

all_compare2_df <- all_compare2_df |>
  mutate(across(
    contains("prob") | contains("lcl") | contains("ucl") | contains("diff"),
    ~ round(.x * 100, 2)))

## ########################################################################
## Change descriptions
## ########################################################################

pval_col <- paste0("p.value_", ref_latest_col)

all_compare2_df <- all_compare2_df |>
  mutate(Predicted_Change_Description =
           case_when(
             is.na(.data[[pval_col]])                                                              ~ "Not calculated",
             .data[[pval_col]] >= 0.05                                                             ~ "No significant difference",
             .data[[pval_col]] < 0.05 & Predicted_diff_rel < -50                                  ~ "Large decrease",
             .data[[pval_col]] < 0.05 & Predicted_diff_rel <  0 & Predicted_diff_rel >= -50       ~ "Small decrease",
             .data[[pval_col]] < 0.05 & Predicted_diff_rel >  0 & Predicted_diff_rel <=  300      ~ "Small increase",
             .data[[pval_col]] < 0.05 & Predicted_diff_rel >  300                                 ~ "Large increase"
           )
  )

## ########################################################################
## Summaries
## ########################################################################

area_summary_count_predicted <- all_compare2_df |>
  group_by(TestName, Predicted_Change_Description) |>
  count() |>
  mutate(Percentage = round(n / 106 * 100, 2))

## ########################################################################
## Prepare map data
## ########################################################################

sub_icb_boundaries <- sub_icb_boundaries |>
  mutate(
    NHSCode_PostMerge = str_trim(coalesce(
      str_match(SICBL23NM, "ICB - (.*)")[, 2],
      str_match(SICBL23NM, "ICB – (.*)")[, 2]
    ))
  )

map_data <- sub_icb_boundaries |>
  left_join(all_compare2_df, by = "NHSCode_PostMerge")

summary(is.na(map_data$Predicted_diff_rel))

# Relative change categories
rel_change_breaks <- c(-Inf, -75, -50, 0, 50, 100, Inf)
rel_change_labels <- c(
  "≥75% reduction",
  "50–75% reduction",
  "0–50% reduction",
  "0–50% increase",
  "50–100% increase",
  "≥100% increase"
)

rel_change_levels <- c(rel_change_labels, "No significant difference")

rel_change_colours <- c(
  "≥75% reduction"            = "#2C7FB8",
  "50–75% reduction"          = "#73B3D8",
  "0–50% reduction"           = "#BDD7E7",
  "No significant difference" = "white",
  "0–50% increase"            = "#FCBBA1",
  "50–100% increase"          = "#FC9272",
  "≥100% increase"            = "#CB181D")

prob_breaks  <- c(0, 10, 20, 30, 40, 50, Inf)
prob_labels  <- c("< 10%", "10% to 19%", "20% to 29%", "30% to 39%", "40% to 49%", "50% +")
prob_colours <- c(
  "< 10%"      = "#73B3D8",
  "10% to 19%" = "#BDD7E7",
  "20% to 29%" = "#E7EFF6",
  "30% to 39%" = "#FCBBA1",
  "40% to 49%" = "#FC9272",
  "50% +"      = "#FB6A4A"
)

categorise_rel_change <- function(rel_change, change_desc) {
  cat <- as.character(cut(rel_change, breaks = rel_change_breaks, labels = rel_change_labels, right = FALSE))
  cat[change_desc == "No significant difference"] <- "No significant difference"
  factor(cat, levels = rel_change_levels)
}

prob_prev_col   <- paste0("Predicted_prob_", ref_prev_col)
prob_latest_col <- paste0("Predicted_prob_", ref_latest_col)

map_data <- map_data |>
  mutate(
    Predicted_diff_rel_cat      = categorise_rel_change(Predicted_diff_rel, Predicted_Change_Description),
    Predicted_prob_prev_cat     = cut(.data[[prob_prev_col]],   breaks = prob_breaks, labels = prob_labels, right = FALSE),
    Predicted_prob_latest_cat   = cut(.data[[prob_latest_col]], breaks = prob_breaks, labels = prob_labels, right = FALSE)
  )

## ########################################################################
## Map helper functions
## ########################################################################

map_theme <- theme_minimal() +
  theme(
    panel.grid   = element_blank(),
    axis.text    = element_blank(),
    axis.title   = element_blank(),
    plot.title   = element_text(size = 8),
    legend.text  = element_text(size = 7),
    legend.title = element_text(size = 8)
  )

make_prob_map <- function(data, fill_col, title) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "black", linewidth = 0.2) +
    scale_fill_manual(
      values   = prob_colours,
      name     = "Proportion (%)",
      na.value = "white",
      drop     = FALSE,
      limits   = prob_labels
    ) +
    map_theme +
    labs(title = title) +
    theme(legend.position = "none")
}

make_rel_map <- function(data, fill_col, title) {
  ggplot(data) +
    geom_sf(aes(fill = .data[[fill_col]]), colour = "black", linewidth = 0.2) +
    scale_fill_manual(
      values   = rel_change_colours,
      name     = "Relative % change",
      na.value = "white",
      drop     = FALSE,
      limits   = rel_change_levels
    ) +
    map_theme +
    labs(title = title) +
    theme(legend.position = "none")
}

## ########################################################################
## Legends
## ########################################################################

get_legend <- function(plot) cowplot::get_legend(plot)

legend_prob <- get_legend(
  ggplot(
    data.frame(category = factor(prob_labels, levels = prob_labels)),
    aes(x = 1, y = 1, fill = category)
  ) +
    geom_point(shape = 22, size = 5, colour = "black") +
    scale_fill_manual(
      values = prob_colours, name = "Proportion (%)", limits = prob_labels, drop = FALSE,
      guide  = guide_legend(ncol = 3, byrow = TRUE,
                            keyheight = unit(4, "mm"), keywidth = unit(6, "mm"))
    ) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9), legend.text = element_text(size = 7))
)

legend_rel_levels <- c(
  "≥75% reduction",
  "50–75% reduction",
  "0–50% reduction",
  "No significant difference",
  "0–50% increase",
  "50–100% increase",
  "≥100% increase")

legend_rel <- get_legend(
  ggplot(
    data.frame(category = factor(legend_rel_levels, levels = legend_rel_levels)),
    aes(x = 1, y = 1, fill = category)
  ) +
    geom_point(shape = 22, size = 5, colour = "black") +
    scale_fill_manual(
      values = rel_change_colours, name = "Relative % change",
      limits = legend_rel_levels, drop = FALSE,
      guide  = guide_legend(ncol = 4, byrow = TRUE,
                            keyheight = unit(4, "mm"), keywidth = unit(6, "mm"))
    ) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(size = 9), legend.text = element_text(size = 7))
)

legend_row <- patchwork::wrap_elements(legend_prob) | patchwork::wrap_elements(legend_rel)

## ########################################################################
## Per-test patchwork maps
## ########################################################################

tests     <- unique(map_data$TestName[!is.na(map_data$TestName)])
test_maps <- list()

for (test in tests) {
  d <- map_data |> filter(TestName == test)
  
  p2 <- make_prob_map(d, "Predicted_prob_prev_cat",   paste(test, "| Predicted", label_prev_nice))
  p4 <- make_prob_map(d, "Predicted_prob_latest_cat", paste(test, "| Predicted", label_latest_nice))
  p6 <- make_rel_map(d,  "Predicted_diff_rel_cat",    paste(test))
  
  test_maps[[test]] <-
    p2 /
    p4 /
    p6 /
    legend_row +
    plot_layout(heights = c(1, 1, 1, 0.18)) +
    plot_annotation(
      title = paste0(
        "Proportion (%) of patients waiting 6+ weeks for a test in ", label_latest_nice, ",\n",
        "and percentage change relative to ", label_prev_nice, ": \n",
        test
      )
    )
  
  ggsave(
    filename = file.path(output, paste0("map_", gsub("[^a-zA-Z0-9]", "_", test), ".png")),
    plot     = test_maps[[test]],
    width    = 16, height = 18
  )
}

## ########################################################################
## Modality summary patchworks
## ########################################################################

prob_map   <- function(test) make_prob_map(map_data |> filter(TestName == test), "Predicted_prob_latest_cat", paste(test, "| %", label_latest_nice))
prob_map_prev <- function(test) make_prob_map(map_data |> filter(TestName == test), "Predicted_prob_prev_cat",   paste(test, "| %", label_prev_nice))
rel_map    <- function(test) make_rel_map( map_data |> filter(TestName == test), "Predicted_diff_rel_cat",    paste(test, ""))

## ---- Endoscopy ----------------------------------------------------------

endo_prob_plot_prev <-
  (prob_map_prev("Gastroscopy")            | prob_map_prev("Colonoscopy")) /
  (prob_map_prev("Flexible Sigmoidoscopy") | prob_map_prev("Cystoscopy")) /
  patchwork::wrap_elements(legend_prob) +
  plot_layout(heights = c(1, 1, 0.18)) +
  plot_annotation(
    title = paste("Endoscopy: Proportion (%) of patients waiting 6+ weeks,", label_prev_nice),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

endo_prob_plot <-
  (prob_map("Gastroscopy")            | prob_map("Colonoscopy")) /
  (prob_map("Flexible Sigmoidoscopy") | prob_map("Cystoscopy")) /
  patchwork::wrap_elements(legend_prob) +
  plot_layout(heights = c(1, 1, 0.18)) +
  plot_annotation(
    title = paste("Endoscopy: Proportion (%) of patients waiting 6+ weeks,", label_latest_nice),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

endo_rel_plot <-
  (rel_map("Gastroscopy")            | rel_map("Colonoscopy")) /
  (rel_map("Flexible Sigmoidoscopy") | rel_map("Cystoscopy")) /
  patchwork::wrap_elements(legend_rel) +
  plot_layout(heights = c(1, 1, 0.18)) +
  plot_annotation(
    title = paste0("Endoscopy: Relative change (%) in proportion waiting 6+ weeks, ", label_prev_nice, " to ", label_latest_nice),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(output, "map_summary_endoscopy_plot_prev.png"), plot = endo_prob_plot_prev, width = 16, height = 10)
ggsave(file.path(output, "map_summary_endoscopy_prob.png"),   plot = endo_prob_plot,   width = 16, height = 10)
ggsave(file.path(output, "map_summary_endoscopy_rel.png"),    plot = endo_rel_plot,    width = 16, height = 10)

## ---- Imaging ------------------------------------------------------------

img_prob_plot_prev <-
  (prob_map_prev("MRI") | prob_map_prev("CT")) /
  (prob_map_prev("Non Obstetric Ultrasound") | prob_map_prev("Echocardiography")) /
  patchwork::wrap_elements(legend_prob) +
  plot_layout(heights = c(1, 1, 0.18)) +
  plot_annotation(
    title = paste("Imaging: Proportion (%) of patients waiting 6+ weeks,", label_prev_nice),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

img_prob_plot <-
  (prob_map("MRI") | prob_map("CT")) /
  (prob_map("Non Obstetric Ultrasound") | prob_map("Echocardiography")) /
  patchwork::wrap_elements(legend_prob) +
  plot_layout(heights = c(1, 1, 0.18)) +
  plot_annotation(
    title = paste("Imaging: Proportion (%) of patients waiting 6+ weeks,", label_latest_nice),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

img_rel_plot <-
  (rel_map("MRI") | rel_map("CT")) /
  (rel_map("Non Obstetric Ultrasound") | rel_map("Echocardiography")) /
  patchwork::wrap_elements(legend_rel) +
  plot_layout(heights = c(1, 1, 0.18)) +
  plot_annotation(
    title = paste0("Imaging: Relative change (%) in proportion waiting 6+ weeks, ", label_prev_nice, " to ", label_latest_nice),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
  )

ggsave(file.path(output, "map_summary_imaging_plot_prev.png"), plot = img_prob_plot_prev, width = 16, height = 10)
ggsave(file.path(output, "map_summary_imaging_prob.png"),   plot = img_prob_plot,   width = 16, height = 10)
ggsave(file.path(output, "map_summary_imaging_rel.png"),    plot = img_rel_plot,    width = 16, height = 10)

# To draw interactively:
endo_prob_plot_prev
endo_prob_plot
endo_rel_plot
img_prob_plot_prev
img_prob_plot
img_rel_plot



## ########################################################################
## Change summary table (automated months)
## ########################################################################

change_summary_table <- all_compare2_df |>
  filter(!is.na(Predicted_Change_Description)) |>
  mutate(Predicted_diff_rel_cat = categorise_rel_change(Predicted_diff_rel, Predicted_Change_Description)) |>
  group_by(TestName, Predicted_diff_rel_cat) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(Percentage = round(n / 106 * 100, 2)) |>
  pivot_wider(
    names_from  = Predicted_diff_rel_cat,
    values_from = c(n, Percentage),
    values_fill = 0,
    names_glue  = "{Predicted_diff_rel_cat}_{.value}"
  ) |>
  select(
    TestName,
    any_of(c(
      "≥75% reduction_n",            "≥75% reduction_Percentage",
      "50–75% reduction_n",          "50–75% reduction_Percentage",
      "0–50% reduction_n",           "0–50% reduction_Percentage",
      "No significant difference_n", "No significant difference_Percentage",
      "0–50% increase_n",            "0–50% increase_Percentage",
      "50–100% increase_n",          "50–100% increase_Percentage",
      "≥100% increase_n",            "≥100% increase_Percentage"
    ))
  ) |>
  rename(Test = TestName)



## ########################################################################
## Add sheet to existing workbook with matched styling
## ########################################################################

label_style <- createStyle(
  fontColour     = "#000000",
  fgFill         = "#DCE6F1",
  textDecoration = "bold",
  wrapText       = TRUE,
  valign         = "top",
  border         = "Bottom",
  borderColour   = "#4472C4"
)

wb <- loadWorkbook(file.path(output, "all_results_combined.xlsx"))

sheet_name <- paste0("Sub-ICB change ", ref_latest_col)

addWorksheet(wb, sheet_name)

clean_labels <- names(change_summary_table) |>
  str_replace_all("_n$", " (n)") |>
  str_replace_all("_Percentage$", " (%)")

writeData(wb, sheet_name,
          as.data.frame(as.list(clean_labels)),
          startRow = 1, colNames = FALSE, rowNames = FALSE)
writeData(wb, sheet_name, change_summary_table,
          startRow = 2, colNames = FALSE, rowNames = FALSE)

addStyle(wb, sheet_name, label_style,
         rows = 1, cols = seq_along(change_summary_table), gridExpand = TRUE)

setColWidths(wb, sheet_name, cols = seq_along(change_summary_table), widths = "auto")
freezePane(wb, sheet_name, firstActiveRow = 2)

saveWorkbook(wb, file.path(output, "all_results_combined.xlsx"), overwrite = TRUE)