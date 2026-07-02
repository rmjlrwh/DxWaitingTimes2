# AN 0 CCG mergers

rm(list=ls())
gc()
# Load packages
library(tidyverse)
library(readxl)
library(stringr)
library(writexl)
library(openxlsx)

source("R/config.R")

# Load  data from data management github. Available at: https://github.com/HOPE-UoM/ACED/tree/main/ICBData 
load(file.path(data_dir, "Metadata_Indicators_ICB.RData"))
load(file.path(data_dir, "ICB_combined.RData"))

# Add metadata for geography
load(file.path(data_dir, "ICB_metadata_ICB.RData"))
load(file.path(data_dir, "ICB_metadata_SUBICB.RData"))
load(file.path(data_dir, "ICB_metadata_CCG.RData"))


## #########################################################################################


# Import merger data succarc file. Available at: https://github.com/HOPE-UoM/ACED/tree/main/ICBData

# Read the sheet named "CCG Mergers"
ccg_mergers <- read_csv(file.path(data_dir, "succarc.csv"),
                        col_names = c("NHSCode", "NHSCode_PostMerge", "Succ_Reason", "Succ_EffectiveDate", "SuccIndicator"))

save(ccg_mergers, file = file.path(data_out,"ccg_mergers.RData"))


## #########################################################################################

# Filter for both 6+ and 13+ week waiting times
icb_combined <- icb_combined %>%
  filter((str_detect(IndicatorName, "Number waiting 6") | str_detect(IndicatorName, "Number waiting 13"))
         & (
           str_detect(IndicatorName, " ct ")
           | str_detect(IndicatorName, " colonoscopy ")
           | str_detect(IndicatorName, " flexi")
           | str_detect(IndicatorName, " gastroscopy ")
           | str_detect(IndicatorName, " mri ")
           | str_detect(IndicatorName, " non obstetric ultrasound ")
           | str_detect(IndicatorName, " cystoscopy ")
           | str_detect(IndicatorName, " echocardiography ")
         )
  ) |>
  filter(AreaType == "ICB sub-locations" | AreaType == "CCGs") |>
  filter(TimeperiodSortable >= 20190000) |>
  mutate(
    Date = ymd(TimeperiodSortable),  # Convert YYYYMMDD to Date
  )



## #########################################################################################

# Merge in new CCG codes
icb_combined <- icb_combined |>
  left_join(ccg_mergers %>% select(`NHSCode`,`NHSCode_PostMerge`), 
            by = "NHSCode") |>
  mutate(NHSCode_PostMerge = case_when(
    is.na(NHSCode_PostMerge) ~ NHSCode,
    TRUE ~ NHSCode_PostMerge)
  ) 

## #########################################################################################

# Add consistent Sub-ICB names first, then fill in with discontinued CCG names

# Rename vars for merge
subicb <- subicb |>
  rename(
    NHSCode_PostMerge = NHSCode
    , AreaCode_PostMerge = AreaCode
    , AreaName_PostMerge = AreaName
    , ParentName_PostMerge = ParentName
    , ParentCode_PostMerge = ParentCode
    , ParentNHSCode_PostMerge = ParentNHSCode
  )

ccg <- ccg |>
  rename(
    NHSCode_PostMerge = NHSCode
    , AreaName_PostMerge = CCGName
  )


# Join sub ICBs data 
# Check if icb_combined contains any rows where AreaType is "subicbs"
if (any(icb_combined$AreaType == "ICB sub-locations" | icb_combined$AreaType == "CCGs")) {
  icb_combined <- icb_combined |>
    left_join(subicb %>% select(`NHSCode_PostMerge`, `AreaCode_PostMerge`, `AreaName_PostMerge`, `ParentName_PostMerge`, `ParentCode_PostMerge`, `ParentNHSCode_PostMerge`), 
              by = "NHSCode_PostMerge")
}


# Join CCGs data - To fill in any area names for CCGs that only exist in CCG lookup (pre sub ICB lookup)
# Can't add other information for CCGs that were merged - only have their name in the ccg lookup
if (any(icb_combined$AreaType == "CCGs" & is.na(icb_combined$AreaName_PostMerge), na.rm = TRUE)) {
  icb_combined <- icb_combined |>
    left_join(ccg %>% select(`NHSCode_PostMerge`, `AreaName_PostMerge`), 
              by = "NHSCode_PostMerge", suffix = c("", "_new")) |>
    mutate(
      AreaName_PostMerge = coalesce(AreaName_PostMerge, AreaName_PostMerge_new)
    ) |>
    select(-ends_with("_new"))
}


## #########################################################################################

# List of CCGs that didn't merge 

# Fine - only non match was NA (England)
missing_NHSCode <- icb_combined %>%
  filter(is.na(NHSCode_PostMerge) & (AreaType == "ICB sub-locations" | AreaType == "CCGs")) %>%
  distinct(NHSCode) %>%
  arrange(NHSCode)



## #########################################################################################

# Collapse by new CCG

icb_combined_mergers <- icb_combined %>%
  group_by(AreaType,NHSCode_PostMerge, AreaName_PostMerge, TimeperiodSortable, Timeperiodrange, SourceID, Year, Date, IndicatorID, IndicatorName,Sex, Age, Timeperiod) %>%
  summarise(
    Count = sum(Count, na.rm = TRUE),
    Denominator = sum(Denominator, na.rm = TRUE),
    .groups = "drop"
  )

# Remove sub iCBs from analysis that underwent changes in April 2026
# 4 sub ICBs should be excluded

# Sub ICB shut - Frimley D4U1Y

# Some LSOAs from Frimley moved into new "Thames Valley " - (code appears to be U2G6B),
# which should be excluded

# Sub ICBs had LSOAs moved into them from Frimley so are bigger now
# NHS Surrey and Sussex 92A
# NHS Hampshire and Isle of Wight D9Y0V

# Also drop 13Q - "national commissioning hub" just exists in 2021 May data

# Define the list of strings to exclude
exclude_list <- c("U2G6B", "92A", "D9Y0V", "D4U1Y", "13Q")

# Filter rows where a column does NOT match any string in the list
icb_combined_mergers <- icb_combined_mergers %>%
  filter(!NHSCode_PostMerge %in% exclude_list)

# Should be 100 sub ICBs in total?
save(icb_combined_mergers, file = file.path(data_out,"icb_combined_mergers.RData"))

## #####################################################################
## Save list of unique sub-ICB codes and names 

subicb2 <- icb_combined_mergers %>%
  distinct(NHSCode_PostMerge, AreaName_PostMerge)

subicb <- subicb |>
  select(NHSCode_PostMerge, AreaName_PostMerge) 

# Find rows in subicb2 that are NOT in subicb = should be 0
extras_in_subicb2 <- anti_join(subicb2, subicb, by = c("NHSCode_PostMerge", "AreaName_PostMerge"))

# Save for later use
save(subicb, file = file.path(data_out,"subicb.RData"))


