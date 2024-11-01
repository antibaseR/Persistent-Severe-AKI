library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
library(ggalluvial)

impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# Connect to the server and download the data
library(DBI)
con = dbConnect(odbc::odbc(), "CCM_EHR", timeout = 10)

kdigo_data = tbl(con, "vw_Hidenic_KDIGO_Scr_RRT_UO_FINAL_New") %>% 
  collect()

# variables that needs to be imputed
imp_var = c('age_at_admission',
            'gender',
            'race',
            'charlson',
            'apache3',
            'TOTAL_INPUTS_BEFORE_AKI',
            'TOTAL_SALINE_BEFORE_AKI',
            'TOTAL_OUTPUTS_BEFORE_AKI',
            'TOTAL_Urine_Output_BEFORE_AKI',
            'CUMULATIVE_BALANCE_BEFORE_AKI',
            'MAX_LACTATE',
            'MAX_CHLORIDE',
            "MAX_SODIUM",
            'MIN_ALBUMIN',
            'MIN_PLATELETS',
            'refCreat',
            'weight',
            'height',
            'BMI',
            'sofa_24',
            'MAX_DBP',
            'MAX_HR',
            'MIN_SpO2',
            'AVG_RR',
            'MAX_TEMP',
            'MIN_HGB',
            'MAX_TOTAL_BILI',
            'MAX_INR',
            'MAX_WBC',
            'MIN_PH',
            'BASELINE_eGFR',
            'MIN_PULSE_PRESS_BEFORE_AKI',
            'MAX_PULSE_PRESS_BEFORE_AKI',
            'MIN_MAP_BEFORE_AKI')

missingness = tibble()
missingness_by_AKI = tibble()

for (i in 1:length(impDat)){
  
  # get each imputed dataset
  imp = impDat[[i]]
  # replace the imputed variables
  #imp = imp[order(c(imp$patientid, imp$encounter_start_date_time_sh)),]
  imp = imp %>%
    arrange(patientid, patientvisitid)
  data = rawData
  data = data %>%
    arrange(patientid, patientvisitid)
  data = data %>%
    select(-all_of(imp_var))
  imp_sub = imp %>%
    select(patientid, patientvisitid, all_of(imp_var))
  data = merge(data, imp_sub, by = c("patientid", "patientvisitid"))
  
  ## exclusion
  #The exclusions we need to have for the persistent AKI primary analysis are:
  #1. Kidney transplant
  #2. End stage Renal disease (by ICD-9/10 codes or ICD9/10 codes equivalent to CKD stage 4 or 5)
  #3. eGFR <15 at admission
  #4. Creatinine ≥ 4 mg/dL at admission
  #5. ECMO
  #6. Only use the first encounter
  data = data %>%
    filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
  data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
  #subset to rows that are not duplicates of a previous encounter for that patient
  data = data[!duplicated(data$patientid),]
  
  data = data %>%
    select(patientid, patientvisitid, AKI_Category_subgroup)
  kdigo_data_sub = kdigo_data %>%
    select(patientid, patientvisitid, kdigo, event_Date_sh) %>%
    filter(patientvisitid %in% data$patientvisitid)
  
  # Assuming df has columns 'patientvisitid' and 'event_Date_sh'
  relative_day_df = kdigo_data_sub %>%
    # Convert 'event_Date_sh' to Date if it is not already in Date format
    mutate(event_Date = as.Date(event_Date_sh)) %>%
    # Group by 'patientvisitid'
    group_by(patientvisitid) %>%
    # Arrange by 'event_Date_sh' to ensure order
    arrange(event_Date) %>%
    # Calculate the relative day for each event_Date_sh
    mutate(relative_day = as.numeric(difftime(event_Date, first(event_Date), units = "days"))) %>%
    # Ungroup the data frame
    ungroup()
  
  relative_day_df$relative_day = relative_day_df$relative_day + 1
  relative_day_df_sub = relative_day_df %>%
    select(patientid, patientvisitid, relative_day, event_Date_sh)
  kdigo_data_sub = merge(kdigo_data_sub, relative_day_df_sub, by = c("patientid", "patientvisitid", "event_Date_sh"))
  kdigo_data_with_aki = merge(kdigo_data_sub, data, by = c("patientid", "patientvisitid"))
  
  # missingness of KDIGO by day and AKI
  missingness_by_AKI_new <- kdigo_data_with_aki %>%
    filter(relative_day <= 8) %>%
    group_by(patientid, patientvisitid, relative_day, AKI_Category_subgroup) %>%
    summarize(kdigo_all_missing = all(is.na(kdigo)), .groups = 'drop') %>%
    group_by(relative_day, AKI_Category_subgroup) %>%
    summarize(
      total_days = n(),
      missing_days = sum(kdigo_all_missing),
      missing_percentage = (missing_days / total_days) * 100
    ) %>%
    mutate(imp = i)
  
  missingness_new <- kdigo_data_with_aki %>%
    filter(relative_day <= 8) %>%
    group_by(patientid, patientvisitid, relative_day) %>%
    summarize(kdigo_all_missing = all(is.na(kdigo)), .groups = 'drop') %>%
    group_by(relative_day) %>%
    summarize(
      total_days = n(),
      missing_days = sum(kdigo_all_missing),
      missing_percentage = (missing_days / total_days) * 100
    ) %>%
    mutate(imp = i)

  missingness = rbind(missingness, missingness_new)
  missingness_by_AKI = rbind(missingness_by_AKI, missingness_by_AKI_new)
  
  gc()
  gc()
  gc()
  gc()
  gc()
  
}

missingness_pooled = missingness %>%
  filter(relative_day <= 8) %>%
  group_by(relative_day) %>%
  summarise(missingness = mean(missing_percentage))

write.csv(missingness_pooled, "G:/Persistent_AKI/Xinlei/Result/KDIGO_Missing.csv")


missingness_by_AKI_pooled = missingness_by_AKI %>%
  filter(relative_day <= 8) %>%
  group_by(AKI_Category_subgroup, relative_day) %>%
  summarise(missingness = mean(missing_percentage))

write.csv(missingness_by_AKI_pooled, "G:/Persistent_AKI/Xinlei/Result/KDIGO_Missing_by_AKI.csv")
