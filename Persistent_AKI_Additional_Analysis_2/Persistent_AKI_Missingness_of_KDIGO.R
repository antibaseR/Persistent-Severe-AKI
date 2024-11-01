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

kdigo_data_sub = kdigo_data %>%
  select(patientid, patientvisitid, kdigo, event_Date_sh) %>%
  filter(patientvisitid %in% rawData$patientvisitid)

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

missingness_df = tibble()
missingness_by_group_df = tibble()

KDIGO_variable = c("day1_kidgo", "day2_kidgo", "day3_kidgo", "day4_kidgo", 
                   "day5_kidgo", "day6_kidgo", "day7_kidgo", "day8_kidgo")


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
  
  KDIGO_data = data %>%
    select(patientid, patientvisitid, 
           day1_kidgo, day2_kidgo, day3_kidgo,
           day4_kidgo, day5_kidgo, day6_kidgo,
           day7_kidgo, day8_kidgo, AKI_Category_subgroup) %>%
    filter(AKI_Category_subgroup != "no_AKI")
  
  # Calculate the missingness percentage for each variable
  missingness = sapply(KDIGO_data[KDIGO_variable], function(x) mean(is.na(x)) * 100)
  
  missingness_new = tibble(KDIGO_day = c(1:8),
                                 Missingness = missingness,
                                 imp = i)
  
  # Calculate the missingness percentage for each variable by AKI_Category_subgroup
  missingness_by_group_new = KDIGO_data %>%
    group_by(AKI_Category_subgroup) %>%
    summarise(across(all_of(KDIGO_variable), ~ mean(is.na(.)) * 100, .names = "{col}"))
  missingness_by_group_new$imp = i
  
  missingness_df = rbind(missingness_df, missingness_new)
  missingness_by_group_df = rbind(missingness_by_group_df, missingness_by_group_new)
  
  # missingness of all kdigo
  data_sub = data %>%
    select(patientid, patientvisitid, AKI_Category_subgroup, AKI_Category_subgroup2)
  data_merge = merge(data_sub, kdigo_data_sub, c("patientid", "patientvisitid"))
                                                      
}

missingness_df_pooled = missingness_df %>%
  group_by(KDIGO_day) %>%
  summarise(missingness = mean(Missingness))

missingness_by_group_df_pooled = missingness_by_group_df %>%
  group_by(AKI_Category_subgroup) %>%
  summarise(missingness = across(all_of(KDIGO_variable), ~ mean(.), .names = "{col}"))

