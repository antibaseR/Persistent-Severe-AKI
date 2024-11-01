library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
library(ggalluvial)

analysis = "_333"

sgData = read.csv("G:/Persistent_AKI/Persistent_AKI_Data/AKI_Subgroups.csv")
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the model fit
n333_list = vector()

# create list to save the results
contrast_list = c("Persistent_Transient", 
                  "Persistent_No",
                  "Persistent_Stage1",
                  "Persistent_No_and_Stage1",
                  "Persistent_All",
                  "Transient_No",
                  "Transient_Stage1",
                  "Transient_No_and_Stage1")

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

# variables we need in the analysis
surv_var = c("age_at_admission",
             "gender",
             "race",
             "charlson",
             "vp_d",
             "apache3",
             "sofa_24",
             "SEPSIS",
             "CARDIAC_SURG",
             "HEART_TRANSP",
             "LUNG_TRANSP",
             "LIVER_TRANSP",
             'MAX_LACTATE',
             "MV_DAYS",
             "rrt_discharge",
             "weight",
             "BMI",
             "TOTAL_INPUTS_BEFORE_AKI",
             "TOTAL_OUTPUTS_BEFORE_AKI",
             "dead",
             "death_date_sh",
             "encounter_start_date_time_sh",
             "encounter_end_date_time_sh",
             "dt_icu_start_sh",
             "FIRST_AKI_DATETIME",
             "AKI_Category",
             "AKI_Category_true",
             "AKI_Category_subgroup",
             "CKD")


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
  
  # merge the subgroup data
  sgData_new = sgData %>%
    select(-c("dt_icu_start_sh", "dt_icu_end_sh", "AKI_Category", "aki_2_3_icu", "aki_1_icu"))
  data = merge(data, sgData_new, by=c("patientid", "patientvisitid"), all.x=T)
  
  ## exclusion
  #The exclusions we need to have for the persistent AKI primary analysis are:
  #1. Kidney transplant
  #2. End stage Renal disease (by ICD-9/10 codes or ICD9/10 codes equivalent to CKD stage 4 or 5)
  #3. eGFR <15 at admission
  #4. Creatinine ≥ 4 mg/dL at admission
  #5. ECMO
  #6. Only use the first encounter
  if (analysis == "_333"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  # reformat the outcome
  data$AKI_Category = as.character(data$AKI_Category)
  data$AKI_Category_true = ifelse((data$aki_1_icu==1 & data$AKI_Category=="no_AKI"), "stage_1_AKI", data$AKI_Category)
  data$AKI_Category = factor(data$AKI_Category, 
                             levels = c("no_AKI", "Transient_AKI", "Persistent_AKI"))
  data$AKI_Category_true = factor(data$AKI_Category_true, 
                                  levels = c("no_AKI", "stage_1_AKI", "Transient_AKI", "Persistent_AKI")) 
  
  #check = data %>%
  #  select(AKI_Category_true, subgroup, AKI_Category_subgroup) %>%
  #  filter(AKI_Category_true == "no_AKI" | AKI_Category_true == "Persistent_AKI")
  data = data %>%
    mutate(AKI_Category_subgroup = ifelse(is.na(subgroup), as.character(AKI_Category_true), subgroup))
  
  data_n = data %>%
    select(day1_kidgo, day2_kidgo, day3_kidgo, day4_kidgo,
           day5_kidgo, day6_kidgo, day7_kidgo,
           AKI_Category_subgroup, AKI_Category)
  
  # create a variable for the KDIGO sequence
  data_n$KDIGO_pattern = paste0(data_n$day1_kidgo, data_n$day2_kidgo, data_n$day3_kidgo, data_n$day4_kidgo, 
                                data_n$day5_kidgo, data_n$day6_kidgo, data_n$day7_kidgo)
  
  subdata = data_n %>%
    filter(AKI_Category=="Transient_AKI")
  subdata$three_three_three = grepl("333", subdata$KDIGO_pattern)
  
  n333_list[i]=sum(grepl("333", subdata$KDIGO_pattern))
  
  write.csv(subdata, paste("G:/Persistent_AKI/Xinlei/kdigo_333_data_", i, ".csv"))
  
}
