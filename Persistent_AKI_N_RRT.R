library(tidyverse)
library(caret)
library(pROC)

# setting
analysis="_RRT"
# analysis="_RRT_sensitivity"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the results
rrt_list = vector()
rrt_pAKI_list = vector()
rrt_tAKI_list = vector()
rrt_s1AKI_list = vector()
rrt_nAKI_list = vector()



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

risk_var = c("age_at_admission",
             "gender",
             "charlson",
             "apache3",
             "race",
             "MV_DAYS",
             "refCreat",
             "AKI_Category_final",
             "rrt_discharge")


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
  if (analysis == "_RRT"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  rrt_list[i] = data %>%
    filter(rrt_discharge==1) %>%
    nrow()
  
  rrt_pAKI_list[i] = data %>%
    filter(rrt_discharge==1) %>%
    filter(AKI_Category_subgroup == "Persistent_Severe_AKI") %>%
    nrow()
    
  rrt_tAKI_list[i] = data %>%
    filter(rrt_discharge==1) %>%
    filter(AKI_Category_subgroup == "No_Persistent_Severe_AKI") %>%
    nrow()
    
  rrt_s1AKI_list[i] = data %>%
    filter(rrt_discharge==1) %>%
    filter(AKI_Category_subgroup == "stage_1_AKI") %>%
    nrow()
    
  rrt_nAKI_list[i] = data %>%
    filter(rrt_discharge==1) %>%
    filter(AKI_Category_subgroup == "no_AKI") %>%
    nrow()  
 
}

mean(rrt_list)/190550
mean(rrt_pAKI_list)/8059
mean(rrt_tAKI_list)/57060
mean(rrt_s1AKI_list)/31472
mean(rrt_nAKI_list)/93959
