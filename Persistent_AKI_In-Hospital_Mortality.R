library(tidyverse)
library(survminer)
library(survival)
library(lubridate)

analysis = "_in_hospital_mortality"

impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the model fit
all_d_list = vector()
pAKI_d_list = vector()
tAKI_d_list = vector()
s1AKI_d_list = vector()
nAKI_d_list = vector()

all_n_list = vector()
pAKI_n_list = vector()
tAKI_n_list = vector()
s1AKI_n_list = vector()
nAKI_n_list = vector()

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
  if (analysis == "_in_hospital_mortality"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  
  coxdata = data %>% select(dead, AKI_Category_subgroup)
  
  all_d = coxdata %>%
    filter(dead == 1) %>%
    nrow()
  
  pAKI_d = coxdata %>%
    filter(AKI_Category_subgroup=="Persistent_Severe_AKI") %>%
    filter(dead == 1) %>%
    nrow()
  
  tAKI_d = coxdata %>%
    filter(AKI_Category_subgroup=="No_Persistent_Severe_AKI") %>%
    filter(dead == 1) %>%
    nrow()
  
  s1AKI_d = coxdata %>%
    filter(AKI_Category_subgroup=="stage_1_AKI") %>%
    filter(dead == 1) %>%
    nrow()
  
  nAKI_d = coxdata %>%
    filter(AKI_Category_subgroup=="no_AKI") %>%
    filter(dead == 1) %>%
    nrow()
  
  all_n = coxdata %>%
    nrow()
  
  pAKI_n = coxdata %>%
    filter(AKI_Category_subgroup=="Persistent_Severe_AKI") %>%
    nrow()
  
  tAKI_n = coxdata %>%
    filter(AKI_Category_subgroup=="No_Persistent_Severe_AKI") %>%
    nrow()
  
  s1AKI_n = coxdata %>%
    filter(AKI_Category_subgroup=="stage_1_AKI") %>%
    nrow()
  
  nAKI_n = coxdata %>%
    filter(AKI_Category_subgroup=="no_AKI") %>%
    nrow()
  
  all_d_list[i] = all_d
  pAKI_d_list[i] = pAKI_d
  tAKI_d_list[i] = tAKI_d
  s1AKI_d_list[i] = s1AKI_d
  nAKI_d_list[i] = nAKI_d
  
  all_n_list[i] = all_n
  pAKI_n_list[i] = pAKI_n
  tAKI_n_list[i] = tAKI_n
  s1AKI_n_list[i] = s1AKI_n
  nAKI_n_list[i] = nAKI_n
  
}

mean(all_d_list)
mean(pAKI_d_list)
mean(tAKI_d_list)
mean(s1AKI_d_list)
mean(nAKI_d_list)

mean(all_n_list)
mean(pAKI_n_list)
mean(tAKI_n_list)
mean(s1AKI_n_list)
mean(nAKI_n_list)


mean(all_d_list)/mean(all_n_list)
mean(pAKI_d_list)/mean(pAKI_n_list)
mean(tAKI_d_list)/mean(tAKI_n_list)
mean(s1AKI_d_list)/mean(s1AKI_n_list)
mean(nAKI_d_list)/mean(nAKI_n_list)



