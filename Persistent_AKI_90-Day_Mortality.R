library(tidyverse)
library(survminer)
library(survival)
library(lubridate)

analysis = "_90day_mortality"

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
  if (analysis == "_90day_mortality"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  
  subdata = data
  
  ## if the death date is earlier than the ICU admission date than replace it with NA
  subdata$death_date_sh_true = ifelse(difftime(subdata$death_date_sh, subdata$dt_icu_start_sh, units="days")<0, NA, subdata$death_date_sh)
  ## if the survived subjects' death date is earlier than or the same as the hospital discharge date, then the subjects died in the hospital
  subdata$discharge_to_death_days = difftime(subdata$death_date_sh_true, subdata$encounter_end_date_time_sh, unit="days")
  subdata$dead_hosp = ifelse((!is.na(subdata$discharge_to_death_days) & subdata$discharge_to_death_days<=0 & subdata$dead==0), 1, as.numeric(as.character(subdata$dead)))
  ## if the dead subjects' death date is later than the hospital discharge date, the subjects survived in hospital
  subdata$dead_hosp = ifelse((!is.na(subdata$discharge_to_death_days) & subdata$discharge_to_death_days>0 & subdata$dead==1), 0, as.numeric(as.character(subdata$dead)))
  
  ## number of days from ICU admission to death
  subdata$icu_to_death_days = difftime(subdata$death_date_sh_true, subdata$dt_icu_start_sh, units="days")
  ## number of days from ICU admission to hospital discharge
  subdata$icu_to_discharge_days = difftime(subdata$encounter_end_date_time_sh, subdata$dt_icu_start_sh, units="days")
  ## replace missing date with discharge date
  subdata$survived_days = ifelse((is.na(subdata$death_date_sh_true)), subdata$icu_to_discharge_days, subdata$icu_to_death_days)
  
  ## 90-day mortality: 
  ### - if the death date is not NA and the icu to death days is >= 90 -- survived
  ### - if the death date is not NA and the icu to death days is < 90 -- dead
  ### - if the death date is NA that means we didn't observe the event time but we know the result, the data is censored, so the survival = dead_hosp and time = icu_to_discharge_days
  subdata$dead_90 = ifelse(!(is.na(subdata$death_date_sh_true)) & subdata$survived_days>=90, 0, 
                           ifelse(!(is.na(subdata$death_date_sh_true)) & subdata$survived_days<90, 1, subdata$dead_hosp))
  subdata$survival_time = ifelse(subdata$survived_days>90, 90, subdata$survived_days)
  
  coxdata = subdata %>% select(dead_90, AKI_Category_subgroup)
  
  all_d = coxdata %>%
    filter(dead_90 == 1) %>%
    nrow()
  
  pAKI_d = coxdata %>%
    filter(AKI_Category_subgroup=="Persistent_Severe_AKI") %>%
    filter(dead_90 == 1) %>%
    nrow()
  
  tAKI_d = coxdata %>%
    filter(AKI_Category_subgroup=="No_Persistent_Severe_AKI") %>%
    filter(dead_90 == 1) %>%
    nrow()
  
  s1AKI_d = coxdata %>%
    filter(AKI_Category_subgroup=="stage_1_AKI") %>%
    filter(dead_90 == 1) %>%
    nrow()
  
  nAKI_d = coxdata %>%
    filter(AKI_Category_subgroup=="no_AKI") %>%
    filter(dead_90 == 1) %>%
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

mean(pAKI_n_list) + mean(tAKI_n_list) + 
  mean(s1AKI_n_list) + mean(nAKI_n_list)

mean(pAKI_d_list+tAKI_d_list+s1AKI_d_list+nAKI_d_list)


mean(all_d_list)/mean(all_n_list)
mean(pAKI_d_list)/mean(pAKI_n_list)
mean(tAKI_d_list)/mean(tAKI_n_list)
mean(s1AKI_d_list)/mean(s1AKI_n_list)
mean(nAKI_d_list)/mean(nAKI_n_list)



