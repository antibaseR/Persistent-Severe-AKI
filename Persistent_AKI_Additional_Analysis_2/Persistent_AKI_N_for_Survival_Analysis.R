library(tidyverse)
library(caret)
library(xgboost)
library(caretEnsemble)
library(pROC)
library(caTools)
library(randomForest)
library(nnet)
library(DALEX)
library(iml)

# setting
for (analysis in c("_survival_sensitivity")){
  
  N = vector()
  N_psAKI = vector()
  
  impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
  rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
  
  # create list to save the model fit
  km_fit_list = list()
  logrank_fit_list = list()
  cox_fit_list = list()
  
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
    
    # subset the data
    if (analysis=="_survival_sensitivity"){
      subdata = data %>%
        filter(!(difftime(dt_icu_end_sh, dt_icu_start_sh, units="days")<=2 & dead==1))
    } else {subdata = data}
    
    subdata = subdata[, c(names(data) %in% surv_var)] %>%
      mutate(TOTAL_BALANCE_BEFORE_AKI = TOTAL_INPUTS_BEFORE_AKI - TOTAL_OUTPUTS_BEFORE_AKI) %>%
      select(-c(TOTAL_INPUTS_BEFORE_AKI,TOTAL_OUTPUTS_BEFORE_AKI))
    
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
    
    # Cox proportional hazard model
    ## variables we need to adjust for the adjusted Cox proportional hazard model
    adj_var = c("age_at_admission",
                "gender",
                "race",
                "charlson",
                "apache3",
                "SEPSIS",
                "MV_DAYS",
                "BMI",
                "TOTAL_BALANCE_BEFORE_AKI",
                "dead_90",
                "survival_time",
                "AKI_Category_subgroup")
    
    coxdata = subdata[, c(names(subdata) %in% adj_var)] %>%
      filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
      mutate(race = ifelse(race==1, 0, 1),
             AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup == "Persistent_Severe_AKI", "Yes", "No"), 
                                            levels=c("No", "Yes")),
             gender_bmi = (as.numeric(as.character(gender)))*BMI)
    
    completeData = coxdata[complete.cases(coxdata), ]
    
    
    N[i] = nrow(completeData) 
    N_psAKI[i] = completeData %>%
      filter(AKI_Category_subgroup == "Yes") %>%
      nrow()
  }
  
  round(mean(N))
  round(mean(N_psAKI))
  
}