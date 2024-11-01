library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
#library(interval)
#library(icenReg)

analysis = "_mortality_ageq"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# check if the death happened after AKI
#time_to_death = difftime(rawData$death_date_sh, rawData$encounter_start_date_time_sh, unit="secs")
#check = tibble(time_to_death = time_to_death, time_to_AKI = rawData$TIME_TO_FIRST_AKI)
#check$death_after_pAKI = ifelse(check$time_to_death>check$time_to_AKI, 1, 0)
#mean(check$death_after_pAKI==1, na.rm=T) # 99% death happened after the development of AKI

# create list to save the model fit
km_fit_list = list()
logrank_fit_list = list()
cox_fit_list = list()
cox_fit_ageq_list = list()

# create list to save the results
Q0_list = vector()
Q1_list = vector()
Q2_list = vector()
Q3_list = vector()
Q4_list = vector()
n1_list = vector()
n2_list = vector()
n3_list = vector()
n4_list = vector()
d1_list = vector()
d2_list = vector()
d3_list = vector()
d4_list = vector()

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
  if (analysis == "_mortality_ageq"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  # subset the data
  subdata = data[, c(names(data) %in% surv_var)] %>%
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
              "AKI_Category_subgroup",
              "dead_90")
  
  coxdata = subdata[, c(names(subdata) %in% adj_var)]
  
  # find the age quartiles
  Q0 = summary(coxdata$age_at_admission)[1]
  Q1 = summary(coxdata$age_at_admission)[2]
  Q2 = summary(coxdata$age_at_admission)[3]
  Q3 = summary(coxdata$age_at_admission)[5]
  Q4 = summary(coxdata$age_at_admission)[6]
  
  # save as a vector
  age_quartiles = c(Q0, Q1, Q2, Q3, Q4)
  
  n1 = coxdata %>%
    filter(age_at_admission>=Q0 & age_at_admission<Q1) %>%
    nrow()
  n2 = coxdata %>%
    filter(age_at_admission>=Q1 & age_at_admission<Q2) %>%
    nrow()
  n3 = coxdata %>%
    filter(age_at_admission>=Q2 & age_at_admission<Q3) %>%
    nrow()
  n4 = coxdata %>%
    filter(age_at_admission>=Q3 & age_at_admission<=Q4) %>%
    nrow()
  
  d1 = coxdata %>%
    filter(age_at_admission>=Q0 & age_at_admission<Q1) %>%
    filter(dead_90 == 1) %>%
    nrow()
  d2 = coxdata %>%
    filter(age_at_admission>=Q1 & age_at_admission<Q2) %>%
    filter(dead_90 == 1) %>%
    nrow()
  d3 = coxdata %>%
    filter(age_at_admission>=Q2 & age_at_admission<Q3) %>%
    filter(dead_90 == 1) %>%
    nrow()
  d4 = coxdata %>%
    filter(age_at_admission>=Q3 & age_at_admission<=Q4) %>%
    filter(dead_90 == 1) %>%
    nrow()
  
  Q0_list = append(Q0_list, Q0)
  Q1_list = append(Q1_list, Q1)
  Q2_list = append(Q2_list, Q2)
  Q3_list = append(Q3_list, Q3)
  Q4_list = append(Q4_list, Q4)
  
  n1_list = append(n1_list, n1)
  n2_list = append(n2_list, n2)
  n3_list = append(n3_list, n3)
  n4_list = append(n4_list, n4)
  
  d1_list = append(d1_list, d1)
  d2_list = append(d2_list, d2)
  d3_list = append(d3_list, d3)
  d4_list = append(d4_list, d4)

}

mean(Q0_list)
mean(Q1_list)
mean(Q2_list)
mean(Q3_list)
mean(Q4_list)


mean(n1_list)
mean(n2_list)
mean(n3_list)
mean(n4_list)


mean(d1_list)
mean(d2_list)
mean(d3_list)
mean(d4_list)

n1_list + n2_list + n3_list + n4_list

mean(n1_list + n2_list + n3_list + n4_list)

mean(d1_list/n1_list)
mean(d2_list/n2_list)
mean(d3_list/n3_list)
mean(d4_list/n4_list)

chisq.dat = t(matrix(c(mean(d1_list),
                   mean(d2_list),
                   mean(d3_list),
                   mean(d4_list),
                   mean(n1_list)-mean(d1_list),
                   mean(n2_list)-mean(d2_list),
                   mean(n3_list)-mean(d3_list),
                   mean(n4_list)-mean(d4_list)), ncol=2))

chisq.test(chisq.dat)