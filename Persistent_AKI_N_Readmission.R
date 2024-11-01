library(tidyverse)
library(survival)
library(cmprsk)

# setting
analysis="_readmission"
#analysis="_readmission_sensitivity"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
rawData_true = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
check2 = rawData %>% select(patientid, patientvisitid, encounter_start_date_time_sh, encounter_end_date_time_sh, readmission2, daydiff) %>% filter(patientid == 62021)

# create list to save the model fit
cox_fit_list = list()

# create list to save the results
n_list = vector()
nAKI_list = vector()
stage1AKI_list = vector()
tAKI_list = vector()
pAKI_list = vector()

readmission_n_list = vector()
readmission_nAKI_list = vector()
readmission_stage1AKI_list = vector()
readmission_tAKI_list = vector()
readmission_pAKI_list = vector()

status1_list = vector()
status2_list = vector()
status0_list = vector()

readmission_status1_list = vector()
readmission_status2_list = vector()
readmission_status0_list = vector()

survivor_pAKI_first_encounter_list = vector()
survivor_tAKI_first_encounter_list = vector()

pAKI_first_encounter_readmitted_list = vector()
tAKI_first_encounter_readmitted_list = vector()

pAKI_first_encounter_readmitted_patient_list = vector()
tAKI_first_encounter_readmitted_patient_list = vector()

first_encounter_n_list = vector()
first_encounter_pAKI_n_list = vector()
first_encounter_tAKI_n_list = vector()

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
surv_var = c("patientid",
             "patientvisitid",
             "encounter_start_date_time_sh",
             "encounter_end_date_time_sh", 
             "age_at_admission",
             "gender",
             "race",
             "charlson",
             "apache3",
             "MV_DAYS",
             "readmission2",
             "daydiff",
             "aki_1_icu",
             "AKI_Category_subgroup")



# calculate days before death
## if the death date is earlier than the encounter date than replace it with NA
rawData$death_date_sh_true = ifelse(difftime(rawData$death_date_sh, rawData$encounter_start_date_time_sh, units="days")<0, NA, rawData$death_date_sh)
## if the survived subjects' death date is earlier than or the same as the hospital discharge date, then the subjects died in the hospital
rawData$discharge_to_death_days = difftime(rawData$death_date_sh_true, rawData$encounter_end_date_time_sh, unit="days")
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days<=0 & rawData$dead==0), 1, as.numeric(as.character(rawData$dead)))
## if the dead subjects' death date is later than the hospital discharge date, the subjects survived in hospital
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days>0 & rawData$dead==1), 0, as.numeric(as.character(rawData$dead)))

## number of days from hospital discharge to death
rawData$discharge_to_death_days = difftime(rawData$death_date_sh_true, rawData$encounter_end_date_time_sh, units="days")
## number of days from hospital admission to hospital discharge
rawData$discharge_to_discharge_days = difftime(rawData$encounter_end_date_time_sh, rawData$encounter_end_date_time_sh, units="days")
## replace missing date with discharge date
rawData$survived_days = ifelse((is.na(rawData$death_date_sh_true)), rawData$discharge_to_discharge_days, rawData$discharge_to_death_days)

check = rawData %>%
  select(patientid, patientvisitid, 
         encounter_start_date_time_sh, encounter_end_date_time_sh, death_date_sh_true,
         dead_hosp,
         discharge_to_death_days, discharge_to_discharge_days, survived_days)


## 30-day mortality: 
### - if the death date is not NA and the discharge to death days is >= 30 -- survived
### - if the death date is not NA and the discharge to death days is < 30 -- dead
### - if the death date is NA that means we didn't observe the event time but we know the result, the data is censored, so the survival = dead_hosp and time = encounter2_to_discharge_days
rawData$dead_30 = ifelse(!(is.na(rawData$death_date_sh_true)) & rawData$survived_days>=30, 0, 
                         ifelse(!(is.na(rawData$death_date_sh_true)) & rawData$survived_days<30, 1, 
                                rawData$dead_hosp))
rawData$survival_time = ifelse(rawData$survived_days>30, 30, rawData$survived_days)

check = rawData %>%
  select(patientid, patientvisitid, death_date_sh_true, dead_30, survived_days, survival_time)

## discharge to readmission days
## for each subject, remove it's first record of readmission days and add a new row of NA for that
daydiff_df = rawData %>%
  arrange(patientid, encounter_start_date_time_sh) %>%
  select(patientid, patientvisitid, encounter_start_date_time_sh, encounter_end_date_time_sh, daydiff) %>%
  group_by(patientid) %>%
  summarise(across(.fns = ~ c(., NA))) %>%
  filter(row_number() > 1) %>%
  ungroup()

rawData = rawData %>%
  arrange(patientid, encounter_start_date_time_sh) %>%
  mutate(daydiff_new = daydiff_df$daydiff,
         readmission_days = ifelse((daydiff_new>=30|is.na(daydiff_new>=30)), 30, daydiff_new),
         readmission_days = ifelse(readmission_days<0, 30, readmission_days),
         readmission = ifelse(readmission_days<30, 1, 0))

## subjects who died cannot be readmitted (survival_time < daydiff_new)
## subjects who are readmitted can died later (survival_time > daydiff_new)
check = rawData %>%
  select(patientid, encounter_start_date_time_sh, encounter_end_date_time_sh, patient_encounters, readmission, daydiff_new) %>%
  filter(readmission==0)

check = rawData %>% 
  select(patientid, patientvisitid, patient_encounters, encounter_start_date_time_sh, encounter_end_date_time_sh, death_date_sh_true,readmission, readmission_days, dead_30, survival_time) %>%
  filter(readmission==1 & dead_30 == 1) %>%
  filter(readmission_days > survival_time)

# get the competing risk event and the event time
rawData = rawData %>%
  mutate(status = ifelse((readmission==1 & dead_30==0), 1, 
                         ifelse((readmission==0 & dead_30==1), 2,
                                ifelse((readmission==0 & dead_30==0), 0, 1))))



check = rawData %>%
  filter(readmission==1 & dead_30==1 & survival_time<readmission_days) %>%
  select(patientid, patientvisitid, encounter_start_date_time_sh, encounter_end_date_time_sh,
         death_date_sh_true, readmission, readmission_days, dead_30, survival_time, status)


check_id = rawData %>%
  filter(patientid %in% check$patientid) %>%
  select(patientid, patientvisitid, encounter_start_date_time_sh, encounter_end_date_time_sh,
         death_date_sh_true, readmission, readmission_days, dead_30, survival_time, status)

check_id_true = rawData_true %>%
  filter(patientid %in% check$patientid) %>%
  select(patientid, patientvisitid, encounter_start_date_time_sh, encounter_end_date_time_sh,
         readmission)


rawData = rawData %>%
  mutate(time = ifelse(status==1, readmission_days,
                       ifelse(status==2, survival_time, 30)))

check = rawData %>%
  select(patientid, patientvisitid, dead_30, readmission, survival_time, readmission_days, status, time)

# find the first encounter of each subject, and if the subject died within the first encounter
death_first_encounter_df = rawData %>%
  group_by(patientid) %>%
  arrange(encounter_start_date_time_sh) %>%
  filter(row_number() == 1) %>%
  mutate(death_at_first_encounter = dead_hosp) %>%
  select(patientid, death_at_first_encounter)
rawData = merge(rawData, death_first_encounter_df, by=c("patientid"), all.x=TRUE)


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
  if (analysis == "_readmission"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  # dataframe only with the first encounter
  #first_encounter_df = data %>%
    #filter(death_at_first_encounter==0) %>%
    #group_by(patientid) %>%
    #arrange(patientid, encounter_start_date_time_sh) %>%
    #filter(row_number() == 1)
  first_encounter_df = data[order(data$patientvisitid),] #order dataframe by encounter id
  #first_encounter_df = first_encounter_df %>% filter(death_at_first_encounter==0)
  first_encounter_df = first_encounter_df[!duplicated(first_encounter_df$patientid),]
  
  first_encounter_n = first_encounter_df %>% nrow()
  first_encounter_n_list = append(first_encounter_n_list, first_encounter_n)
  
  first_encounter_pAKI = first_encounter_df %>%
    filter(AKI_Category_subgroup == "Persistent_Severe_AKI")
  first_encounter_pAKI_n = first_encounter_pAKI %>% nrow()
  first_encounter_pAKI_n_list = append(first_encounter_pAKI_n_list, first_encounter_pAKI_n)
  
  first_encounter_tAKI = first_encounter_df %>%
    filter(AKI_Category_subgroup == "No_Persistent_Severe_AKI")
  first_encounter_tAKI_n = first_encounter_tAKI %>% nrow()
  first_encounter_tAKI_n_list = append(first_encounter_tAKI_n_list, first_encounter_tAKI_n)
  
  pAKI_first_encounter_readmitted = data %>%
    filter(patientid %in% first_encounter_pAKI$patientid) %>%
    filter(readmission==1) %>%
    nrow()
  
  tAKI_first_encounter_readmitted = data %>%
    filter(patientid %in% first_encounter_tAKI$patientid) %>%
    filter(readmission==1) %>%
    nrow()
  
  pAKI_first_encounter_readmitted_patient = data %>%
    filter(patientid %in% first_encounter_pAKI$patientid) %>%
    filter(readmission==1) %>%
    select(patientid) %>%
    unique() %>%
    nrow()
  
  tAKI_first_encounter_readmitted_patient = data %>%
    filter(patientid %in% first_encounter_tAKI$patientid) %>%
    filter(readmission==1) %>%
    select(patientid) %>%
    unique() %>%
    nrow()
  
  pAKI_first_encounter_readmitted_list = append(pAKI_first_encounter_readmitted_list, pAKI_first_encounter_readmitted)
  tAKI_first_encounter_readmitted_list = append(tAKI_first_encounter_readmitted_list, tAKI_first_encounter_readmitted)
  pAKI_first_encounter_readmitted_patient_list = append(pAKI_first_encounter_readmitted_patient_list, pAKI_first_encounter_readmitted_patient)
  tAKI_first_encounter_readmitted_patient_list = append(tAKI_first_encounter_readmitted_patient_list, tAKI_first_encounter_readmitted_patient)
  
  
  # only include those who survived the first encounter
  data = data %>%
    filter(death_at_first_encounter==0)
  
  n = data %>% nrow()
  nAKI = data %>% filter(AKI_Category_subgroup == "no_AKI") %>% nrow()
  stage1AKI = data %>% filter(AKI_Category_subgroup == "stage_1_AKI") %>% nrow()
  tAKI = data %>% filter(AKI_Category_subgroup == "No_Persistent_Severe_AKI") %>% nrow()
  pAKI = data %>% filter(AKI_Category_subgroup == "Persistent_Severe_AKI") %>% nrow()
  
  readmission_n = data %>% filter(readmission==1) %>% nrow()
  readmission_nAKI = data %>% filter(readmission==1) %>% filter(AKI_Category_subgroup == "no_AKI") %>% nrow()
  readmission_stage1AKI = data %>% filter(readmission==1) %>% filter(AKI_Category_subgroup == "stage_1_AKI") %>% nrow()
  readmission_tAKI = data %>% filter(readmission==1) %>% filter(AKI_Category_subgroup == "No_Persistent_Severe_AKI") %>% nrow()
  readmission_pAKI = data %>% filter(readmission==1) %>% filter(AKI_Category_subgroup == "Persistent_Severe_AKI") %>% nrow()
  
  n_list = append(n_list, n)
  nAKI_list = append(nAKI_list, nAKI)
  stage1AKI_list = append(stage1AKI_list, stage1AKI)
  tAKI_list = append(tAKI_list, tAKI)
  pAKI_list = append(pAKI_list, pAKI)
  
  readmission_n_list = append(readmission_n_list, readmission_n)
  readmission_nAKI_list = append(readmission_nAKI_list, readmission_nAKI)
  readmission_stage1AKI_list = append(readmission_stage1AKI_list, readmission_stage1AKI)
  readmission_tAKI_list = append(readmission_tAKI_list, readmission_tAKI)
  readmission_pAKI_list = append(readmission_pAKI_list, readmission_pAKI)
  
  status1 = data %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
    filter(status==1) %>%
    #filter(readmission==1) %>%
    nrow()
  
  status2 = data %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
    filter(status==2) %>%
    nrow()
  
  status0 = data %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
    filter(status==0) %>%
    nrow()
  
  readmission_status1 = data %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
    filter(readmission==1) %>%
    filter(status==1) %>%
    nrow()
  
  readmission_status2 = data %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
    filter(readmission==1) %>%
    filter(status==2) %>%
    nrow()
  
  readmission_status0 = data %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
    filter(readmission==1) %>%
    filter(status==0) %>%
    nrow()
  
  status0_list = append(status0_list, status0)
  status1_list = append(status1_list, status1)
  status2_list = append(status2_list, status2)
  
  readmission_status0_list = append(readmission_status0_list, readmission_status0)
  readmission_status1_list = append(readmission_status1_list, readmission_status1)
  readmission_status2_list = append(readmission_status2_list, readmission_status2)
  
  data_survivor_AKI =  data %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) 
  survivor_AKI_id = data_survivor_AKI$patientid
  
  data_survivor_all_encounter = data %>%
    filter(patientid %in% survivor_AKI_id) %>%
    arrange(patientid, encounter_start_date_time_sh)
  
  survivor_first_encounter_data = data_survivor_all_encounter %>%
    group_by(patientid) %>%
    arrange(patientid, encounter_start_date_time_sh) %>%
    filter(row_number() == 1)
  
  survivor_pAKI_first_encounter = survivor_first_encounter_data %>%
    filter(AKI_Category_subgroup=="Persistent_Severe_AKI") %>%
    nrow()
  
  survivor_tAKI_first_encounter = survivor_first_encounter_data %>%
    filter(AKI_Category_subgroup=="No_Persistent_Severe_AKI") %>%
    nrow()
  
  survivor_pAKI_first_encounter_list = append(survivor_pAKI_first_encounter_list, survivor_pAKI_first_encounter)
  survivor_tAKI_first_encounter_list = append(survivor_tAKI_first_encounter_list, survivor_tAKI_first_encounter)
  
  
  
  
}


mean(n_list)
mean(nAKI_list)
mean(stage1AKI_list)
mean(tAKI_list)
mean(pAKI_list)

mean(readmission_n_list)
mean(readmission_nAKI_list)
mean(readmission_stage1AKI_list)
mean(readmission_tAKI_list)
mean(readmission_pAKI_list)

mean(readmission_tAKI_list)/mean(tAKI_list)
mean(readmission_pAKI_list)/mean(pAKI_list)

mean(status0_list)
mean(status1_list)
mean(status2_list)

mean(status0_list)+mean(status1_list)+mean(status2_list)

mean(readmission_status0_list)
mean(readmission_status1_list)
mean(readmission_status2_list)

mean(readmission_status0_list) + mean(readmission_status1_list) + mean(readmission_status2_list)

# how many among 81K having Persistent/Transient in the first encounter
mean(survivor_pAKI_first_encounter_list)
mean(survivor_tAKI_first_encounter_list)

# number of patients in the first encounter by AKI
mean(first_encounter_n_list)
mean(first_encounter_pAKI_n_list)
mean(first_encounter_tAKI_n_list)


# what proportion of patients that has pAKI during the 1st admission ended up readmitted
mean(pAKI_first_encounter_readmitted_list)
mean(pAKI_first_encounter_readmitted_patient_list)/mean(first_encounter_pAKI_n_list)
# what proportion of patients that has tAKI during the 1st admission ended up readmitted
mean(tAKI_first_encounter_readmitted_list)
mean(tAKI_first_encounter_readmitted_patient_list)/mean(first_encounter_tAKI_n_list)
