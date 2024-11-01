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
readmission_n_list = vector()
readmission_nAKI_list = vector()
readmission_stage1AKI_list = vector()
readmission_tAKI_list = vector()
readmission_pAKI_list = vector()

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
             "HOSP_LOS",
             "AKI_Category")



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

readmission_df_new = rawData %>%
  arrange(patientid, encounter_start_date_time_sh) %>%
  group_by(patientid) %>%
  summarise(patientvisitid = last(99999999)) %>%
  mutate(encounter_start_date_time_sh=NA, readmission2=0)

readmission_df_orig = rawData %>% 
  arrange(patientid, encounter_start_date_time_sh) %>%
  select(patientid, patientvisitid, encounter_start_date_time_sh, readmission2)

readmission_df = rbind(readmission_df_orig, readmission_df_new) %>%
  arrange(patientid, encounter_start_date_time_sh) %>%
  group_by(patientid) %>%
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
    #data = data[!duplicated(data$patientid),]
  }

  # Cox proportional hazard model
  ## variables we need to adjust for the adjusted Cox proportional hazard model
  adj_var = c("age_at_admission",
              "gender",
              "race",
              "charlson",
              "apache3",
              "MV_DAYS",
              "HOSP_LOS",
              "AKI_Category_final")
  
  coxdata = data %>%
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    filter(death_at_first_encounter==0) %>%
    mutate(AKI_Category_final = as.numeric(factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "1", "0")))-1,
           race = ifelse(race==1, 0, 1),
           AKI_charlson = AKI_Category_final * charlson,
           AKI_apache3 = AKI_Category_final * apache3,
           AKI_Category_final = as.factor(AKI_Category_final)) %>%
    select(status, time, all_of(adj_var), AKI_charlson, AKI_apache3)
  
  covs = coxdata %>%
    select(all_of(adj_var), AKI_charlson, AKI_apache3) %>%
    as.data.frame()
  
  cox_fit = crr(ftime=coxdata$time, fstatus=coxdata$status, cov1=covs, failcode=1, cencode=0)
  summary(cox_fit)
  
  # append the fitted object to the list
  cox_fit_list[[i]] = cox_fit
}

# save the list as rds
saveRDS(cox_fit_list, paste("G:/Persistent_AKI/Xinlei/Data/cox_fit_list", analysis, ".rds", sep = ""))





# read the saved results
cox_fit_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/cox_fit_list", analysis, ".rds", sep = ""))




# combine the results from the Cox ph model
cox_fit_df = data.frame()
for (i in 1:length(cox_fit_list)){
  cox_fit = cox_fit_list[[i]]
  cox_fit_df_new = data.frame(variable = names(cox_fit$coef),
                              imp = i,
                              beta = cox_fit$coef,
                              se = sqrt(diag(cox_fit$var)))
  cox_fit_df = rbind(cox_fit_df, cox_fit_df_new)
}
rownames(cox_fit_df) = seq(1, nrow(cox_fit_df), 1)

beta_df = cox_fit_df %>%
  group_by(variable) %>%
  summarise(beta_bar = mean(beta))

cox_fit_df = merge(cox_fit_df, beta_df, by="variable")

Vw_df = cox_fit_df %>% 
  group_by(variable) %>%
  summarise(Vw = mean(se^2))

Vb_df = cox_fit_df %>%
  group_by(variable) %>%
  summarise(Vb = sum((beta-beta_bar)^2)/(length(unique(cox_fit_df$imp))-1))

beta_df = merge(merge(beta_df, Vw_df, by="variable"), Vb_df, by="variable")
beta_df$Vt = beta_df$Vb + beta_df$Vw
beta_df$se_pool = sqrt(beta_df$Vt)
beta_df$z = beta_df$beta_bar/beta_df$se_pool
beta_df$p_value = 2*pnorm(-abs(beta_df$z))

beta_df$variable_true = c("Age at admission",
                          "Persistent severe AKI: yes * Apache 3",
                          "Persistent severe AKI: yes",
                          "Persistent severe AKI: yes * Charlson",
                          "Apache 3",
                          "Charlson",
                          "Gender: male",
                          "Hospital length of stay",
                          "Days of mechanical ventilation",
                          "Race: non-white")

beta_df$hazard_ratio = exp(beta_df$beta_bar)

order = c("Age at admission",
          "Gender: male",
          "Race: non-white",
          "Apache 3",
          "Charlson",
          "Hospital length of stay",
          "Days of mechanical ventilation",
          "Persistent severe AKI: yes",
          "Persistent severe AKI: yes * Apache 3",
          "Persistent severe AKI: yes * Charlson")

cox_res_df = beta_df %>%
  mutate(variable_true = factor(variable_true, levels = order)) %>%
  arrange(variable_true) %>%
  select(variable_true, beta_bar, hazard_ratio, se_pool, z, p_value)
colnames(cox_res_df) = c("Variable", "Coefficient", "Hazard Ratio", "Std Error", "Z-Statistics", "P-Value")
write.csv(cox_res_df, paste("G:/Persistent_AKI/Xinlei/Result/cox_res_df", analysis, ".csv", sep = ""))


