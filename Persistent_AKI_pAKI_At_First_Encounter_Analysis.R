library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
#library(interval)
#library(icenReg)


# setting
analysis="_pAKI_at_first_encounter_method1.2"
#analysis="_pAKI_at_first_encounter_method1.2_sensitivity"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the model fit
cox_fit_list = list()

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
             "readmission",
             "daydiff",
             "aki_1_icu",
             "AKI_Category")

# if the subject has pAKI in the first encounter
pAKI_first_encounter_df = rawData %>%
  group_by(patientid) %>%
  arrange(encounter_start_date_time_sh) %>%
  filter(row_number() == 1) %>%
  select(patientid, AKI_Category_subgroup) %>%
  mutate(pAKI_at_first_encounter = ifelse(AKI_Category_subgroup=="Persistent_Severe_AKI", 1, 0)) %>%
  select(patientid, pAKI_at_first_encounter)
rawData = merge(rawData, pAKI_first_encounter_df, by=c("patientid"), all.x=TRUE)

# method 1: include only subjects who have >= 2 encounters
## time 0: time of the 2nd encounter
## y1: 1=pAKI, 2=death, 3=none of them by day 90 (multi-nomial logistic regression)
## y2: pAKI or death = 1, none of them = 0 by day 90 (cox model)
# method 2: include everyone
## time 0: time of the first encounter

## find the date of the 2nd encounter
time0_df = rawData %>%
  group_by(patientid) %>%
  arrange(encounter_start_date_time_sh) %>%
  filter(row_number() == 2) %>%
  select(patientid, encounter_start_date_time_sh)
colnames(time0_df)[ncol(time0_df)] = "secound_encounter_date_time"
rawData = merge(rawData, time0_df, by=c("patientid"), all.x=TRUE)

#check = tibble(rawData$patientid, 
#               rawData$patient_encounters, 
#               rawData$encounter_start_date_time_sh, 
#               rawData$readmission, 
#               rawData$daydiff,
#               rawData$secound_encounter_date_time)

# calculate days before death
## if the death date is earlier than the ICU admission date than replace it with NA
rawData$death_date_sh_true = ifelse(difftime(rawData$death_date_sh, rawData$dt_icu_start_sh, units="days")<0, NA, rawData$death_date_sh)
## if the survived subjects' death date is earlier than or the same as the hospital discharge date, then the subjects died in the hospital
rawData$discharge_to_death_days = difftime(rawData$death_date_sh_true, rawData$encounter_end_date_time_sh, unit="days")
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days<=0 & rawData$dead==0), 1, as.numeric(as.character(rawData$dead)))
## if the dead subjects' death date is later than the hospital discharge date, the subjects survived in hospital
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days>0 & rawData$dead==1), 0, as.numeric(as.character(rawData$dead)))

## number of days from 2nd encounter to death
rawData$encounter2_to_death_days = difftime(rawData$death_date_sh_true, rawData$secound_encounter_date_time, units="days")
## number of days from 2nd encounter to hospital discharge
rawData$encounter2_to_discharge_days = difftime(rawData$encounter_end_date_time_sh, rawData$secound_encounter_date_time, units="days")
## replace missing date with discharge date
rawData$survived_days = ifelse((is.na(rawData$death_date_sh_true)), rawData$encounter2_to_discharge_days, rawData$encounter2_to_death_days)

#check = tibble(rawData$patientid, 
#               rawData$patient_encounters, 
#               rawData$encounter_start_date_time_sh, 
#               rawData$encounter_end_date_time_sh, 
#               rawData$readmission, 
#               rawData$secound_encounter_date_time,
#               rawData$survived_days,
#               rawData$death_date_sh_true)

rawData_sub = rawData %>%
  filter(readmission==1)

#check = tibble(rawData_sub$patientid, 
#               rawData_sub$patient_encounters, 
#               rawData_sub$encounter_start_date_time_sh, 
#               rawData_sub$encounter_end_date_time_sh, 
#               rawData_sub$readmission, 
#               rawData_sub$secound_encounter_date_time,
#               rawData_sub$survived_days,
#               rawData_sub$death_date_sh_true)

## 90-day mortality: 
### - if the death date is not NA and the 2nd encounter to death days is >= 90 -- survived
### - if the death date is not NA and the 2nd encounter to death days is < 90 -- dead
### - if the death date is NA that means we didn't observe the event time but we know the result, the data is censored, so the survival = dead_hosp and time = encounter2_to_discharge_days
rawData_sub$dead_90 = ifelse(!(is.na(rawData_sub$death_date_sh_true)) & rawData_sub$survived_days>=90, 0, 
                         ifelse(!(is.na(rawData_sub$death_date_sh_true)) & rawData_sub$survived_days<90, 1, 
                                rawData_sub$dead_hosp))
rawData_sub$survival_time = ifelse(rawData_sub$survived_days>90, 90, rawData_sub$survived_days)


## 2nd encounter to first AKI days
rawData_sub$encounter2_to_aki_1_days = difftime(rawData_sub$aki_1_icu_dt_sh, rawData_sub$secound_encounter_date_time, unit="days")
rawData_sub$encounter2_to_aki_2_3_days = difftime(rawData_sub$aki_2_3_icu_dt_sh, rawData_sub$secound_encounter_date_time, unit="days")

#rawData$encounter2_to_aki_2_3_days = difftime(rawData$aki_2_3_icu_dt_sh, rawData$secound_encounter_date_time, unit="days")
#check = rawData %>%
#  select(patientid, patientvisitid, 
#         encounter_start_date_time_sh, encounter_end_date_time_sh, secound_encounter_date_time,
#         aki_2_3_icu,
#         aki_2_3_icu_dt_sh,
#         encounter2_to_aki_2_3_days,
#         patient_encounters,
#         readmission) %>%
#  filter(patientid==717836)
  # filter(encounter2_to_aki_2_3_days<0)
#sum(check$readmission==1)


## 90-day pAKI
### - for those who has AKI_Category == "Persistent_Severe_AKI", if the encounter2_to_aki_2_3_days >= 90, then pAKI_90 = 0
### - for those who has AKI_Category == "Persistent_Severe_AKI", if the encounter2_to_aki_2_3_days < 90, then pAKI_90 = 1
### - for those who has AKI_Category != "Persistent_Severe_AKI", pAKI_90 = 0
rawData_sub$pAKI_90 = ifelse((rawData_sub$AKI_Category_subgroup == "Persistent_Severe_AKI" & rawData_sub$encounter2_to_aki_2_3_days >= 90), 0,
                          ifelse((rawData_sub$AKI_Category_subgroup == "Persistent_Severe_AKI" & rawData_sub$encounter2_to_aki_2_3_days < 90), 1, 0))
rawData_sub$pAKI_time = ifelse(rawData_sub$pAKI_90==1, rawData_sub$encounter2_to_aki_2_3_days, 90)

check = tibble(rawData_sub$patientid, 
               rawData_sub$patient_encounters, 
               rawData_sub$encounter2_to_aki_2_3_days,
               rawData_sub$AKI_Category,
               rawData_sub$pAKI_time,
               rawData_sub$pAKI_90,
               rawData_sub$survival_time,
               rawData_sub$dead_90)

# the composite outcome and time
rawData_sub$dead_pAKI_90 = ifelse((rawData_sub$dead_90==1 | rawData_sub$pAKI_90 == 1), 1, 0)
rawData_sub$dead_pAKI_time = ifelse((rawData_sub$dead_90==1 & rawData_sub$pAKI_90==0), rawData_sub$survival_time,
                                    ifelse((rawData_sub$dead_90==0 & rawData_sub$pAKI_90==1), rawData_sub$pAKI_time,
                                           ifelse((rawData_sub$dead_90==1 & rawData_sub$pAKI_90==1 & rawData_sub$survival_time<rawData_sub$pAKI_time), rawData_sub$survival_time, 
                                                  ifelse((rawData_sub$dead_90==1 & rawData_sub$pAKI_90==1 & rawData_sub$survival_time>rawData_sub$pAKI_time), rawData_sub$pAKI_time, 90))))

check = rawData_sub %>%
  select(patientid, patient_encounters, 
         pAKI_90, dead_90, dead_pAKI_90, 
         pAKI_time, survival_time, death_date_sh_true,
         dead_pAKI_90,
         dead_pAKI_time)

visitid = unique(rawData_sub$patientvisitid)


for (i in 1:length(impDat)){
  
  # get each imputed dataset
  imp = impDat[[i]]
  # replace the imputed variables
  imp = imp[imp$patientvisitid %in% visitid, ]
  imp = imp %>%
    arrange(patientid, patientvisitid)
  data = rawData_sub %>%
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
  if(analysis == "_pAKI_at_first_encounter_method1.2"){
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
              "dead_pAKI_90",
              "dead_pAKI_time",
              "pAKI_at_first_encounter")
  
  coxdata = data %>%
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    mutate(race = ifelse(race==1, 0, 1)) %>%
    select(all_of(adj_var)) %>%
    mutate(charlson_age = charlson * age_at_admission)
  
  #set.seed(960725)
  #obj = rfsrc(Surv(dead_pAKI_time, dead_pAKI_90) ~ ., coxdata,
  #ntree = 1, nodesize = 5, nsplit = 10, importance = TRUE)
  #tree = get.tree(obj, 1)
  #tree
  
  cox_fit = coxph(Surv(dead_pAKI_time, dead_pAKI_90) ~ ., data=coxdata)
  #cox_fit = icenReg::ic_sp(Surv(time=coxdata$time, time2=coxdata$time2, type="interval2") ~ AKI_Category_true, data=coxdata, model="ph")
  #summary(cox_fit)
  
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
  cox_fit_df_new = data.frame(variable = names(cox_fit$coefficients),
                              imp = i,
                              beta = cox_fit$coefficients,
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
beta_df$Vt = beta_df$Vb + beta_df$Vw + (beta_df$Vb/length(unique(cox_fit_df$imp)))
beta_df$se_pool = sqrt(beta_df$Vt)
beta_df$z = beta_df$beta_bar/beta_df$se_pool
beta_df$p_value = 2*pnorm(-abs(beta_df$z))

beta_df$variable_true = c("Age at admission",
                          "Apache 3",
                          "Charlson",
                          "Charlson * Age at admission",
                          "Gender: male",
                          "Days of mechanical ventilation",
                          "Persistent AKI at the first encounter: yes",
                          "Race: non-white")

beta_df$hazard_ratio = exp(beta_df$beta_bar)

order =  c("Age at admission",
           "Gender: male",
           "Race: non-white",
           "Apache 3",
           "Charlson",
           "Charlson * Age at admission",
           "Days of mechanical ventilation",
           "Persistent AKI at the first encounter: yes")

cox_res_df = beta_df %>%
  mutate(variable_true = factor(variable_true, levels = order)) %>%
  arrange(variable_true) %>%
  select(variable_true, beta_bar, hazard_ratio, se_pool, z, p_value)
colnames(cox_res_df) = c("Variable", "Coefficient", "Hazard Ratio", "Std Error", "Z-Statistics", "P-Value")
write.csv(cox_res_df, paste("G:/Persistent_AKI/Xinlei/Result/cox_res_df", analysis, ".csv", sep = ""))


