library(tidyverse)
library(survival)
library(survminer)

# setting
#analysis = c("_sensitivity")
analysis=""

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create a empty dataframe
risk_table_all = data.frame()
risk_table_psAKI = data.frame()
risk_table_npsAKI = data.frame()
risk_table_stage1 = data.frame()
risk_table_noAKI = data.frame()

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
             "race",
             "charlson",
             # comorbidity start
             "mi",
             "chf",
             "pvd",
             "cevd",
             "dementia",
             "cpd",
             "rheumd",
             "pud",
             "mld",
             "diab", 
             "diabwc",
             "diabetes",
             "hp",
             "rend",
             "canc",
             "msld",
             "metacanc",
             'aids',
             # comorbidity end
             "NSAIDS_BEFORE_AKI",
             "ACE_BEFORE_AKI",
             "ARB_BEFORE_AKI",
             "DIURETICS_BEFORE_AKI",
             "CALCINEURIN_BEFORE_AKI",
             "ANTIBIOTICS_BEFORE_AKI",
             "VANCOMYCIN_BEFORE_AKI",
             "OTHER_BEFORE_AKI",
             "apache3",
             "TOTAL_INPUTS_BEFORE_AKI",
             "TOTAL_SALINE_BEFORE_AKI",
             "TOTAL_OUTPUTS_BEFORE_AKI",
             "TOTAL_Urine_Output_BEFORE_AKI",
             "TOTAL_BLOOD_PROD_BEFORE_AKI",
             "TOTAL_ALBUMIN_BEFORE_AKI",
             "CUMULATIVE_BALANCE_BEFORE_AKI",
             "TOTAL_DOSE_BEFORE_AKI_norepinephrine",
             "TOTAL_DOSE_BEFORE_AKI_phenylephrine",
             "TOTAL_DOSE_BEFORE_AKI_vasopressin",
             "TOTAL_DOSE_BEFORE_AKI_dopamine",
             "TOTAL_DOSE_BEFORE_AKI_epinephrine",
             "TOTAL_DOSE_BEFORE_AKI_dobutamine",
             "TOTAL_DOSE_BEFORE_AKI_milrinone",
             "MAX_LACTATE",
             "MAX_CHLORIDE",
             "MAX_SODIUM",
             "MIN_ALBUMIN",
             "MIN_PLATELETS",
             "MV_DAYS",
             "SEPSIS",
             "SEPTIC_SHOCK",
             "CARDIAC_SURG",
             'HEART_TRANSP',
             "LUNG_TRANSP",
             "LIVER_TRANSP",
             "refCreat",
             "FIRST_AKI_STAGE",
             "TIME_TO_FIRST_AKI",
             "BMI",
             "sofa_24",
             "ALL_CKD_ESRD",
             "MAX_DBP",
             "MAX_HR",
             "MIN_SpO2",
             "AVG_RR",
             "MAX_TEMP",
             "MIN_HGB",
             "MAX_TOTAL_BILI",
             "MAX_INR",
             "MAX_WBC",
             "MIN_PH",
             'MAX_PEEP',
             "MIN_CARDIAC_INDEX",
             "MAX_CVP",
             "MIN_CVP",
             "MEAN_DIASTOLIC_PRESS",
             "MIN_SCVO2",
             "MIN_SVO2",
             "BASELINE_eGFR",
             "MIN_PULSE_PRESS_BEFORE_AKI",
             "MIN_MAP_BEFORE_AKI",
             "AKI_Category_final",
             "dead_90")


for (i in 1:length(impDat)){
  
  # get each imputed dataset
  imp = impDat[[i]]
  # replace the imputed variables
  # imp = imp[order(c(imp$patientid, imp$encounter_start_date_time_sh)),]
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
  
  data = data %>%
    filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
  data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter id
  # subset to rows that are not duplicates of a previous encounter for that patient
  data = data[!duplicated(data$patientid),]
  
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
  
  # for everyone
  km.fit = survfit(Surv(survival_time, dead_90) ~ 1, data = subdata)
  summary_data = surv_summary(km.fit)
  risk_table_new = summary_data[c("time", "n.risk")]
  risk_table_new = risk_table_new %>%
    filter(time %in% seq(0, 90, 10))
  risk_table_new$imp = i
  risk_table_all = rbind(risk_table_all, risk_table_new)
  
  # for persistent severe AKI group
  subdata_psAKI = subdata %>%
    filter(AKI_Category_subgroup == "Persistent_Severe_AKI")
  km.fit = survfit(Surv(survival_time, dead_90) ~ 1, data = subdata_psAKI)
  summary_data = surv_summary(km.fit)
  risk_table_new = summary_data[c("time", "n.risk")]
  risk_table_new = risk_table_new %>%
    filter(time %in% seq(0, 90, 10))
  risk_table_new$imp = i
  risk_table_psAKI = rbind(risk_table_psAKI, risk_table_new)
  
  # for no persistent severe AKI group
  subdata_npsAKI = subdata %>%
    filter(AKI_Category_subgroup == "No_Persistent_Severe_AKI")
  km.fit = survfit(Surv(survival_time, dead_90) ~ 1, data = subdata_npsAKI)
  summary_data = surv_summary(km.fit)
  risk_table_new = summary_data[c("time", "n.risk")]
  risk_table_new = risk_table_new %>%
    filter(time %in% seq(0, 90, 10))
  risk_table_new$imp = i
  risk_table_npsAKI = rbind(risk_table_npsAKI, risk_table_new)
  
  # for no stage 1 AKI group
  subdata_stage1 = subdata %>%
    filter(AKI_Category_subgroup == "stage_1_AKI")
  km.fit = survfit(Surv(survival_time, dead_90) ~ 1, data = subdata_stage1)
  summary_data = surv_summary(km.fit)
  risk_table_new = summary_data[c("time", "n.risk")]
  risk_table_new = risk_table_new %>%
    filter(time %in% seq(0, 90, 10))
  risk_table_new$imp = i
  risk_table_stage1 = rbind(risk_table_stage1, risk_table_new)
  
  # for no stage 1 AKI group
  subdata_noAKI = subdata %>%
    filter(AKI_Category_subgroup == "no_AKI")
  km.fit = survfit(Surv(survival_time, dead_90) ~ 1, data = subdata_noAKI)
  summary_data = surv_summary(km.fit)
  risk_table_new = summary_data[c("time", "n.risk")]
  risk_table_new = risk_table_new %>%
    filter(time %in% seq(0, 90, 10))
  risk_table_new$imp = i
  risk_table_noAKI = rbind(risk_table_noAKI, risk_table_new)


}


risk_table_all_sum = risk_table_all %>%
  group_by(time) %>%
  summarise(all_risk = round(mean(n.risk)))

risk_table_psAKI_sum = risk_table_psAKI %>%
  group_by(time) %>%
  summarise(psAKI_risk = round(mean(n.risk)))

risk_table_npsAKI_sum = risk_table_npsAKI %>%
  group_by(time) %>%
  summarise(npsAKI_risk = round(mean(n.risk)))

risk_table_stage1_sum = risk_table_stage1 %>%
  group_by(time) %>%
  summarise(stage1_risk = round(mean(n.risk)))

risk_table_noAKI_sum = risk_table_noAKI %>%
  group_by(time) %>%
  summarise(noAKI_risk = round(mean(n.risk)))

risk_table_sum = risk_table_all_sum
risk_table_sum$psAKI_risk = risk_table_psAKI_sum$psAKI_risk
risk_table_sum$npsAKI_risk = risk_table_npsAKI_sum$npsAKI_risk
risk_table_sum$stage1_risk = risk_table_stage1_sum$stage1_risk
risk_table_sum$noAKI_risk = risk_table_noAKI_sum$noAKI_risk

names(risk_table_sum) = c("Time", "All", "Persistent Severe AKI", "No Persistent Severe AKI", "Stage 1 AKI", "No AKI")

write.csv(risk_table_sum, "./Persistent_AKI_Additional_Analysis/Result/risk_table.csv")
