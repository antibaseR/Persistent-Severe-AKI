library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
#library(interval)
#library(icenReg)


# setting
analysis="_intervention_by_each"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the model fit
cox_fit_list = list()
cox_fit_int_list = list()

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
             "charlson",
             "apache3",
             "MV_DAYS",
             "refCreat",
             "AKI_Category_subgroup",
             "NEPHROTOXIN_INTERVENTION",
             "DIURETIC_INTERVENTION",
             "VASOPRESSOR_INTERVENTION",
             "SALINE_INTERVENTION")


## variables we need to adjust for the adjusted Cox proportional hazard model
adj_var = list()
adj_var_true = list()
adj_var_order = list()


vaso_cov = c("age_at_admission",
             "gender",
             "apache3",
             "MIN_MAP_BEFORE_AKI",
             "SEPTIC_SHOCK",
             "TOTAL_INPUTS_BEFORE_AKI",
             "dead_pAKI_90",
             "dead_pAKI_time")

vaso_cov_true = c("Age at admission",
                  "Apache3",
                  "Gender: male",
                  "Minimum MAP before AKI",
                  "Septiic shock: yes",
                  "Total fluid input before AKI",
                  "Early vasopressor/inotrope"
)

vaso_cov_order = c("Age at admission",
                   "Gender: male",
                   "Apache3",
                   "Minimum MAP before AKI",
                   "Septiic shock: yes",
                   "Total fluid input before AKI",
                   "Early vasopressor/inotrope"
)


saline_cov = c("age_at_admission",
               "gender",
               "apache3",
               "MIN_MAP_BEFORE_AKI",
               "TOTAL_INPUTS_BEFORE_AKI",
               "dead_pAKI_90",
               "dead_pAKI_time")

saline_cov_true = c("Age at admission",
                    "Apache 3",
                    "Gender: male",
                    "Minimum MAP before AKI",
                    "Saline avoidance before AKI",
                    "Total fluid input before AKI")

saline_cov_order = c("Age at admission",
                     "Gender: male",
                     "Apache 3",
                     "Minimum MAP before AKI",
                     "Total fluid input before AKI",
                     "Saline avoidance before AKI")

diu_cov = c("age_at_admission",
            "gender",
            "apache3",
            "chf",
            "mld",
            "msld",
            "dead_pAKI_90",
            "dead_pAKI_time"
)

diu_cov_true = c("Age at admission",
                 "Apache 3",
                 "Congestive heart failure: yes",
                 "Diuretic avoidance before AKI",
                 "Gender: male",
                 "Mild liver disease: yes",
                 "Moderate or severe liver disease: yes"
)

diu_cov_order = c("Age at admission",
                  "Gender: male",
                  "Apache 3",
                  "Congestive heart failure: yes",
                  "Mild liver disease: yes",
                  "Moderate or severe liver disease: yes",
                  "Diuretic avoidance before AKI")

nephro_cov = c("age_at_admission",
               "gender",
               "apache3",
               "mv",
               "dead_pAKI_90",
               "dead_pAKI_time")

nephro_cov_true = c("Age at admission",
                    "Apache 3",
                    "Gender: male",
                    "Machanical ventilation: yes",
                    "Nephrotoxin avoidance before AKI")

nephro_cov_order = c("Age at admission",
                     "Gender: male",
                     "Apache 3",
                     "Machanical ventilation: yes",
                     "Nephrotoxin avoidance before AKI")

adj_var[[1]] = nephro_cov
adj_var[[2]] = diu_cov
adj_var[[3]] = vaso_cov
adj_var[[4]] = saline_cov

adj_var_true[[1]] = nephro_cov_true
adj_var_true[[2]] = diu_cov_true
adj_var_true[[3]] = vaso_cov_true
adj_var_true[[4]] = saline_cov_true

adj_var_order[[1]] = nephro_cov_order
adj_var_order[[2]] = diu_cov_order
adj_var_order[[3]] = vaso_cov_order
adj_var_order[[4]] = saline_cov_order

int_var = c("NEPHROTOXIN_INTERVENTION",
            "DIURETIC_INTERVENTION",
            "VASOPRESSOR_INTERVENTION",
            "SALINE_INTERVENTION")

int_var_label = c("Nephrotoxin avoidance before AKI",
                  "Diuretic aviodance before AKI",
                  "Early vasopressor/inotrope before AKI",
                  "Saline avoidance before AKI")


# calculate days before death
## if the death date is earlier than the ICU admission date than replace it with NA
rawData$death_date_sh_true = ifelse(difftime(rawData$death_date_sh, rawData$dt_icu_start_sh, units="days")<0, NA, rawData$death_date_sh)
## if the survived subjects' death date is earlier than or the same as the hospital discharge date, then the subjects died in the hospital
rawData$discharge_to_death_days = difftime(rawData$death_date_sh_true, rawData$encounter_end_date_time_sh, unit="days")
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days<=0 & rawData$dead==0), 1, as.numeric(as.character(rawData$dead)))
## if the dead subjects' death date is later than the hospital discharge date, the subjects survived in hospital
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days>0 & rawData$dead==1), 0, as.numeric(as.character(rawData$dead)))

## number of days from ICU admission to death
rawData$icu_to_death_days = difftime(rawData$death_date_sh_true, rawData$dt_icu_start_sh, units="days")
## number of days from ICU admission to hospital discharge
rawData$icu_to_discharge_days = difftime(rawData$encounter_end_date_time_sh, rawData$dt_icu_start_sh, units="days")
## replace missing date with discharge date
rawData$survived_days = ifelse((is.na(rawData$death_date_sh_true)), rawData$icu_to_discharge_days, rawData$icu_to_death_days)


## 90-day mortality: 
### - if the death date is not NA and the icu to death days is >= 90 -- survived
### - if the death date is not NA and the icu to death days is < 90 -- dead
### - if the death date is NA that means we didn't observe the event time but we know the result, the data is censored, so the survival = dead_hosp and time = icu_to_discharge_days
rawData$dead_90 = ifelse(!(is.na(rawData$death_date_sh_true)) & rawData$survived_days>=90, 0, 
                         ifelse(!(is.na(rawData$death_date_sh_true)) & rawData$survived_days<90, 1, rawData$dead_hosp))
rawData$survival_time = ifelse(rawData$survived_days>90, 90, rawData$survived_days)


## ICU admission to first AKI days
rawData$icu_to_aki_1_days = difftime(rawData$aki_1_icu_dt_sh, rawData$dt_icu_start_sh, unit="days")
rawData$icu_to_aki_2_3_days = difftime(rawData$aki_2_3_icu_dt_sh, rawData$dt_icu_start_sh, unit="days")


## 90-day pAKI
### - for those who has AKI_Category == "Persistent_Severe_AKI", if the icu_to_aki_2_3_days >= 90, then pAKI_90 = 0
### - for those who has AKI_Category == "Persistent_Severe_AKI", if the icu_to_aki_2_3_days < 90, then pAKI_90 = 1
### - for those who has AKI_Category != "Persistent_Severe_AKI", pAKI_90 = 0
rawData$pAKI_90 = ifelse((rawData$AKI_Category_subgroup == "Persistent_Severe_AKI" & rawData$icu_to_aki_2_3_days >= 90), 0,
                         ifelse((rawData$AKI_Category_subgroup == "Persistent_Severe_AKI" & rawData$icu_to_aki_2_3_days < 90), 1, 0))
rawData$pAKI_time = ifelse(rawData$pAKI_90==1, rawData$icu_to_aki_2_3_days, 90)


# the composite outcome and time
rawData$dead_pAKI_90 = ifelse((rawData$dead_90==1 | rawData$pAKI_90 == 1), 1, 0)
rawData$dead_pAKI_time = ifelse((rawData$dead_90==1 & rawData$pAKI_90==0), rawData$survival_time,
                                ifelse((rawData$dead_90==0 & rawData$pAKI_90==1), rawData$pAKI_time,
                                       ifelse((rawData$dead_90==1 & rawData$pAKI_90==1 & rawData$survival_time<rawData$pAKI_time), rawData$survival_time, 
                                              ifelse((rawData$dead_90==1 & rawData$pAKI_90==1 & rawData$survival_time>rawData$pAKI_time), rawData$pAKI_time, 90))))



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
  if(analysis == "_intervention_by_each"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  data = data %>%
    mutate(mv = ifelse(MV_DAYS>0, 1, 0))
  
  # Cox proportional hazard model
  for (int in 1:length(int_var)){
    
    coxdata = data %>%
      filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
      mutate(race = ifelse(race==1, 0, 1)) %>%
      select(all_of(adj_var[[int]]), all_of(int_var[int]))

    cox_fit = coxph(Surv(dead_pAKI_time, dead_pAKI_90) ~ ., data=coxdata)
    cox_fit_int_list[[int]] = cox_fit
    
  }
  # append the fitted object to the list
  cox_fit_list[[i]] = cox_fit_int_list
}

# save the list as rds
saveRDS(cox_fit_list, paste("G:/Persistent_AKI/Xinlei/Data/cox_fit_list", analysis, ".rds", sep = ""))


# read the saved results
cox_fit_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/cox_fit_list", analysis, ".rds", sep = ""))



for (j in 1:4){
  cox_fit_df = data.frame()
  for (i in 1:length(cox_fit_list)){
    cox_fit = cox_fit_list[[i]][[j]]
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
  beta_df$Vt = beta_df$Vb + beta_df$Vw
  beta_df$se_pool = sqrt(beta_df$Vt)
  beta_df$z = beta_df$beta_bar/beta_df$se_pool
  beta_df$p_value = 2*pnorm(-abs(beta_df$z))
  
  beta_df$variable_true = adj_var_true[[j]]
  
  beta_df$hazard_ratio = exp(beta_df$beta_bar)
  
  order = adj_var_order[[j]]
  
  cox_res_df = beta_df %>%
    mutate(variable_true = factor(variable_true, levels = order)) %>%
    arrange(variable_true) %>%
    select(variable_true, beta_bar, hazard_ratio, se_pool, z, p_value)
  colnames(cox_res_df) = c("Variable", "Coefficient", "Hazard Ratio", "Std Error", "Z-Statistics", "P-Value")
  write.csv(cox_res_df, paste("G:/Persistent_AKI/Xinlei/Result/cox_res_df", analysis, "_", int_var[j], ".csv", sep = ""))

}
