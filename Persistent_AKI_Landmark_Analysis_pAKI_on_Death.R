library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
library(Landmarking)
options(scipen=999)

# setting
analysis="_landmark_pAKI_on_death"
# analysis="_landmark_pAKI_on_death_sensitivity"


# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
landmarkData = read.csv('G:/Persistent_AKI/Persistent_AKI_Data/Landmark_Data.csv')


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
surv_var = c("patientid",
             "patientvisitid",
             "age_at_admission",
             "gender",
             "charlson",
             "apache3",
             "SEPSIS",
             "BMI",
             "lactate",
             "WBC_Max",
             "Bili_max",
             "creat_max",
             "MAX_MAP",
             "dead",
             "death_date_sh",
             "encounter_start_date_time_sh",
             "encounter_end_date_time_sh",
             "dt_icu_start_sh",
             "FIRST_AKI_DATETIME",
             "AKI_Category",
             "AKI_Category_true",
             "AKI_Category_subgroup",
             "aki_2_3_icu",
             "aki_2_3_icu_dt_sh")


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
  if (analysis == "_landmark_pAKI_on_death"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  # subset the data
  subdata = data[, c(names(data) %in% surv_var)]
  
  ## if the death date is earlier than the ICU admission date than replace it with NA
  subdata$death_date_sh_true = ifelse(difftime(subdata$death_date_sh, subdata$dt_icu_start_sh, units="days")<0, NA, subdata$death_date_sh)
  ## if the survived subjects' death date is earlier than or the same as the hospital discharge date, then the subjects died in the hospital
  subdata$discharge_to_death_days = difftime(subdata$death_date_sh_true, subdata$encounter_end_date_time_sh, unit="days")
  subdata$dead_hosp = ifelse((!is.na(subdata$discharge_to_death_days) & subdata$discharge_to_death_days<=0 & subdata$dead==0), 1, as.numeric(as.character(subdata$dead)))
  ## if the dead subjects' death date is later than the hospital discharge date, the subjects survived in hospital
  subdata$dead_hosp = ifelse((!is.na(subdata$discharge_to_death_days) & subdata$discharge_to_death_days>0 & subdata$dead==1), 0, as.numeric(as.character(subdata$dead)))
  
  ## number of days from hospital admission to death
  subdata$admission_to_death_days = difftime(subdata$death_date_sh_true, subdata$encounter_start_date_time_sh, units="days")
  ## number of days from hospital admission to hospital discharge
  subdata$admission_to_discharge_days = difftime(subdata$encounter_end_date_time_sh, subdata$encounter_start_date_time_sh, units="days")
  ## replace missing date with discharge date
  subdata$survived_days = ifelse((is.na(subdata$death_date_sh_true)), subdata$admission_to_discharge_days, subdata$admission_to_death_days)
  
  ## 90-day mortality: 
  ### - if the death date is not NA and the icu to death days is >= 90 -- survived
  ### - if the death date is not NA and the icu to death days is < 90 -- dead
  ### - if the death date is NA that means we didn't observe the event time but we know the result, the data is censored, so the survival = dead_hosp and time = icu_to_discharge_days
  subdata$dead_90 = ifelse(!(is.na(subdata$death_date_sh_true)) & subdata$survived_days>=90, 0, 
                           ifelse(!(is.na(subdata$death_date_sh_true)) & subdata$survived_days<90, 1, subdata$dead_hosp))
  subdata$survival_time = ifelse(subdata$survived_days>90, 90, subdata$survived_days)
  
  ## number of days from ICU admission to the development of pAKI
  subdata$icu_to_AKI_days = difftime(subdata$aki_2_3_icu_dt_sh, subdata$dt_icu_start_sh, unit="days")
  ## NA means no pAKI time, means no pAKI, means icu to pAKI days > 90
  subdata$icu_to_AKI_days = ifelse((is.na(subdata$icu_to_AKI_days)), 90, subdata$icu_to_AKI_days)
  check = subdata %>%
    select(dt_icu_start_sh, aki_2_3_icu, aki_2_3_icu_dt_sh, icu_to_AKI_days)
  
  ## 90-day mortality: 
  ### - if the death date is not NA and the icu to death days is >= 90 -- survived
  ### - if the death date is not NA and the icu to death days is < 90 -- dead
  ### - if the death date is NA that means we didn't observe the event time but we know the result, the data is censored, so the survival = dead_hosp and time = icu_to_discharge_days
  subdata$dead_90 = ifelse(!(is.na(subdata$death_date_sh_true)) & subdata$survived_days>=90, 0, 
                           ifelse(!(is.na(subdata$death_date_sh_true)) & subdata$survived_days<90, 1, subdata$dead_hosp))
  subdata$survival_time = ifelse(subdata$survived_days>90, 90, subdata$survived_days)
  
  # subset the landmark data and merge with the subdata
  landmarkData_sub = landmarkData %>%
    select(-c(encounter_start_date_time_sh, encounter_end_date_time_sh))
  subdata_long = merge(subdata, landmarkData_sub, by=c("patientid", "patientvisitid"))
  
  # days from icu admission date to measurement date
  #subdata_long = subdata_long %>%
  #  mutate(relative_day_icu = difftime(Day_sh, dt_icu_start_sh, unit = "days"))
  
  # Landmark analysis: LOCF model
  ## variables we need to adjust for the model
  adj_var = c("patientid",
              "patientvisitid",
              "age_at_admission",
              "gender",
              "charlson",
              "apache3",
              "SEPSIS",
              "lactate",
              "WBC_Max",
              "Bili_max",
              "creat_max",
              "MAX_MAP",
              "dead_90",
              "survival_time",
              "relative_day",
              "AKI_Category_subgroup")
  
  coxdata = subdata_long[, c(names(subdata_long) %in% adj_var)] %>%
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup=="Persistent_Severe_AKI", 1, 0)))
  
  last_relative_day = coxdata %>%
    group_by(patientid, patientvisitid) %>%
    summarize(last_relative_day = max(relative_day))
  
  coxdata = merge(coxdata, last_relative_day, c("patientid", "patientvisitid"))
  coxdata$last_relative_day = ifelse(coxdata$dead_90==0, 999, coxdata$last_relative_day)
  
  cox_fit = list()
  for (t in 1:14){
    coxdata_t = coxdata %>%
      filter(last_relative_day>=t) %>%
      select(-c("patientid", "patientvisitid", "relative_day", "last_relative_day" ))
    cox_fit[[t]] = coxph(Surv(survival_time, dead_90) ~ ., data=coxdata_t)
  }
  
  # append the fitted object to the list
  cox_fit_list[[i]] = cox_fit
}

# save the list as rds
saveRDS(cox_fit_list, paste("G:/Persistent_AKI/Xinlei/Data/cox_fit_list", analysis, ".rds", sep = ""))


# read the list
cox_fit_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/cox_fit_list", analysis, ".rds", sep = ""))


# combine the odds ratio table
landmark_data = tibble()
for (i in 1:length(cox_fit_list)){
  landmark_data_new = tibble()
  cox_fit = cox_fit_list[[i]]
  for (t in seq(1, 14, 1)){
    new_coef = cox_fit[[t]]$coefficients
    new_se = sqrt(diag(cox_fit[[t]]$var))
    variable = names(new_coef)
    new_data = tibble(variable = variable,
                      beta = new_coef,
                      se = new_se,
                      time = t)
    landmark_data_new = rbind(landmark_data_new, new_data)
  }
  landmark_data_new$imp = i
  landmark_data = rbind(landmark_data, landmark_data_new)
}


beta_df = landmark_data %>%
  group_by(variable, time) %>%
  summarise(beta_bar = mean(beta))

landmark_data = merge(landmark_data, beta_df, by=c("variable", "time"))

Vw_df = landmark_data %>% 
  group_by(variable, time) %>%
  summarise(Vw = mean(se^2))

Vb_df = landmark_data %>%
  group_by(variable, time) %>%
  summarise(Vb = sum((beta-beta_bar)^2)/(length(unique(landmark_data$imp))-1))

beta_df = merge(merge(beta_df, Vw_df, by=c("variable", "time")), Vb_df, by=c("variable", "time"))
beta_df$Vt = beta_df$Vb + beta_df$Vw + (beta_df$Vb/length(unique(landmark_data$imp)))
beta_df$se_pool = sqrt(beta_df$Vt)
beta_df$z = beta_df$beta_bar/beta_df$se_pool
beta_df$p_value = 2*pnorm(-abs(beta_df$z))

landmark_death_odds_ratio = beta_df %>%
  filter(variable == "AKI_Category_subgroup1") %>%
  select(time, beta_bar, se_pool, p_value) %>%
  mutate(`Hazard Ratio` = exp(beta_bar),
         ci_lower = exp(beta_bar - 1.96*se_pool),
         ci_upper = exp(beta_bar + 1.96*se_pool)
  ) %>%
  arrange(time)

landmark_death_odds_ratio %>%
  filter(time > 2) %>%
  ggplot(aes(x = time, y = `Hazard Ratio`)) +
  geom_point(color = "black", aes(size = `Hazard Ratio`)) +
  geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper),
                width = .2, color = "black",
                linetype = "twodash",
                position = position_dodge(0.05)) +
  geom_hline(yintercept = 1, linetype="dashed", color = "grey") + 
  theme_classic() +
  labs(x = "Landmark Time from Risk Exposure", y = "Hazard Ratio (HR) of Death") +
  theme(axis.title.x = element_text(size=14, face="bold"),
        axis.title.y = element_text(size=14, face="bold")) +
  scale_x_continuous(breaks = seq(3,14,1), labels = seq(1, 12, 1)) +
  scale_y_continuous(breaks = seq(1,1.5,0.1), limits = c(0.9, 1.56))

ggsave("G:/Persistent_AKI/Xinlei/Result/plot_landmark.png", width=10, height=5)




