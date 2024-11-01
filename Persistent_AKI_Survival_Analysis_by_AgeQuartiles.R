library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
#library(interval)
#library(icenReg)


# setting
analysis="_survival_ageq"
# analysis = c("_survival_ageq_sensitivity")

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
contrast_list = c("Persistent_Transient", 
                  "Persistent_No",
                  "Persistent_Stage1",
                  "Persistent_No_and_Stage1",
                  "Persistent_All",
                  "Transient_No",
                  "Transient_Stage1",
                  "Transient_No_and_Stage1")

Q0_list = vector()
Q1_list = vector()
Q2_list = vector()
Q3_list = vector()
Q4_list = vector()
n1_list = vector()
n2_list = vector()
n3_list = vector()
n4_list = vector()

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
  if (analysis == "_survival_ageq"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  # reformat the outcome
  data$AKI_Category = as.character(data$AKI_Category)
  data$AKI_Category_true = ifelse((data$aki_1_icu==1 & data$AKI_Category=="no_AKI"), "stage_1_AKI", data$AKI_Category)
  data$AKI_Category = factor(data$AKI_Category, 
                             levels = c("no_AKI", "Transient_AKI", "Persistent_AKI"))
  data$AKI_Category_true = factor(data$AKI_Category_true, 
                                  levels = c("no_AKI", "stage_1_AKI", "Transient_AKI", "Persistent_AKI")) 
  
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
  subdata$status = ifelse((is.na(subdata$death_date_sh_true) & subdata$dead_90==1), 3, subdata$dead_90)
  subdata$time = ifelse(subdata$status==3, 0, subdata$survival_time)
  subdata$time2 = subdata$survival_time
  
  # KM estimator and curve
  kmdata = subdata %>% select(survival_time, dead_90, status, AKI_Category_true, time, time2)
  if(i==1){saveRDS(kmdata, paste("G:/Persistent_AKI/Xinlei/Data/kmdata", analysis, ".rds", sep = ""))}
  # surve_obj = Surv(time = kmdata$time, time2 = kmdata$time2, event = kmdata$status, type = "interval")
  # km_fit = survfit(surve_obj ~ AKI_Category_true, data = kmdata)
  # ggsurvplot(km_fit, data = kmdata, pval = TRUE)
  # log-rank test for testing the differences between curves
  # logrank_fit = survdiff(surve_obj ~ AKI_Category_true, data=kmdata)
  
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
              #"time", "time2",
              "dead_90",
              "survival_time",
              #"status",
              "AKI_Category_true")
  
  coxdata = subdata[, c(names(subdata) %in% adj_var)] %>%
    filter(AKI_Category_true %in% c("Persistent_AKI", "Transient_AKI")) %>%
    mutate(race = ifelse(race==1, 0, 1),
           AKI_Category_true = factor(AKI_Category_true, levels=c("Transient_AKI", "Persistent_AKI")),
           gender_bmi = (as.numeric(as.character(gender)))*BMI)
  
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
    filter(age_at_admission>=Q3 & age_at_admission<Q4) %>%
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
  
  for (ageq in 1:4) {
    
    coxdata_sub = coxdata %>%
      filter(age_at_admission < age_quartiles[(ageq+1)]) %>%
      filter(age_at_admission >= age_quartiles[(ageq)]) %>%
      select(-c(age_at_admission))
    
    cox_fit = coxph(Surv(survival_time, dead_90) ~ ., data=coxdata_sub)

    # append the fitted object to the list
    cox_fit_ageq_list[[ageq]] = cox_fit
    
  }
  cox_fit_list[[i]] = cox_fit_ageq_list
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

n1_list + n2_list + n3_list + n4_list

# save the list as rds
saveRDS(cox_fit_list, paste("G:/Persistent_AKI/Xinlei/Data/cox_fit_list", analysis, ".rds", sep = ""))


# read the saved results
cox_fit_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/cox_fit_list", analysis, ".rds", sep = ""))


# for the 5 imputed dataset, the KM curve is the same... Since it's only about survival and time...
kmdata = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/kmdata", analysis, ".rds", sep = ""))
#surve_obj = Surv(time = kmdata$time, time2 = kmdata$time2, event = kmdata$status, type = "interval")
surve_obj = Surv(time = kmdata$survival_time, event = kmdata$dead_90)
km_fit = survfit(surve_obj ~ AKI_Category_true, data = kmdata)
km_plot = ggsurvplot(km_fit, data = kmdata,
                     size = 0.1,                # change line size
                     xlim = c(0,90),
                     break.x.by = 10,
                     palette =
                       c("#71771e", "#2E9FDF", "#E7B800", "#cd4b4b"),# custom color palettes
                     conf.int = TRUE,          # Add confidence interval
                     pval = T,              # Add p-value
                     risk.table = TRUE,        # Add risk table
                     risk.table.col = "strata",# Risk table color by groups
                     legend.labs =
                       c("No AKI", "Stage-1 AKI", "Transient AKI", "Persistent AKI"),    # Change legend labels
                     risk.table.height = 0.25, # Useful to change when you have multiple groups
                     ggtheme = theme_bw()           
)
pdf(paste("G:/Persistent_AKI/Xinlei/Result/KM_Plot/km_plot", ".pdf", sep=""), width = 15, height = 8)
print(km_plot$plot, newpage = FALSE)
dev.off()
#logrank_fit = interval::ictest(surve_obj ~ AKI_Category_true, score="logrank1", data=kmdata)
logrank_fit = survdiff(surve_obj ~ AKI_Category_true, data=kmdata)


## pair-wise comparison of the KM curve between AKI categories.
pair_list = list(c("Persistent_AKI", "Transient_AKI"),
                 c("Persistent_AKI", "stage_1_AKI"),
                 c("Persistent_AKI", "no_AKI"),
                 c("Transient_AKI", "stage_1_AKI"),
                 c("Transient_AKI", "no_AKI"),
                 c("stage_1_AKI", "no_AKI"))

labs_list = list(c("Transient AKI", "Persistent AKI"),
                 c("Stage-1 AKI", "Persistent AKI"),
                 c("No AKI", "Persistent AKI"),
                 c("Stage-1 AKI", "Transient AKI"),
                 c("No AKI", "Transient AKI"),
                 c("No AKI", "Stage-1 AKI"))

color_list = list(c("#E7B800", "#cd4b4b"),
                  c("#2E9FDF", "#cd4b4b"),
                  c("#71771e", "#cd4b4b"),
                  c("#2E9FDF", "#E7B800"),
                  c("#71771e", "#E7B800"),
                  c("#71771e", "#2E9FDF"))


km_plot_list = list()
log_rank_list = list()

for (i in 1:length(pair_list)){
  pair = pair_list[[i]]
  kmdata_sub = kmdata %>% filter(AKI_Category_true %in% pair)
  surve_obj = Surv(kmdata_sub$survival_time, kmdata_sub$dead_90)
  km_fit = survfit(surve_obj ~ AKI_Category_true, data = kmdata_sub)
  km_plot = ggsurvplot(km_fit, data = kmdata_sub,
                       size = 0.1,                # change line size
                       xlim = c(0,90),
                       break.x.by = 10,
                       palette = color_list[[i]],# custom color palettes
                       conf.int = TRUE,          # Add confidence interval
                       pval = TRUE,              # Add p-value
                       risk.table = TRUE,        # Add risk table
                       risk.table.col = "strata",# Risk table color by groups
                       legend.labs = labs_list[[i]],# Change legend labels
                       risk.table.height = 0.25, # Useful to change when you have multiple groups
                       ggtheme = theme_bw())
  pdf(paste("G:/Persistent_AKI/Xinlei/Result/KM_Plot/km_plot_", pair[1], "_", pair[2], ".pdf", sep=""), width = 15, height = 8)
  print(km_plot$plot, newpage = FALSE)
  dev.off()
  #logrank_fit = interval::ictest(surve_obj ~ AKI_Category_true, score="logrank1", data=kmdata_sub)
  logrank_fit = survdiff(surve_obj ~ AKI_Category_true, data=kmdata_sub)
  
  km_plot_list[[i]] = km_plot
  log_rank_list[[i]] = logrank_fit
}

saveRDS(km_plot_list, paste("G:/Persistent_AKI/Xinlei/Data/km_plot_list", analysis, ".rds", sep = ""))
saveRDS(log_rank_list, paste("G:/Persistent_AKI/Xinlei/Data/log_rank_list", analysis, ".rds", sep = ""))




for (j in 1:4){
  
  # combine the results from the Cox ph model
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
  
  beta_df$variable_true = c("Persistent AKI: yes",
                            "Apache 3",
                            "BMI",
                            "Charlson",
                            "Gender: male * BMI",
                            "Gender: male",
                            "Days of mechanical ventilation",
                            "Race: non-white",
                            "Sepsis: yes",
                            "Total fluid balance before AKI")
  
  beta_df$hazard_ratio = exp(beta_df$beta_bar)
  
  order = c("Gender: male",
            "Race: non-white",
            "Apache 3",
            "BMI",
            "Gender: male * BMI",
            "Charlson",
            "Days of mechanical ventilation",
            "Sepsis: yes",
            "Total fluid balance before AKI",
            "Persistent AKI: yes")
  
  cox_res_df = beta_df %>%
    mutate(variable_true = factor(variable_true, levels = order)) %>%
    arrange(variable_true) %>%
    select(variable_true, beta_bar, hazard_ratio, se_pool, z, p_value)
  colnames(cox_res_df) = c("Variable", "Coefficient", "Hazard Ratio", "Std Error", "Z-Statistics", "P-Value")
  write.csv(cox_res_df, paste("G:/Persistent_AKI/Xinlei/Result/cox_res_df", analysis, "_", j, ".csv", sep = ""))
  
}


