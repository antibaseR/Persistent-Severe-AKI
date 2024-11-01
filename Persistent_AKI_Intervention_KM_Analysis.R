library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
#library(interval)
#library(icenReg)


# setting
analysis="_intervention_KM"
# analysis = c("_intervention_KM_sensitivity")

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

# create list to save the results
pvalue_nephrotoxin_list = list()
pvalue_diuretic_list = list()
pvalue_vasopressor_list = list()
pvalue_saline_list = list()

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
             "AKI_Category_subgroup",
             "dead",
             "death_date_sh",
             "encounter_start_date_time_sh",
             "encounter_end_date_time_sh",
             "dt_icu_start_sh",
             "NEPHROTOXIN_INTERVENTION",
             "DIURETIC_INTERVENTION",
             "VASOPRESSOR_INTERVENTION",
             "SALINE_INTERVENTION")



i=1
for (i in 1:1){
  
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
  if (analysis == "_survival"){
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
  
  
  intervention_list = c("NEPHROTOXIN_INTERVENTION",
                        "DIURETIC_INTERVENTION",
                        "VASOPRESSOR_INTERVENTION",
                        "SALINE_INTERVENTION")
  
  
  # KM estimator and curve
  # NEPHROTOXIN_INTERVENTION
  int = "NEPHROTOXIN_INTERVENTION"
  kmdata = subdata %>%
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup=="Persistent_Severe_AKI", 1, 0))) %>%
    select(survival_time, dead_90, AKI_Category_subgroup, all_of(int))
  saveRDS(kmdata, paste("G:/Persistent_AKI/Xinlei/Data/kmdata", analysis, ".rds", sep = ""))
  
  kmdata_sub = kmdata %>% filter(AKI_Category_subgroup == "1")
  surve_obj = Surv(time = kmdata_sub$survival_time, event = kmdata_sub$dead_90)
  logrank_fit_NEPHROTOXIN_INTERVENTION_p = survdiff(surve_obj ~ NEPHROTOXIN_INTERVENTION, data=kmdata_sub)
  kmdata_sub = kmdata %>% filter(AKI_Category_subgroup == "0")
  surve_obj = Surv(time = kmdata_sub$survival_time, event = kmdata_sub$dead_90)
  logrank_fit_NEPHROTOXIN_INTERVENTION_t = survdiff(surve_obj ~ NEPHROTOXIN_INTERVENTION, data=kmdata_sub)
  
  surve_obj = Surv(time = kmdata$survival_time, event = kmdata$dead_90)
  km_fit = survfit(surve_obj ~ AKI_Category_subgroup + NEPHROTOXIN_INTERVENTION, data = kmdata)
  km_plot = ggsurvplot(km_fit, data = kmdata,
                       size = 0.1,                # change line size
                       xlim = c(0,90),
                       break.x.by = 10,
                       palette = c("#E7B800", "#A88600", "#cd4b4b", "#812E2E"),# custom color palettes
                       conf.int = TRUE,          # Add confidence interval
                       pval = F,              # Add p-value
                       risk.table = F,        # Add risk table
                       risk.table.col = "strata",# Risk table color by groups
                       legend.labs = c("No persistent severe AKI without nephrotoxin avoidance",
                                       "No persistent severe AKI with nephrotoxin avoidance",
                                       "Persistent severe AKI without nephrotoxin avoidance",
                                       "Persistent severe AKI with nephrotoxin avoidance"),    # Change legend labels
                       risk.table.height = 0.25,
                       font.x = c(22), 
                       font.y = c(22), 
                       font.tickslab = c(18), 
                       font.legend = c(16),
                       legend = "right")
  
  km_plot$plot = km_plot$plot + 
    ggplot2::annotate("text",
                      size = 8,
                      x = 15, y = 1,
                      label = paste("P-value: ", ifelse(logrank_fit_NEPHROTOXIN_INTERVENTION_t$pvalue<0.001, "< 0.001", round(logrank_fit_NEPHROTOXIN_INTERVENTION_t$pvalue, 3)))
    ) +
    ggplot2::annotate("text",
                      size = 8,
                      x = 10, y = 0.5,
                      label = paste("P-value: ", ifelse(logrank_fit_NEPHROTOXIN_INTERVENTION_p$pvalue<0.001, "< 0.001", round(logrank_fit_NEPHROTOXIN_INTERVENTION_p$pvalue, 3)))
                      )
  
  pdf(paste("G:/Persistent_AKI/Xinlei/Result/KM_Plot/km_plot", analysis, "_", int, ".pdf", sep=""), width = 20, height = 10)
  print(km_plot$plot, newpage = FALSE)
  dev.off()
  
  
  
  
  
  int = "DIURETIC_INTERVENTION"
  kmdata = subdata %>%
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup=="Persistent_Severe_AKI", 1, 0))) %>%
    select(survival_time, dead_90, AKI_Category_subgroup, all_of(int))
  saveRDS(kmdata, paste("G:/Persistent_AKI/Xinlei/Data/kmdata", analysis, ".rds", sep = ""))
  
  kmdata_sub = kmdata %>% filter(AKI_Category_subgroup == "1")
  surve_obj = Surv(time = kmdata_sub$survival_time, event = kmdata_sub$dead_90)
  logrank_fit_DIURETIC_INTERVENTION_p = survdiff(surve_obj ~ DIURETIC_INTERVENTION, data=kmdata_sub)
  kmdata_sub = kmdata %>% filter(AKI_Category_subgroup == "0")
  surve_obj = Surv(time = kmdata_sub$survival_time, event = kmdata_sub$dead_90)
  logrank_fit_DIURETIC_INTERVENTION_t = survdiff(surve_obj ~ DIURETIC_INTERVENTION, data=kmdata_sub)
  
  surve_obj = Surv(time = kmdata$survival_time, event = kmdata$dead_90)
  km_fit = survfit(surve_obj ~ AKI_Category_subgroup + DIURETIC_INTERVENTION, data = kmdata)
  km_plot = ggsurvplot(km_fit, data = kmdata,
                       size = 0.1,                # change line size
                       xlim = c(0,90),
                       break.x.by = 10,
                       palette = c("#E7B800", "#A88600", "#cd4b4b", "#812E2E"),
                       conf.int = TRUE,          # Add confidence interval
                       pval = F,              # Add p-value
                       risk.table = F,        # Add risk table
                       risk.table.col = "strata",# Risk table color by groups
                       legend.labs = c("No persistent severe AKI without diuretic avoidance",
                                       "No persistent severe AKI with diuretic avoidance",
                                       "Persistent severe AKI without diuretic avoidance",
                                       "Persistent severe AKI with diuretic avoidance"),    # Change legend labels
                       risk.table.height = 0.25,
                       font.x = c(22), 
                       font.y = c(22), 
                       font.tickslab = c(18), 
                       font.legend = c(16),
                       legend = "right")
  
  km_plot$plot = km_plot$plot + 
    ggplot2::annotate("text",
                      size = 8,
                      x = 15, y = 1,
                      label = paste("P-value: ", ifelse(logrank_fit_DIURETIC_INTERVENTION_t$pvalue<0.001, "< 0.001", round(logrank_fit_DIURETIC_INTERVENTION_t$pvalue, 3)))
    ) +
    ggplot2::annotate("text",
                      size = 8,
                      x = 10, y = 0.5,
                      label = paste("P-value: ", ifelse(logrank_fit_DIURETIC_INTERVENTION_p$pvalue<0.001, "< 0.001", round(logrank_fit_DIURETIC_INTERVENTION_p$pvalue, 3)))
                      )
  
  pdf(paste("G:/Persistent_AKI/Xinlei/Result/KM_Plot/km_plot", analysis, "_", int, ".pdf", sep=""), width = 20, height = 10)
  print(km_plot$plot, newpage = FALSE)
  dev.off()
  
  
  
  
  
  
  int = "VASOPRESSOR_INTERVENTION"
  kmdata = subdata %>%
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup=="Persistent_Severe_AKI", 1, 0))) %>%
    select(survival_time, dead_90, AKI_Category_subgroup, all_of(int))
  saveRDS(kmdata, paste("G:/Persistent_AKI/Xinlei/Data/kmdata", analysis, ".rds", sep = ""))
  
  kmdata_sub = kmdata %>% filter(AKI_Category_subgroup == "1")
  surve_obj = Surv(time = kmdata_sub$survival_time, event = kmdata_sub$dead_90)
  logrank_fit_VASOPRESSOR_INTERVENTION_p = survdiff(surve_obj ~ VASOPRESSOR_INTERVENTION, data=kmdata_sub)
  kmdata_sub = kmdata %>% filter(AKI_Category_subgroup == "0")
  surve_obj = Surv(time = kmdata_sub$survival_time, event = kmdata_sub$dead_90)
  logrank_fit_VASOPRESSOR_INTERVENTION_t = survdiff(surve_obj ~ VASOPRESSOR_INTERVENTION, data=kmdata_sub)
  
  surve_obj = Surv(time = kmdata$survival_time, event = kmdata$dead_90)
  km_fit = survfit(surve_obj ~ AKI_Category_subgroup + VASOPRESSOR_INTERVENTION, data = kmdata)
  km_plot = ggsurvplot(km_fit, data = kmdata,
                       size = 0.1,                # change line size
                       xlim = c(0,90),
                       break.x.by = 10,
                       palette = c("#E7B800", "#A88600", "#cd4b4b", "#812E2E"),
                       conf.int = TRUE,          # Add confidence interval
                       pval = F,              # Add p-value
                       risk.table = F,        # Add risk table
                       risk.table.col = "strata",# Risk table color by groups
                       legend.labs = c("No persistent severe AKI with late vasopressor/ineotrope intervention",
                                       "No persistent severe AKI with early vasopressor/ineotrope intervention",
                                       "Persistent severe AKI with late vasopressor/ineotrope intervention",
                                       "Persistent severe AKI with early vasopressor/ineotrope intervention"),    # Change legend labels
                       risk.table.height = 0.25,
                       font.x = c(22), 
                       font.y = c(22), 
                       font.tickslab = c(18), 
                       font.legend = c(16),
                       legend = "right")
  
  km_plot$plot = km_plot$plot + 
    ggplot2::annotate("text",
                      size = 8,
                      x = 15, y = 1,
                      label = paste("P-value: ", ifelse(logrank_fit_VASOPRESSOR_INTERVENTION_t$pvalue<0.001, "< 0.001", round(logrank_fit_VASOPRESSOR_INTERVENTION_t$pvalue, 3)))
    ) +
    ggplot2::annotate("text",
                      size = 8,
                      x = 10, y = 0.35,
                      label = paste("P-value: ", ifelse(logrank_fit_VASOPRESSOR_INTERVENTION_p$pvalue<0.001, "< 0.001", round(logrank_fit_VASOPRESSOR_INTERVENTION_p$pvalue, 3)))
                      )
  
  pdf(paste("G:/Persistent_AKI/Xinlei/Result/KM_Plot/km_plot", analysis, "_", int, ".pdf", sep=""), width = 20, height = 10)
  print(km_plot$plot, newpage = FALSE)
  dev.off()
  
  
  
  
  
  
  int = "SALINE_INTERVENTION"
  kmdata = subdata %>%
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup=="Persistent_Severe_AKI", 1, 0))) %>%
    select(survival_time, dead_90, AKI_Category_subgroup, all_of(int))
  saveRDS(kmdata, paste("G:/Persistent_AKI/Xinlei/Data/kmdata", analysis, ".rds", sep = ""))
  
  # log-rank test for the p-value
  kmdata_sub = kmdata %>% filter(AKI_Category_subgroup == "1")
  surve_obj = Surv(time = kmdata_sub$survival_time, event = kmdata_sub$dead_90)
  logrank_fit_SALINE_INTERVENTION_p = survdiff(surve_obj ~ SALINE_INTERVENTION, data=kmdata_sub)
  kmdata_sub = kmdata %>% filter(AKI_Category_subgroup == "0")
  surve_obj = Surv(time = kmdata_sub$survival_time, event = kmdata_sub$dead_90)
  logrank_fit_SALINE_INTERVENTION_t = survdiff(surve_obj ~ SALINE_INTERVENTION, data=kmdata_sub)
  
  surve_obj = Surv(time = kmdata$survival_time, event = kmdata$dead_90)
  km_fit = survfit(surve_obj ~ AKI_Category_subgroup + SALINE_INTERVENTION, data = kmdata)
  km_plot = ggsurvplot(km_fit, data = kmdata,
                       size = 0.1,                # change line size
                       xlim = c(0,90),
                       break.x.by = 10,
                       palette = c("#E7B800", "#A88600", "#cd4b4b", "#812E2E"),
                       conf.int = TRUE,          # Add confidence interval
                       pval = F,              # Add p-value
                       risk.table = F,        # Add risk table
                       risk.table.col = "strata",# Risk table color by groups
                       legend.labs = c("No persistent severe AKI without saline avoidance",
                                       "No persistent severe AKI with saline avoidance",
                                       "Persistent severe AKI without saline avoidance",
                                       "Persistent severe AKI with saline avoidance"),    # Change legend labels
                       risk.table.height = 0.25,
                       font.x = c(22), 
                       font.y = c(22), 
                       font.tickslab = c(18), 
                       font.legend = c(16),
                       legend = "right")
  
  km_plot$plot = km_plot$plot + 
    ggplot2::annotate("text",
                      size = 8,
                      x = 15, y = 1,
                      label = paste("P-value: ", ifelse(logrank_fit_SALINE_INTERVENTION_t$pvalue<0.001, "< 0.001", round(logrank_fit_SALINE_INTERVENTION_t$pvalue, 3)))
                      ) +
    ggplot2::annotate("text",
                      size = 8,
                      x = 10, y = 0.4,
                      label = paste("P-value: ", ifelse(logrank_fit_SALINE_INTERVENTION_p$pvalue<0.001, "< 0.001", round(logrank_fit_SALINE_INTERVENTION_p$pvalue, 3)))
                      )
  
  pdf(paste("G:/Persistent_AKI/Xinlei/Result/KM_Plot/km_plot", analysis, "_", int, ".pdf", sep=""), width = 20, height = 10)
  print(km_plot$plot, newpage = FALSE)
  dev.off()
  
  
  
  
  pvalue_nephrotoxin_list[[i]] = c(logrank_fit_NEPHROTOXIN_INTERVENTION_p$pvalue, logrank_fit_NEPHROTOXIN_INTERVENTION_t$pvalue)
  pvalue_diuretic_list[[i]] = c(logrank_fit_DIURETIC_INTERVENTION_p$pvalue, logrank_fit_DIURETIC_INTERVENTION_t$pvalue)
  pvalue_vasopressor_list[[i]] = c(logrank_fit_VASOPRESSOR_INTERVENTION_p$pvalue, logrank_fit_VASOPRESSOR_INTERVENTION_t$pvalue)
  pvalue_saline_list[[i]] = c(logrank_fit_SALINE_INTERVENTION_p$pvalue, logrank_fit_SALINE_INTERVENTION_t$pvalue)

}