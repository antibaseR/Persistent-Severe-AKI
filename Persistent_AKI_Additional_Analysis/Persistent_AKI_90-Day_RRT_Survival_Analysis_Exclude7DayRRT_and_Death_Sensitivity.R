library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
library(randomForestSRC)



for (analysis in c("_RRT_excluding_7days_death_and_rrt")){
  
  impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
  rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
  rrtData = read.csv('G:/Persistent_AKI/Persistent_AKI_Data/RRT_Data.csv')
  rrtData = rrtData %>%
    select(patientid, patientvisitid, first_rrt_date, rrt_90)
  
  # create list to save the model fit
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
               "CKD",
               "refCreat",
               "rrt_90",
               "rrt_days")
  
  
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
    data = merge(merge(data, imp_sub, by = c("patientid", "patientvisitid")), rrtData, by = c("patientid", "patientvisitid"))
    
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
    
    data$first_rrt_date = as.Date(data$first_rrt_date, "%m/%d/%Y")
    data$aki_time = as.Date(data$FIRST_AKI_DATETIME)
    data$rrt_days = difftime(data$first_rrt_date, data$aki_time, units = "days")
    data$rrt_days = ifelse(is.na(data$rrt_days), 90, data$rrt_days)
    data = data %>%
      filter(rrt_days <= 90 & rrt_days >= 0)
    
    subdata = data[, c(names(data) %in% surv_var)] %>%
      mutate(TOTAL_BALANCE_BEFORE_AKI = TOTAL_INPUTS_BEFORE_AKI - TOTAL_OUTPUTS_BEFORE_AKI) %>%
      select(-c(TOTAL_INPUTS_BEFORE_AKI,TOTAL_OUTPUTS_BEFORE_AKI))
    
    ## if the death date is earlier than the ICU admission date than replace it with NA
    subdata$death_date_sh_true = ifelse(difftime(subdata$death_date_sh, subdata$dt_icu_start_sh, units="days")<0, NA, subdata$death_date_sh)
    ## flag for death
    subdata$dead_true = ifelse(subdata$dead==1, 1, ifelse((subdata$dead==0 & is.na(subdata$death_date_sh_true)), 0, 1))
    ## if a subject has death date as NA, replace the NA with the discharge date
    subdata$death_date_sh_true = ifelse(is.na(subdata$death_date_sh_true), subdata$encounter_end_date_time_sh, subdata$death_date_sh_true)
    ## death observation time
    subdata$dtime = difftime(subdata$death_date_sh_true, subdata$FIRST_AKI_DATETIME, units="days")
    subdata$dtime = ifelse(subdata$dtime<0, 0, subdata$dtime)
    ## 90-day mortality
    subdata$dead_90 = ifelse((subdata$dead_true==1 & subdata$dtime <= 90), 1, 0)
    ## reformat dtime
    subdata$dtime = ifelse(subdata$dtime > 90, 90, subdata$dtime)
    
    # excluding paitents who died or had RRT in the first 7 days
    subdata = subdata %>%
      filter(!(dead_90==1 & dtime<=7)) %>%
      filter(rrt_days >= 7)
    
    # Cox proportional hazard model
    ## variables we need to adjust for the adjusted Cox proportional hazard model
    adj_var = c("age_at_admission",
                "gender",
                "race",
                "charlson",
                "apache3",
                "MV_DAYS",
                "refCreat",
                "rrt_90",
                "rrt_days",
                "AKI_Category_subgroup")
    
    coxdata = subdata[, c(names(subdata) %in% adj_var)] %>%
      filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
      mutate(race = ifelse(race==1, 0, 1),
             AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No"), 
                                            levels=c("No", "Yes")),
             apache_refCreat = apache3 * refCreat)
    
    #set.seed(960725)
    #obj = rfsrc(Surv(rrt_days, rrt_90) ~ ., coxdata,
    #ntree = 1, nodesize = 5, nsplit = 10, importance = TRUE)
    #tree = get.tree(obj, 1)
    #tree
    
    cox_fit = coxph(Surv(rrt_days, rrt_90) ~ ., data=coxdata)
    # append the fitted object to the list
    cox_fit_list[[i]] = cox_fit
  }
  
  # save the list as rds
  saveRDS(cox_fit_list, paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/cox_fit_list", analysis, ".rds", sep = ""))
  
  
  # read the saved results
  cox_fit_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/cox_fit_list", analysis, ".rds", sep = ""))
  
  # combine the results from the Cox ph model
  cox_fit_df = data.frame()
  cov_gender_bmi = c()
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
                            "Persistent severe AKI: yes",
                            "Apache 3 * Reference creatinine",
                            "Apache 3",
                            "Charlson",
                            "Gender: male",
                            "Days of mechanical ventilation",
                            "Race: non-white",
                            "Reference creatinine")
  
  beta_df$hazard_ratio = exp(beta_df$beta_bar)
  
  order = c("Age at admission",
            "Gender: male",
            "Race: non-white",
            "Apache 3",
            "Charlson",
            "Days of mechanical ventilation",
            "Reference creatinine",
            "Apache 3 * Reference creatinine",
            "Persistent severe AKI: yes"
  )
  
  cox_res_df = beta_df %>%
    mutate(variable_true = factor(variable_true, levels = order)) %>%
    arrange(variable_true) %>%
    select(variable_true, beta_bar, hazard_ratio, se_pool, z, p_value)
  colnames(cox_res_df) = c("Variable", "Coefficient", "Hazard Ratio", "Std Error", "Z-Statistics", "P-Value")
  
  
  write.csv(cox_res_df, paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/cox_res_df", analysis, ".csv", sep = ""))
  
  rm(list = ls())
  
  
}

