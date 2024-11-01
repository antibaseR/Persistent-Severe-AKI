library(tidyverse)
library(survminer)
library(survival)
library(lubridate)

analysis = "_survival_analysis"
#analysis = "_survival_sensitivity"

for (analysis in c("_survival_analysis", "_survival_sensitivity")){
  
  impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
  rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
  
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
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
    
    # subset the data
    if (analysis=="_survival_sensitivity"){
      subdata = data %>%
        filter(!(difftime(dt_icu_end_sh, dt_icu_start_sh, units="days")<=2 & dead==1))
    } else {subdata = data}
    
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
                "gender",
                "race",
                "charlson",
                "apache3",
                "SEPSIS",
                "MV_DAYS",
                "BMI",
                "TOTAL_BALANCE_BEFORE_AKI",
                "dead_90",
                "survival_time",
                "AKI_Category_subgroup")
    
    coxdata = subdata[, c(names(subdata) %in% adj_var)] %>%
      filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
      mutate(race = ifelse(race==1, 0, 1),
             AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No"), 
                                            levels=c("No", "Yes")),
             apache_sepsis = apache3 * as.numeric(SEPSIS))
    
    #set.seed(960725)
    #obj = rfsrc(Surv(survival_time, dead_90) ~ ., coxdata,
    #ntree = 1, nodesize = 5, nsplit = 10, importance = TRUE)
    #tree = get.tree(obj, 1)
    #tree
    
    cox_fit = coxph(Surv(survival_time, dead_90) ~ ., data=coxdata)
    
    if (i==1){
      
      # model diagnosis
      ## Martingale residuals
      coxdata$resid_mart = residuals(cox_fit, type = "martingale")
      ## Cox-Snell residuals
      coxdata$resid_coxsnell <- -(coxdata$resid_mart - coxdata$dead_90)
      ## Fit model on Cox-Snell residuals (Approximately Expo(1) distributed under correct model)
      fit_coxsnell <- coxph(formula = Surv(resid_coxsnell, dead_90) ~ 1,
                            data    = coxdata,
                            ties    = c("efron","breslow","exact")[1])
      ## Nelson-Aalen estimator for baseline hazard (all covariates zero)
      df_base_haz <- basehaz(fit_coxsnell, centered = FALSE)
      ## Cox ph assumption
      test.ph = cox.zph(cox_fit)
      
      # Plot
      ## Martingale plot
      plt1 = ggplot(data = coxdata, mapping = aes(x = MV_DAYS, y = resid_mart)) +
        geom_point() +
        geom_smooth() +
        labs(title = "Days of mechanical ventilation",
             x = "Days of mechanical ventilation",
             y = "Martingale residual") +
        theme_bw() + theme(legend.key = element_blank())
      
      plt2 = ggplot(data = coxdata, mapping = aes(x = age_at_admission, y = resid_mart)) +
        geom_point() +
        geom_smooth() +
        labs(title = "Age",
             x = "Age",
             y = "Martingale residual") +
        theme_bw() + theme(legend.key = element_blank())
      
      plt3 = ggplot(data = coxdata, mapping = aes(x = charlson, y = resid_mart)) +
        geom_point() +
        geom_smooth() +
        labs(title = "Charlson",
             x = "Charlson",
             y = "Martingale residual") +
        theme_bw() + theme(legend.key = element_blank())
      
      plt4 = ggplot(data = coxdata, mapping = aes(x = apache3, y = resid_mart)) +
        geom_point() +
        geom_smooth() +
        labs(title = "APACHE 3",
             x = "APACHE 3",
             y = "Martingale residual") +
        theme_bw() + theme(legend.key = element_blank())
      
      plt5 = ggplot(data = coxdata, mapping = aes(x = BMI, y = resid_mart)) +
        geom_point() +
        geom_smooth() +
        labs(title = "BMI",
             x = "BMI",
             y = "Martingale residual") +
        theme_bw() + theme(legend.key = element_blank())
      
      plt6 = ggplot(data = coxdata, mapping = aes(x = TOTAL_BALANCE_BEFORE_AKI, y = resid_mart)) +
        geom_point() +
        geom_smooth() +
        labs(title = "Total fluid balance before AKI",
             x = "Total fluid balance before AKI",
             y = "Martingale residual") +
        theme_bw() + theme(legend.key = element_blank())
      
      (plt1 + plt2 + plt3)/(plt4 + plt5 + plt6)
      
      ggsave(paste("G:/Persistent_AKI/Xinlei/Result/Diagnosis/martingale_plot", analysis, ".png", sep = ""), width = 10, height = 8)
      
      ## Overall fit
      ggplot(data = df_base_haz, mapping = aes(x = time, y = hazard)) +
        geom_point() +
        scale_x_continuous(limit = c(0,3.5)) +
        scale_y_continuous(limit = c(0,3.5)) +
        labs(x = "Cox-Snell residuals as pseudo observed times",
             y = "Estimated cumulative hazard at pseudo observed times") +
        theme_bw() + theme(legend.key = element_blank())
      
      ### C-index
      cindex = survConcordance(Surv(coxdata$survival_time, coxdata$dead_90) ~ predict(cox_fit))
      print(cindex) # 0.73
      
      ggsave(paste("G:/Persistent_AKI/Xinlei/Result/Diagnosis/cox_snell_plot", analysis, ".png", sep = ""), width = 8, height = 8)
      
      ## Ph assumption
      ph.plt = ggcoxzph(test.ph)
      
      pb1 = ggplot_build(ph.plt$`1`)
      pb1$plot$labels$y = "Beta(t) for days of mechanical ventilation"
      pb1 = ggplot_gtable(pb1)
      
      pb2 = ggplot_build(ph.plt$`2`)
      pb2$plot$labels$y = "Beta(t) for sepsis"
      pb2 = ggplot_gtable(pb2)
      
      pb3 = ggplot_build(ph.plt$`3`)
      pb3$plot$labels$y = "Beta(t) for persistent severe AKI"
      pb3 = ggplot_gtable(pb3)
      
      pb4 = ggplot_build(ph.plt$`4`)
      pb4$plot$labels$y = "Beta(t) for age"
      pb4 = ggplot_gtable(pb4)
      
      pb5 = ggplot_build(ph.plt$`5`)
      pb5$plot$labels$y = "Beta(t) for gender"
      pb5 = ggplot_gtable(pb5)
      
      pb6 = ggplot_build(ph.plt$`6`)
      pb6$plot$labels$y = "Beta(t) for race"
      pb6 = ggplot_gtable(pb6)
      
      pb7 = ggplot_build(ph.plt$`7`)
      pb7$plot$labels$y = "Beta(t) for charlson"
      pb7 = ggplot_gtable(pb7)
      
      pb8 = ggplot_build(ph.plt$`8`)
      pb8$plot$labels$y = "Beta(t) for Apache 3"
      pb8 = ggplot_gtable(pb8)
      
      pb9 = ggplot_build(ph.plt$`9`)
      pb9$plot$labels$y = "Beta(t) for BMI"
      pb9 = ggplot_gtable(pb9)
      
      pb10 = ggplot_build(ph.plt$`10`)
      pb10$plot$labels$y = "Beta(t) for total fluid balance before AKI"
      pb10 = ggplot_gtable(pb10)
      
      pb11 = ggplot_build(ph.plt$`11`)
      pb11$plot$labels$y = "Beta(t) for apache 3 * sepsis"
      pb11 = ggplot_gtable(pb11)
      
      pb.plt1 = ggplotify::as.ggplot(pb1)
      pb.plt2 = ggplotify::as.ggplot(pb2)
      pb.plt3 = ggplotify::as.ggplot(pb3)
      pb.plt4 = ggplotify::as.ggplot(pb4)
      pb.plt5 = ggplotify::as.ggplot(pb5)
      pb.plt6 = ggplotify::as.ggplot(pb6)
      pb.plt7 = ggplotify::as.ggplot(pb7)
      pb.plt8 = ggplotify::as.ggplot(pb8)
      pb.plt9 = ggplotify::as.ggplot(pb9)
      pb.plt10 = ggplotify::as.ggplot(pb10)
      pb.plt11 = ggplotify::as.ggplot(pb11)
      
      (pb.plt1 + pb.plt2 + pb.plt3 + pb.plt4 + pb.plt5 + pb.plt6 + pb.plt7 + pb.plt8 + pb.plt9 + pb.plt10 + pb.plt11) + 
        plot_layout(ncol = 3) +
        plot_annotation(
          title = 'Global Schoenfeld Test p < 0.001',
          theme=theme(plot.title=element_text(hjust=0.5, size = 20))
        )
      
      ggsave(paste("G:/Persistent_AKI/Xinlei/Result/Diagnosis/ph_plot", analysis, ".png", sep = ""), width = 15, height = 15)
      
      ## Testing influential observations
      ggcoxdiagnostics(cox_fit, type = "deviance", sline.se = F,
                       linear.predictions = FALSE, ggtheme = theme_bw())
      
      ggsave(paste("G:/Persistent_AKI/Xinlei/Result/Diagnosis/deviance_plot", analysis, ".png", sep = ""), width = 8, height = 8)
      
    }
    
    
    
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
  beta_df$Vt = beta_df$Vb + beta_df$Vw + + (beta_df$Vb/length(unique(cox_fit_df$imp)))
  beta_df$se_pool = sqrt(beta_df$Vt)
  beta_df$z = beta_df$beta_bar/beta_df$se_pool
  beta_df$p_value = 2*pnorm(-abs(beta_df$z))
  
  beta_df$variable_true = c("Age at admission",
                            "Persistent severe AKI: yes",
                            "Apache 3 * Sepsis: yes",
                            "Apache 3",
                            "BMI",
                            "Charlson",
                            "Gender: male",
                            "Days of mechanical ventilation",
                            "Race: non-white",
                            "Sepsis: yes",
                            "Total fluid balance before AKI")
  
  beta_df$hazard_ratio = exp(beta_df$beta_bar)
  
  order = c("Age at admission",
            "Gender: male",
            "Race: non-white",
            "BMI",
            "Apache 3",
            "Charlson",
            "Days of mechanical ventilation",
            "Sepsis: yes",
            "Apache 3 * Sepsis: yes",
            "Total fluid balance before AKI",
            "Persistent severe AKI: yes"
  )
  
  cox_res_df = beta_df %>%
    mutate(variable_true = factor(variable_true, levels = order)) %>%
    arrange(variable_true) %>%
    select(variable_true, beta_bar, hazard_ratio, se_pool, z, p_value)
  colnames(cox_res_df) = c("Variable", "Coefficient", "Hazard Ratio", "Std Error", "Z-Statistics", "P-Value")
  write.csv(cox_res_df, paste("G:/Persistent_AKI/Xinlei/Result/cox_res_df", analysis, ".csv", sep = ""))
  
  
  rm(list = ls())
  
  
}

