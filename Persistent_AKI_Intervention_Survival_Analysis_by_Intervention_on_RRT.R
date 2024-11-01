library(tidyverse)
library(caret)
library(pROC)



# setting
analysis_list= c("_intervention_by_each_rrt_only", "_intervention_by_each_rrt_death")

for (analysis in analysis_list){
  
  # load the imputed data
  impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
  rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
  
  # create the correct death within admission variable
  ## if the death date is earlier than the ICU admission date than replace it with NA
  rawData$death_date_sh_true = ifelse(difftime(rawData$death_date_sh, rawData$dt_icu_start_sh, units="days")<0, NA, rawData$death_date_sh)
  ## if the survived subjects' death date is earlier than or the same as the hospital discharge date, then the subjects died in the hospital
  rawData$discharge_to_death_days = difftime(rawData$death_date_sh_true, rawData$encounter_end_date_time_sh, unit="days")
  rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days<=0 & rawData$dead==0), 1, as.numeric(as.character(rawData$dead)))
  ## if the dead subjects' death date is later than the hospital discharge date, the subjects survived in hospital
  rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days>0 & rawData$dead==1), 0, as.numeric(as.character(rawData$dead)))
  
  # create list to save the results
  train_sum_list = list()
  test_sum_list = list()
  var_imp_list = list()
  train_dat_list = list()
  test_dat_list = list()
  cart_mod_list = list()
  logi_mod_list = list()
  
  train_sum_int_list = list()
  test_sum_int_list = list()
  var_imp_int_list = list()
  train_dat_int_list = list()
  test_dat_int_list = list()
  cart_mod_int_list = list()
  logi_mod_int_list = list()
  
  
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
  risk_var = c("patientid",
               "patientvisitid",
               "encounter_start_date_time_sh",
               "encounter_end_date_time_sh", 
               "age_at_admission",
               "gender",
               "charlson",
               "apache3",
               "MV_DAYS",
               "refCreat",
               "rrt_discharge",
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
               "rrt_discharge")
  
  vaso_cov_true = c("Intercept",
                    "Age at admission",
                    "Apache3",
                    "Gender: male",
                    "Minimum MAP before AKI",
                    "Septiic shock: yes",
                    "Total fluid input before AKI",
                    "Early vasopressor/inotrope"
  )
  
  vaso_cov_order = c("Intercept",
                     "Age at admission",
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
                 "rrt_discharge")
  
  saline_cov_true = c("Intercept",
                      "Age at admission",
                      "Apache 3",
                      "Gender: male",
                      "Minimum MAP before AKI",
                      "Saline avoidance before AKI",
                      "Total fluid input before AKI")
  
  saline_cov_order = c("Intercept",
                       "Age at admission",
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
              "rrt_discharge")
  
  diu_cov_true = c("Intercept",
                   "Age at admission",
                   "Apache 3",
                   "Congestive heart failure: yes",
                   "Diuretic avoidance before AKI",
                   "Gender: male",
                   "Mild liver disease: yes",
                   "Moderate or severe liver disease: yes"
  )
  
  diu_cov_order = c("Intercept",
                    "Age at admission",
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
                 "rrt_discharge")
  
  nephro_cov_true = c("Intercept",
                      "Age at admission",
                      "Apache 3",
                      "Gender: male",
                      "Machanical ventilation: yes",
                      "Nephrotoxin avoidance before AKI")
  
  nephro_cov_order = c("Intercept",
                       "Age at admission",
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
    if(analysis %in% c("_intervention_by_each_rrt_only", "_intervention_by_each_rrt_death")){
      data = data %>%
        filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
      data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
      #subset to rows that are not duplicates of a previous encounter for that patient
      data = data[!duplicated(data$patientid),]
    }
    
    data = data %>%
      mutate(mv = ifelse(MV_DAYS>0, 1, 0))
    
    # transform categorical variables into factors
    unique_value = apply(data, 2, function(x) length(unique(x)))
    cat_cols = unique_value<=4
    cat_var = data[, cat_cols] %>% colnames()
    cat_var = cat_var[cat_var!="MIN_SCVO2"]
    data[, cat_var] = lapply(data[, cat_var], factor)
    
    # normalize the continuous variables
    cont_var = names(select_if(data, is.numeric))
    cont_var = cont_var[!(cont_var %in% c("patientid", "patientvisitid"))]
    data[, names(data)%in%cont_var] = scale(data[, names(data)%in%cont_var])
    
    if (analysis == "_intervention_by_each_rrt_only"){
      data = data %>%
        mutate(rrt_discharge = factor(ifelse(rrt_discharge==1, "Yes", "No")))
    }
    
    if (analysis == "_intervention_by_each_rrt_death"){
      data = data %>%
        mutate(rrt_discharge = factor(ifelse((rrt_discharge==1|dead_hosp==1), "Yes", "No")))
    }
    
    # Cox proportional hazard model
    for (int in 1:length(int_var)){
      
      subdata = data %>%
        filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
        mutate(race = ifelse(race==1, 0, 1)) %>%
        select(all_of(adj_var[[int]]), all_of(int_var[int]))
      
      # take the complete cases of the dataset
      completeData = subdata[complete.cases(subdata), ]
      
      # split the data into the training and testing set
      set.seed(960725)
      trainIndex <- createDataPartition(completeData$rrt_discharge, 
                                        p = .8, 
                                        list = FALSE, 
                                        times = 1)
      train = completeData[trainIndex,]
      test = completeData[-trainIndex,]
      
      trainLabel = train$rrt_discharge
      trainData = train[, -(which("rrt_discharge"==colnames(train)))]
      testLabel = test$rrt_discharge
      testData = test[, -(which("rrt_discharge"==colnames(test)))]
      
      # model to predict
      control <- trainControl(method = "repeatedcv", number = 5, repeats = 1, 
                              search = "grid", savePredictions = "final", 
                              summaryFunction = twoClassSummary, classProbs = TRUE, 
                              verboseIter = TRUE)
      
      # run a single rpart tree to find interactions
      #cart.model = train(x = trainData,
      #                   y = trainLabel,
      #                   method = "rpart", 
      #                   trControl = control,
      #                   metric = "ROC")
      #plot(cart.model$finalModel, uniform=TRUE,
      #     main="Classification Tree")
      #text(cart.model$finalModel, use.n.=TRUE, all=TRUE, cex=.8)
      #cart_mod_list[[i]] = cart.model
      
      # run a logistic regression
      logi.model = train(x = trainData,
                         y = trainLabel,
                         method = "glm", 
                         family = "binomial",
                         trControl = control,
                         metric = "ROC")
      train_sum_int_list[[int]] = logi.model
      var_imp_int_list[[int]] = varImp(logi.model)
      
      # evaluating on the testing set
      pred_prob = predict(logi.model, newdata = testData, type = "prob")
      test_sum_int_list[[int]] = roc(testLabel, pred_prob$Yes) # Draw ROC curve.
      
      train_dat_int_list[[int]] = train
      test_dat_int_list[[int]] = test
      
    }
    # append the fitted object to the list
    train_sum_list[[i]] = train_sum_int_list
    test_sum_list[[i]] = test_sum_int_list
  }
  
  # combine the AUC from the model
  for (j in 1:4){
    AUC_df = data.frame(glm=as.numeric())
    for (i in 1:length(test_sum_list)){
      AUC_df[i,] = round(test_sum_list[[i]][[j]]$auc, 3)
    }
    write.csv(AUC_df, paste("G:/Persistent_AKI/Xinlei/Result/AUC_df", analysis, "_", int_var[j], ".csv", sep = ""))
  }
  
  # combine the results from the glm model
  for (j in 1:4){
    glm_fit_df = data.frame()
    for (i in 1:length(train_sum_list)){
      glm_fit = summary(train_sum_list[[i]][[j]])
      glm_fit_df_new = data.frame(variable = rownames(glm_fit$coefficients),
                                  imp = i,
                                  beta = glm_fit$coefficients[,c("Estimate")],
                                  se = glm_fit$coefficients[,c("Std. Error")])
      glm_fit_df = rbind(glm_fit_df, glm_fit_df_new)
    }
    rownames(glm_fit_df) = seq(1, nrow(glm_fit_df), 1)
    
    beta_df = glm_fit_df %>%
      group_by(variable) %>%
      summarise(beta_bar = mean(beta))
    
    glm_fit_df = merge(glm_fit_df, beta_df, by="variable")
    
    Vw_df = glm_fit_df %>% 
      group_by(variable) %>%
      summarise(Vw = mean(se^2))
    
    Vb_df = glm_fit_df %>%
      group_by(variable) %>%
      summarise(Vb = sum((beta-beta_bar)^2)/(length(unique(glm_fit_df$imp))-1))
    
    beta_df = merge(merge(beta_df, Vw_df, by="variable"), Vb_df, by="variable")
    beta_df$Vt = beta_df$Vb + beta_df$Vw
    beta_df$se_pool = sqrt(beta_df$Vt)
    beta_df$z = beta_df$beta_bar/beta_df$se_pool
    beta_df$p_value = 2*pnorm(-abs(beta_df$z))
    
    beta_df$variable_true = adj_var_true[[j]]
    
    beta_df$odds_ratio = exp(beta_df$beta_bar)
    
    order = adj_var_order[[j]]
    
    glm_res_df = beta_df %>%
      mutate(variable_true = factor(variable_true, levels = order)) %>%
      arrange(variable_true) %>%
      select(variable_true, beta_bar, odds_ratio, se_pool, z, p_value)
    colnames(glm_res_df) = c("Variable", "Coefficient", "Odds Ratio", "Std Error", "Z-Statistics", "P-Value")
    write.csv(glm_res_df, paste("G:/Persistent_AKI/Xinlei/Result/glm_res_df", analysis, "_", int_var[j], ".csv", sep = ""))
  }
  rm(list = ls())
}