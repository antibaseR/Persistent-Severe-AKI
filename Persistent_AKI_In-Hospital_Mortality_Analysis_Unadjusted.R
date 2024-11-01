library(tidyverse)
library(caret)
library(pROC)

# setting
analysis = c("_hosp_mortality_unadjusted")

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# check if the death happened after AKI
#time_to_death = difftime(rawData$death_date_sh, rawData$encounter_start_date_time_sh, unit="secs")
#check = tibble(time_to_death = time_to_death, time_to_AKI = rawData$TIME_TO_FIRST_AKI)
#check$death_after_pAKI = ifelse(check$time_to_death>check$time_to_AKI, 1, 0)
#mean(check$death_after_pAKI==1, na.rm=T) # 99% death happened after the development of AKI

# create list to save the results
train_sum_list = list()
test_sum_list = list()
var_imp_list = list()
train_dat_list = list()
test_dat_list = list()
cart_mod_list = list()
logi_mod_list = list()

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
             "AKI_Category_final",
             "CKD",
             "gender_bmi")


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
  if (analysis == "_hosp_mortality_unadjusted"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  # reformat the outcome
  data = data %>%
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    mutate(AKI_Category_final = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "1", "0")))
  
  # subset the data
  if (analysis=="_hosp_mortality_unadjusted_sensitivity"){
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
  
  # Cox proportional hazard model
  ## variables we need to adjust for the adjusted Cox proportional hazard model
  adj_var = c("AKI_Category_final",
              "dead_hosp")
  
  subdata = subdata[, c(names(subdata) %in% adj_var)]
  
  
  # take the complete cases of the dataset
  completeData = subdata[complete.cases(subdata), ]
  completeData = completeData %>%
    mutate(dead_hosp = ifelse(dead_hosp==1, "Yes", "No"))
  
  # split the data into the training and testing set
  set.seed(960725)
  trainIndex <- createDataPartition(completeData$dead_hosp, 
                                    p = .8, 
                                    list = FALSE, 
                                    times = 1)
  train = completeData[trainIndex,]
  test = completeData[-trainIndex,]
  
  trainLabel = train$dead_hosp
  trainData = as.data.frame(train[, -(which("dead_hosp"==colnames(train)))])
  colnames(trainData) = colnames(train)[-(which("dead_hosp"==colnames(train)))]
  testLabel = test$dead_hosp
  testData = as.data.frame(test[, -(which("dead_hosp"==colnames(test)))])
  colnames(testData) = colnames(test)[-(which("dead_hosp"==colnames(test)))]
  
  # model to predict
  control <- trainControl(method = "repeatedcv", number = 5, repeats = 1, 
                          search = "grid", savePredictions = "final", 
                          summaryFunction = twoClassSummary, classProbs = TRUE, 
                          verboseIter = TRUE)
  
  # run a logistic regression
  logi.model = train(x = trainData,
                     y = trainLabel,
                     method = "glm", 
                     family = "binomial",
                     trControl = control,
                     metric = "ROC")
  train_sum_list[[i]] = logi.model
  var_imp_list[[i]] = varImp(logi.model)
  
  # evaluating on the testing set
  pred_prob = predict(logi.model, newdata = testData, type = "prob")
  test_sum_list[[i]] = roc(testLabel, pred_prob$Yes) # Draw ROC curve.
  
  train_dat_list[[i]] = train
  test_dat_list[[i]] = test
  
}



# combine the results from the glm model
glm_fit_df = data.frame()
for (i in 1:length(train_sum_list)){
  glm_fit = summary(train_sum_list[[i]])
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

beta_df$variable_true = c("Intercept",
                          "Persistent severe AKI: yes")

beta_df$odds_ratio = exp(beta_df$beta_bar)

order = c("Intercept",
          "Persistent severe AKI: yes")

glm_res_df = beta_df %>%
  mutate(variable_true = factor(variable_true, levels = order)) %>%
  arrange(variable_true) %>%
  select(variable_true, beta_bar, odds_ratio, se_pool, z, p_value)
colnames(glm_res_df) = c("Variable", "Coefficient", "Odds Ratio", "Std Error", "Z-Statistics", "P-Value")
write.csv(glm_res_df, paste("G:/Persistent_AKI/Xinlei/Result/glm_res_df", analysis, ".csv", sep = ""))




