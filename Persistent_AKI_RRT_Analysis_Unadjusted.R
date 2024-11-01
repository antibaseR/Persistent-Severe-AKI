library(tidyverse)
library(caret)
library(pROC)

# setting
analysis = "_RRT_with_death_unadjusted"


# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the results
train_sum_list = list()
test_sum_list = list()
var_imp_list = list()
train_dat_list = list()
test_dat_list = list()
cart_mod_list = list()
logi_mod_list = list()

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

risk_var = c("AKI_Category_final",
             "rrt_discharge")



# create the correct death within admission variable
## if the death date is earlier than the ICU admission date than replace it with NA
rawData$death_date_sh_true = ifelse(difftime(rawData$death_date_sh, rawData$dt_icu_start_sh, units="days")<0, NA, rawData$death_date_sh)
## if the survived subjects' death date is earlier than or the same as the hospital discharge date, then the subjects died in the hospital
rawData$discharge_to_death_days = difftime(rawData$death_date_sh_true, rawData$encounter_end_date_time_sh, unit="days")
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days<=0 & rawData$dead==0), 1, as.numeric(as.character(rawData$dead)))
## if the dead subjects' death date is later than the hospital discharge date, the subjects survived in hospital
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days>0 & rawData$dead==1), 0, as.numeric(as.character(rawData$dead)))



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
  if (analysis == "_RRT" | analysis == "_RRT_with_death_unadjusted"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
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
  
  # reformat the outcome
  data = data %>%
    filter(AKI_Category_true == "Persistent_AKI" | AKI_Category_true == "Transient_AKI") %>%
    mutate(AKI_Category_final = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "1", "0")))
  
  data = data %>%
    mutate(race = ifelse(race==1, 0, 1))
  
  if (analysis == "_RRT" | analysis == "_RRT_unadjusted"){
    data = data %>%
      mutate(rrt_discharge = factor(ifelse(rrt_discharge==1, "Yes", "No")))
  }
  
  if (analysis == "_RRT_with_death"| analysis == "_RRT_with_death_unadjusted"){
    data = data %>%
      mutate(rrt_discharge = factor(ifelse((rrt_discharge==1|dead_hosp==1), "Yes", "No")))
  }
  
  # subset the data
  
  ## if it's the sensitivity analysis, remove the subjects who died within 48 hours of ICU admission
  if (analysis=="_RRT_sensitivity"){
    subdata = data %>%
      filter(!(difftime(dt_icu_end_sh, dt_icu_start_sh, units="days")<=2 & dead==1))
  } else {subdata = data}
  
  ## subset the data to include only the risk factors we're interested in
  subdata = subdata[, names(data) %in% risk_var]
  
  # check the missingness of the subdata
  missing = data.frame(colSums(is.na(subdata)) / nrow(subdata))
  colnames(missing) = c("missingness")
  
  # remove variables with high missingness
  highmiss_var = rownames(missing)[missing$missingness>=0.8]
  subdata = subdata[, !(names(subdata) %in% highmiss_var)]
  
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
  trainData = as.data.frame(train[, -(which("rrt_discharge"==colnames(train)))])
  colnames(trainData) = colnames(train)[-(which("rrt_discharge"==colnames(train)))]
  testLabel = test$rrt_discharge
  testData = as.data.frame(test[, -(which("rrt_discharge"==colnames(test)))])
  colnames(testData) = colnames(test)[-(which("rrt_discharge"==colnames(test)))]
  
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

saveRDS(train_sum_list, paste("G:/Persistent_AKI/Xinlei/Data/train_sum_list", analysis, ".rds", sep=""))
saveRDS(test_sum_list, paste("G:/Persistent_AKI/Xinlei/Data/test_sum_list", analysis, ".rds", sep=""))
saveRDS(var_imp_list, paste("G:/Persistent_AKI/Xinlei/Data/var_imp_list", analysis, ".rds", sep=""))
saveRDS(cart_mod_list, paste("G:/Persistent_AKI/Xinlei/Data/cart_mod_list", analysis, ".rds", sep=""))
saveRDS(train_dat_list, paste("G:/Persistent_AKI/Xinlei/Data/train_dat_list", analysis, ".rds", sep=""))
saveRDS(test_dat_list, paste("G:/Persistent_AKI/Xinlei/Data/test_dat_list", analysis, ".rds", sep=""))


# read the result
library(readxl)
library(tidyverse)

train_sum_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/train_sum_list", analysis, ".rds", sep=""))
test_sum_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/test_sum_list", analysis, ".rds", sep=""))
var_imp_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/var_imp_list", analysis, ".rds", sep=""))
cart_mod_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/cart_mod_list", analysis, ".rds", sep=""))
train_dat_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/train_dat_list", analysis, ".rds", sep=""))
test_dat_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Data/test_dat_list", analysis, ".rds", sep=""))


# the testing AUC
AUC_df = data.frame(glm=as.numeric())
for (i in 1:length(test_sum_list)){
  AUC_df[i,] = round(test_sum_list[[i]]$auc, 3)
}

colnames(AUC_df) = c("Logistic Regression")
write.csv(AUC_df, paste("G:/Persistent_AKI/Xinlei/Result/AUC_df", analysis, ".csv", sep=""))


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
