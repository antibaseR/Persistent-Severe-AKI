library(tidyverse)
library(caret)
library(xgboost)
library(caretEnsemble)
library(pROC)
library(caTools)
library(randomForest)
library(nnet)
library(DALEX)
library(iml)

# setting
for (analysis in c("", "_sensitivity")){
  
  # load the imputed data
  impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
  rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
  
  # create list to save the results
  train_sum_list = list()
  test_sum_list = list()
  var_imp_list = list()
  shap_val_list = list()
  caret_mod_list = list()
  ensemble_mod_list = list()
  train_dat_list = list()
  test_dat_list = list()
  
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
               "AKI_Category_final")
  
  
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
    
    # normalize the continuous variables
    cont_var = names(select_if(data, is.numeric))
    cont_var = cont_var[!(cont_var %in% c("patientid", "patientvisitid"))]
    data[, names(data)%in%cont_var] = scale(data[, names(data)%in%cont_var])
    
    data = data %>%
      filter(AKI_Category_true == "Persistent_AKI" | AKI_Category_true == "Transient_AKI") %>%
      mutate(AKI_Category_final = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No")))
    
    factor_names = names(data)[sapply(data, is.factor)]
    factor_names = factor_names[!(factor_names %in% c("AKI_Category", "AKI_Category_true", "AKI_Category_final"))]
    data[, factor_names] = data.frame(lapply(data[, factor_names], function(x) as.numeric(as.character(x))))
    
    # subset the data
    
    ## if it's the sensitivity analysis, remove the subjects who died within 48 hours of ICU admission
    nrow(data)
    if (analysis=="_sensitivity"){
      subdata = data %>%
        filter(!(difftime(dt_icu_end_sh, dt_icu_start_sh, units="days")<=2 & dead==1))
    } else {subdata = data}
    
    nrow(data)-nrow(subdata)
    nrow(subdata)
    table(subdata$AKI_Category_final)
    
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
    trainIndex = createDataPartition(completeData$AKI_Category_final, 
                                     p = .8, 
                                     list = FALSE, 
                                     times = 1)
    
    train = completeData[trainIndex,]
    test = completeData[-trainIndex,]
    
    trainLabel = train$AKI_Category_final
    trainData = train[, -(which("AKI_Category_final"==colnames(train)))]
    testLabel = test$AKI_Category_final
    testData = test[, -(which("AKI_Category_final"==colnames(test)))]
    
    # model to predict
    control <- trainControl(method = "repeatedcv", number = 5, repeats = 1, 
                            search = "grid", savePredictions = "final", 
                            summaryFunction = twoClassSummary, classProbs = TRUE, 
                            verboseIter = TRUE)
    # list of algorithms to use in ensemble
    glmnet_grid = expand.grid(alpha=1,
                              #lambda=2^runif(5, min = -10, 3))
                              lambda = 0.00377)
    
    xgb_grid = expand.grid(nrounds=100, 
                           max_depth=c(4), 
                           eta=c(0.2), 
                           gamma=100, 
                           colsample_bytree=c(0.4), 
                           min_child_weight=100, 
                           subsample=c(0.75))
    
    model_list = caretList(x = trainData, y = trainLabel, trControl = control, 
                           metric = "ROC", 
                           tuneList=list(
                             glmnet=caretModelSpec(method="glmnet", tuneGrid=data.frame(glmnet_grid)),
                             xgb=caretModelSpec(method="xgbTree", tuneGrid=data.frame(xgb_grid))
                           ))
    caret_mod_list[[i]] = model_list
    
    # results
    #xyplot(resamples(model_list))
    #res = resamples(model_list)
    #summary(res)
    
    # stack 
    ensembleControl = trainControl(method = "repeatedcv", number = 5, repeats = 1, 
                                   savePredictions = TRUE, summaryFunction=twoClassSummary, 
                                   classProbs = TRUE, verboseIter = TRUE)
    ensemble = caretEnsemble(model_list, metric = "ROC", trControl = ensembleControl)
    ensemble_mod_list[[i]] = ensemble
    train_sum_list[[i]] = summary(ensemble)
    var_imp_list[[i]] = varImp(ensemble)
    
    # evaluating on the testing set
    model_preds = lapply(model_list, predict, newdata=test, type="prob")
    model_preds = lapply(model_preds, function(x) x[, "Yes"])
    model_preds = data.frame(model_preds)
    ens_preds = predict(ensemble, newdata=test, type="prob")
    model_preds$ensemble = ens_preds
    test_sum_list[[i]] = caTools::colAUC(model_preds, test$AKI_Category_final)
    
    train_dat_list[[i]] = train
    test_dat_list[[i]] = test
    
    # calculate the SHAPLEY value 
    ## this is done in the cluster with R.script in my laptop
    #SHAP_df = data.frame()
    #predictor = Predictor$new(ensemble, data = testData, y = testLabel)
    #testData_sample = testData
    #for (j in 1:nrow(testData_sample)){
    #  shapley = Shapley$new(predictor, x.interest = testData_sample[j,])
    #  shapley_values = shapley$results
    #  SHAP_df_new = shapley_values %>% mutate(subject=j)
    #  SHAP_df = rbind(SHAP_df, SHAP_df_new)
    #  print(j)
    #}
    #shap_val_list[[i]] = SHAP_df
    
  }
  
  saveRDS(train_sum_list, paste("G:/Persistent_AKI/Xinlei/Data/train_sum_list", analysis, ".rds", sep=""))
  saveRDS(test_sum_list, paste("G:/Persistent_AKI/Xinlei/Data/test_sum_list", analysis, ".rds", sep=""))
  saveRDS(var_imp_list, paste("G:/Persistent_AKI/Xinlei/Data/var_imp_list", analysis, ".rds", sep=""))
  saveRDS(caret_mod_list, paste("G:/Persistent_AKI/Xinlei/Data/caret_mod_list", analysis, ".rds", sep=""))
  saveRDS(ensemble_mod_list, paste("G:/Persistent_AKI/Xinlei/Data/ensemble_mod_list", analysis, ".rds", sep=""))
  saveRDS(train_dat_list, paste("G:/Persistent_AKI/Xinlei/Data/train_dat_list", analysis, ".rds", sep=""))
  saveRDS(test_dat_list, paste("G:/Persistent_AKI/Xinlei/Data/test_dat_list", analysis, ".rds", sep=""))
  
  rm(list = ls())
}