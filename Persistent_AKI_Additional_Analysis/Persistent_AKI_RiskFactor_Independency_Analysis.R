library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
library(caret)
library(pROC)
library(corrplot)



# setting
analysis="_riskfactor_independency"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the model fit
train_sum_list = list()
test_sum_list = list()
train_sum_var_list = list()
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


## variables we need to adjust for the adjusted logistic regression model
adj_var = list()
adj_var_true = list()
adj_var_order = list()

#-	Cumulative balance before AKI
#-	Fluid input
#-	Minimum platelets
#-	Min arterial pH

cov = c("apache3",
        "MAX_LACTATE",
        "MV_DAYS",
        "sofa_24"
        )

cov_true = c("Apache3 * Max SOFA",
             "Apache 3",
             "Max lactate",
             "Days of mechanical ventilation",
             "Max SOFA")

adj_var[[1]] = cov
adj_var[[2]] = cov
adj_var[[3]] = cov
adj_var[[4]] = cov

cum_balance_true = c("Intercept",
                     "Apache3 * Max SOFA",
                     "Apache 3",
                     "Cumulative balance before AKI",
                     "Max lactate",
                     "Days of mechanical ventilation",
                     "Max SOFA")

fluid_input_true = c("Intercept",
                     "Apache3 * Max SOFA",
                     "Apache 3",
                     "Max lactate",
                     "Days of mechanical ventilation",
                     "Max SOFA",
                     "Total fluids input")

plat_true = c("Intercept",
              "Apache3 * Max SOFA",
              "Apache 3",
              "Max lactate",
              "Min platelets",
              "Days of mechanical ventilation",
              "Max SOFA")

ph_true = c("Intercept",
            "Apache3 * Max SOFA",
            "Apache 3",
            "Max lactate",
            "Min arterial pH",
            "Days of mechanical ventilation",
            "Max SOFA")

adj_var_true[[1]] = cum_balance_true
adj_var_true[[2]] = fluid_input_true
adj_var_true[[3]] = plat_true
adj_var_true[[4]] = ph_true


cum_balance_order = c("Intercept",
                      "Apache 3",
                     "Max lactate",
                     "Days of mechanical ventilation",
                     "Max SOFA",
                     "Apache3 * Max SOFA",
                     "Cumulative balance before AKI")

fluid_input_order = c("Intercept",
                      "Apache 3",
                     "Max lactate",
                     "Days of mechanical ventilation",
                     "Max SOFA",
                     "Apache3 * Max SOFA",
                     "Total fluids input")

plat_order = c("Intercept",
               "Apache 3",
              "Max lactate",
              "Days of mechanical ventilation",
              "Max SOFA",
              "Apache3 * Max SOFA",
              "Min platelets")

ph_order = c("Intercept",
             "Apache 3",
            "Max lactate",
            "Days of mechanical ventilation",
            "Max SOFA",
            "Apache3 * Max SOFA",
            "Min arterial pH")

adj_var_order[[1]] = cum_balance_order
adj_var_order[[2]] = fluid_input_order
adj_var_order[[3]] = plat_order
adj_var_order[[4]] = ph_order

var = c("CUMULATIVE_BALANCE_BEFORE_AKI",
        "TOTAL_INPUTS_BEFORE_AKI",
        "MIN_PLATELETS",
        "MIN_PH")

var_label = c("Cumulative balance before AKI",
              "Total fluids input",
              "Min platelets",
              "Min arterial pH")


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
    filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
    mutate(AKI_Category_final = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No"), 
                                       levels = c("No", "Yes")))
  
  #corrdata = data %>%
  #  select(all_of(var), all_of(cov), AKI_Category_final) %>%
  #  mutate(AKI_Category_final = as.numeric(AKI_Category_final))
  #corr_plot = corrplot(cor(corrdata), method = "number")
  # No collinearity
  
  for (i_var in 1:length(var)){
    
    # subset the data to the variables we need
    subdata = data %>%
      select(all_of(adj_var[[i_var]]), all_of(var[i_var]), AKI_Category_final)
    
    # take the complete cases of the dataset
    completeData = subdata[complete.cases(subdata), ]
    
    # split the data into the training and testing set
    set.seed(960725)
    trainIndex <- createDataPartition(completeData$AKI_Category_final, 
                                      p = .8, 
                                      list = FALSE, 
                                      times = 1)
    #train = completeData[trainIndex,]
    #test = completeData[-trainIndex,]
    
    #trainLabel = train$AKI_Category_final
    #trainData = train[, -(which("AKI_Category_final"==colnames(train)))]
    #testLabel = test$AKI_Category_final
    #testData = test[, -(which("AKI_Category_final"==colnames(test)))]
    
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
    
    # include interaction into the model
    completeData$apache_sofa = completeData$apache3 * completeData$sofa_24
    
    train = completeData[trainIndex,]
    test = completeData[-trainIndex,]
    
    trainLabel = train$AKI_Category_final
    trainData = train[, -(which("AKI_Category_final"==colnames(train)))]
    testLabel = test$AKI_Category_final
    testData = test[, -(which("AKI_Category_final"==colnames(test)))]
    
    # run a logistic regression
    logi.model = train(x = trainData,
                       y = trainLabel,
                       method = "glm", 
                       family = "binomial",
                       trControl = control,
                       metric = "ROC")
    train_sum_var_list[[i_var]] = logi.model
  }
  
  # append the fitted object to the list
  train_sum_list[[i]] = train_sum_var_list
  
}


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
  beta_df$Vt = beta_df$Vb + beta_df$Vw + (beta_df$Vb/length(unique(glm_fit_df$imp)))
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
  write.csv(glm_res_df, paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/glm_res_df", analysis, "_", var_label[j], ".csv", sep = ""))
  
}
