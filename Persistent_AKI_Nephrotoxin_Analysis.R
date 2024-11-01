library(tidyverse)
library(caret)
library(pROC)

analysis = "_nephro"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the results
table_list = list()
chisq_test_list = list()

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
             "charlson",
             "apache3",
             "MV_DAYS",
             "refCreat",
             "AKI_Category_final",
             "NEPHRO_BEFORE_AKI",
             "NEPHROTOXIN_INTERVENTION")


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
  if (analysis == "_nephro"){
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
  data$AKI_Category = as.character(data$AKI_Category)
  data$AKI_Category_true = ifelse((!is.na(data$aki_1_icu) & data$AKI_Category=="no_AKI"), "stage_1_AKI", data$AKI_Category)
  data$AKI_Category = factor(data$AKI_Category, 
                             levels = c("no_AKI", "Transient_AKI", "Persistent_AKI"))
  data$AKI_Category_true = factor(data$AKI_Category_true, 
                                  levels = c("no_AKI", "stage_1_AKI", "Transient_AKI", "Persistent_AKI"))
  data = data %>%
    filter(AKI_Category_true == "Persistent_AKI" | AKI_Category_true == "Transient_AKI") %>%
    mutate(AKI_Category_final = factor(ifelse(AKI_Category_true =="Persistent_AKI", "Yes", "No")))
  
  # subset the data
  
  ## if it's the sensitivity analysis, remove the subjects who died within 48 hours of ICU admission
  if (analysis=="_xxx_sensitivity"){
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
  
  sum(completeData$AKI_Category_final=="Yes")
  sum(completeData$AKI_Category_final=="No")
  
  library(data.table)
  nephro_pAKI_table = table(completeData$NEPHRO_BEFORE_AKI, completeData$AKI_Category_final)
  chisq_test = chisq.test(nephro_pAKI_table)
  
  table_list[[i]] = nephro_pAKI_table
  chisq_test_list[[i]] = chisq_test
  
}
