library(tidyverse)

# create a list to save all the 5 dataset
Final_Data_List = list()

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

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
            'MIN_CARDIAC_INDEX',
            'MIN_CARDIAC_OUTPUT',
            'MAX_CVP',
            'MIN_CVP',
            'MEAN_DIASTOLIC_PRESS',
            'MAX_DIASTOLIC_PRESS',
            'MEAN_SYSTOLIC_PRESS',
            'MAX_SYSTOLIC_PRESS',
            'MIN_SCVO2',
            'MIN_SVO2',
            'BASELINE_eGFR',
            'MIN_PULSE_PRESS',
            'MAX_PULSE_PRESS',
            'MIN_MAP')


for (i in 1:length(impDat)){
  # get each imputed dataset
  imp = impDat[[i]]
  
  # replace the imputed variables
  data = rawData
  data[, names(imp)[names(imp)%in%imp_var]] = imp[, names(imp)[names(imp)%in%imp_var]]
  
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
  data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter id
  data = data[!duplicated(data$patientid),]

  # reformat the outcome
  data$AKI_Category = as.character(data$AKI_Category)
  data$AKI_Category_true = ifelse((data$aki_1_icu==1 & data$AKI_Category=="no_AKI"), "stage_1_AKI", data$AKI_Category)
  data$AKI_Category = factor(data$AKI_Category, 
                             levels = c("no_AKI", "Transient_AKI", "Persistent_AKI"))
  data$AKI_Category_true = factor(data$AKI_Category_true, 
                                  levels = c("no_AKI", "stage_1_AKI", "Transient_AKI", "Persistent_AKI")) 
  Final_Data_List[[i]] = data
}

saveRDS(Final_Data_List, "G:/Persistent_AKI/Final_Imputed_Data.rds")












