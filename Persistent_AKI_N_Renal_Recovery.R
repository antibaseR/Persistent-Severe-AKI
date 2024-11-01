library(tidyverse)
library(caret)
library(pROC)

# setting
analysis="_renal_recovery"
#analysis="_renal_recovery_sensitivity"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

renal_recovery_data = list()
for (i in 1:5){
  renal_recovery_data[[i]] = read_csv(paste("G:/Persistent_AKI/Xinlei/renal_recovery_data_imputed_", i, ".csv", sep=""))
}

# create list to save the results
renal_pAKI_list = vector()
renal_tAKI_list = vector()
pAKI_list = vector()
tAKI_list = vector()

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
             "refCreat",
             "AKI_Category_final",
             "renal_recov")


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
  
  # merge the renal recovery data
  data = merge(data, renal_recovery_data[[i]], by=c("patientid", "patientvisitid"), all.x=T)
  
  ## exclusion
  #The exclusions we need to have for the persistent AKI primary analysis are:
  #1. Kidney transplant
  #2. End stage Renal disease (by ICD-9/10 codes or ICD9/10 codes equivalent to CKD stage 4 or 5)
  #3. eGFR <15 at admission
  #4. Creatinine ≥ 4 mg/dL at admission
  #5. ECMO
  #6. Only use the first encounter
  if (analysis == "_renal_recovery"){
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
  
  # reformat the outcome
  data = data %>%
    filter(AKI_Category_true == "Persistent_AKI" | AKI_Category_true == "Transient_AKI") %>%
    mutate(AKI_Category_final = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", 1, 0)))
  
  data = data %>%
    mutate(renal_recov = factor(ifelse(renal_recov == 1, "Yes", "No")))
  
  ## if it's the sensitivity analysis, remove the subjects who died within 48 hours of ICU admission
  if (analysis=="_renal_recovery_sensitivity"){
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
  
  # the number of encounters 
  ## with persistent AKI
  pAKI = completeData %>%
    filter(AKI_Category_final == 1) %>%
    nrow()
  pAKI_list[i] = pAKI
  ## with transient AKI
  tAKI = completeData %>%
    filter(AKI_Category_final == 0) %>%
    nrow()
  tAKI_list[i] = tAKI
  ## with both the renal recovery and the persistent AKI
  renal_pAKI = completeData %>%
    filter(renal_recov=="Yes" & AKI_Category_final == 1) %>%
    nrow()
  renal_pAKI_list[i] = renal_pAKI
  ## with both the renal recovery and the transient AKI
  renal_tAKI = completeData %>%
    filter(renal_recov=="Yes" & AKI_Category_final == 0) %>%
    nrow()
  renal_tAKI_list[i] = renal_tAKI
  
}

# % of readmission among patients with pAKI
mean(renal_pAKI_list)/mean(pAKI_list)
# % of readmission among patients with tAKI
mean(renal_tAKI_list)/mean(tAKI_list)


