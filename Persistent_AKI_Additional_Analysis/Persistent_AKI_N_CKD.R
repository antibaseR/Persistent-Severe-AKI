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
#analysis = c("_sensitivity")
analysis=""

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the results
CKD_list = vector()
CKD_wo_refCreat_list = vector()
refCreat_list = vector()


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
            #'refCreat',
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
             "AKI_Category_final",
             "dead_90")


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
  
  data = data %>%
    filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & CKD_stage_4==0 & CKD_stage_5==0)
  data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter id
  #subset to rows that are not duplicates of a previous encounter for that patient
  data = data[!duplicated(data$patientid),]
  
  CKD_new = data %>%
    filter(CKD == 1) %>%
    nrow()
  
  CKD_wo_refCreat_new = data %>%
    filter(CKD == 1 & is.na(refCreat)) %>%
    nrow()
  
  refCreat_new = data %>%
    filter(!is.na(refCreat)) %>%
    nrow()
 
  CKD_list = append(CKD_list, CKD_new)
  CKD_wo_refCreat_list = append(CKD_wo_refCreat_list, CKD_wo_refCreat_new)
  refCreat_list = append(refCreat_list, refCreat_new)
}

mean(CKD_list)
mean(CKD_wo_refCreat_list)
mean(refCreat_list)

