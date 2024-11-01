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
#analysis = c("_ageq_sensitivity")
analysis="_ageq"

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the results
pAKI_risk_among_AKI_ageq = vector()
pAKI_risk_among_all_ageq = vector()
pAKI_risk_among_AKI = vector()
pAKI_risk_among_all = vector()
numerator_list = vector()
denominator_AKI_list = vector()
denominator_all_list = vector()

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
             "AKI_Category_subgroup"
)



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
  
  nrow(data[data$ESRD==1,])
  nrow(data[data$ecmo==1,])
  nrow(data[data$BASELINE_eGFR<15,])
  nrow(data[data$refCreat>=4,])
  
  data = data %>%
    filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
  nrow(imp)-nrow(data)
  nrow(data)
  data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
  #subset to rows that are not duplicates of a previous encounter for that patient
  data = data[!duplicated(data$patientid),]
  nrow(data)
  
  # normalize the continuous variables
  cont_var = names(select_if(data, is.numeric))
  cont_var = cont_var[!(cont_var %in% c("patientid", "patientvisitid"))]
  data[, names(data)%in%cont_var] = scale(data[, names(data)%in%cont_var])
  
  factor_names = names(data)[sapply(data, is.factor)]
  #factor_names = factor_names[!(factor_names %in% c("AKI_Category", "AKI_Category_true", "AKI_Category_final"))]
  factor_names = factor_names[!(factor_names %in% c("AKI_Category", "AKI_Category_true", "AKI_Category_subgroup"))]
  data[, factor_names] = data.frame(lapply(data[, factor_names], function(x) as.numeric(as.character(x))))
  
  # subset the data
  
  ## if it's the sensitivity analysis, remove the subjects who died within 48 hours of ICU admission
  nrow(data)
  if (analysis=="_ageq_sensitivity"){
    subdata = data %>%
      filter(!(difftime(dt_icu_end_sh, dt_icu_start_sh, units="days")<=2 & dead==1))
  } else {subdata = data}
  
  #nrow(data)-nrow(subdata)
  #nrow(subdata)
  #table(subdata$AKI_Category_final)
  
  ## subset the data to include only the risk factors we're interested in
  subdata = subdata[, names(data) %in% risk_var]
  
  # check the missingness of the subdata
  missing = data.frame(colSums(is.na(subdata)) / nrow(subdata))
  colnames(missing) = c("missingness")
  
  # remove variables with high missingness
  highmiss_var = rownames(missing)[missing$missingness>=0.8]
  subdata = subdata[, !(names(subdata) %in% highmiss_var)]
  
  # take the complete cases of the dataset
  #completeData = subdata[complete.cases(subdata), ]
  completeData = subdata
  
  # risk of pAKI
  denominator_all = completeData %>%
    nrow()
  denominator_AKI =  completeData %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
    nrow()
  
  numerator = completeData %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI")) %>%
    nrow()
  
  pAKI_risk_among_all_new = numerator/denominator_all
  pAKI_risk_among_AKI_new = numerator/denominator_AKI
  
  pAKI_risk_among_all = append(pAKI_risk_among_all, pAKI_risk_among_all_new)
  pAKI_risk_among_AKI = append(pAKI_risk_among_AKI, pAKI_risk_among_AKI_new)
  
  # find the age quartiles
  Q0 = summary(completeData$age_at_admission)[1]
  Q1 = summary(completeData$age_at_admission)[2]
  Q2 = summary(completeData$age_at_admission)[3]
  Q3 = summary(completeData$age_at_admission)[5]
  Q4 = summary(completeData$age_at_admission)[6]
  
  # save as a vector
  age_quartiles = c(Q0, Q1, Q2, Q3, Q4)
  
  for (ageq in 1:4){
    
    completeData_sub = completeData %>%
      filter(age_at_admission < age_quartiles[(ageq+1)]) %>%
      filter(age_at_admission >= age_quartiles[(ageq)])
    
    denominator_all_ageq = completeData_sub %>%
      nrow()
    denominator_AKI_ageq =  completeData_sub %>%
      filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
      nrow()
    
    numerator_ageq = completeData_sub %>%
      filter(AKI_Category_subgroup == "Persistent_Severe_AKI") %>%
      nrow()
    
    numerator_list = append(numerator_list, numerator_ageq)
    denominator_AKI_list = append(denominator_AKI_list, denominator_AKI_ageq)
    denominator_all_list = append(denominator_all_list, denominator_all_ageq)
    
    #pAKI_risk_among_AKI_ageq_new = numerator_ageq/denominator_AKI_ageq
    #pAKI_risk_among_all_ageq_new = numerator_ageq/denominator_all_ageq
    
    #pAKI_risk_among_AKI_ageq = append(pAKI_risk_among_AKI_ageq, pAKI_risk_among_AKI_ageq_new)
    #pAKI_risk_among_all_ageq = append(pAKI_risk_among_all_ageq, pAKI_risk_among_all_ageq_new)
    
  
  }
  
}

# reformat the result
#pAKI_risk_among_AKI_ageq_mat = matrix(pAKI_risk_among_AKI_ageq, nrow=4, ncol=5)
#pAKI_risk_among_all_ageq_mat = matrix(pAKI_risk_among_all_ageq, nrow=4, ncol=5)

# pool the result
#pAKI_risk_among_AKI_ageq_pooled = rowMeans(pAKI_risk_among_AKI_ageq_mat)
#pAKI_risk_among_all_ageq_pooled = rowMeans(pAKI_risk_among_all_ageq_mat)
#pAKI_risk_among_AKI_pooled = mean(pAKI_risk_among_AKI)
#pAKI_risk_among_all_pooled = mean(pAKI_risk_among_all)

numerator_mat = matrix(numerator_list, nrow=4, ncol=5)
numerator_mat_pooled = rowMeans(numerator_mat)

denominator_AKI_mat = matrix(denominator_AKI_list, nrow=4, ncol=5)
denominator_AKI_mat_pooled = rowMeans(denominator_AKI_mat)
rest_AKI = denominator_AKI_mat_pooled - numerator_mat_pooled
chisq_dat_AKI = t(matrix(c(numerator_mat_pooled, rest_AKI), ncol=2))
chisq.test(chisq_dat_AKI)

denominator_all_mat = matrix(denominator_all_list, nrow=4, ncol=5)
denominator_all_mat_pooled = rowMeans(denominator_all_mat)
rest_all = denominator_all_mat_pooled - numerator_mat_pooled
chisq_dat_all = t(matrix(c(numerator_mat_pooled, rest_all), ncol=2))
chisq.test(chisq_dat_all)



