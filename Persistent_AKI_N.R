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
Q0_list = vector()
Q1_list = vector()
Q2_list = vector()
Q3_list = vector()
Q4_list = vector()
n1_list = vector()
n2_list = vector()
n3_list = vector()
n4_list = vector()
Q0_dead_list = vector()
Q1_dead_list = vector()
Q2_dead_list = vector()
Q3_dead_list = vector()
Q4_dead_list = vector()
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

## if the death date is earlier than the ICU admission date than replace it with NA
rawData$death_date_sh_true = ifelse(difftime(rawData$death_date_sh, rawData$dt_icu_start_sh, units="days")<0, NA, rawData$death_date_sh)
## if the survived subjects' death date is earlier than or the same as the hospital discharge date, then the subjects died in the hospital
rawData$discharge_to_death_days = difftime(rawData$death_date_sh_true, rawData$encounter_end_date_time_sh, unit="days")
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days<=0 & rawData$dead==0), 1, as.numeric(as.character(rawData$dead)))
## if the dead subjects' death date is later than the hospital discharge date, the subjects survived in hospital
rawData$dead_hosp = ifelse((!is.na(rawData$discharge_to_death_days) & rawData$discharge_to_death_days>0 & rawData$dead==1), 0, as.numeric(as.character(rawData$dead)))

## number of days from ICU admission to death
rawData$icu_to_death_days = difftime(rawData$death_date_sh_true, rawData$dt_icu_start_sh, units="days")
## number of days from ICU admission to hospital discharge
rawData$icu_to_discharge_days = difftime(rawData$encounter_end_date_time_sh, rawData$dt_icu_start_sh, units="days")
## replace missing date with discharge date
rawData$survived_days = ifelse((is.na(rawData$death_date_sh_true)), rawData$icu_to_discharge_days, rawData$icu_to_death_days)

## 90-day mortality: 
### - if the death date is not NA and the icu to death days is >= 90 -- survived
### - if the death date is not NA and the icu to death days is < 90 -- dead
### - if the death date is NA that means we didn't observe the event time but we know the result, the data is censored, so the survival = dead_hosp and time = icu_to_discharge_days
rawData$dead_90 = ifelse(!(is.na(rawData$death_date_sh_true)) & rawData$survived_days>=90, 0, 
                         ifelse(!(is.na(rawData$death_date_sh_true)) & rawData$survived_days<90, 1, rawData$dead_hosp))
rawData$survival_time = ifelse(rawData$survived_days>90, 90, rawData$survived_days)
rawData$status = ifelse((is.na(rawData$death_date_sh_true) & rawData$dead_90==1), 3, rawData$dead_90)
rawData$time = ifelse(rawData$status==3, 0, rawData$survival_time)
rawData$time2 = rawData$survival_time


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
  nrow(data[data$CKD_stage_4==1,])
  nrow(data[data$CKD_stage_5==1,])
  
  data = data %>%
    filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
  nrow(imp)-nrow(data)
  nrow(data)
  
  # number of pAKI and tAKI
  pAKI_new = data %>%
    filter(AKI_Category=="Persistent_AKI") %>%
    nrow()
  pAKI_list[i] = pAKI_new
  
  tAKI_new = data %>%
    filter(AKI_Category=="Transient_AKI") %>%
    nrow()
  tAKI_list[i] = tAKI_new
  
  data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter id
  #subset to rows that are not duplicates of a previous encounter for that patient
  data = data[!duplicated(data$patientid),]
  nrow(data)
  
  
  # normalize the continuous variables
  cont_var = names(select_if(data, is.numeric))
  cont_var = cont_var[!(cont_var %in% c("patientid", "patientvisitid"))]

  # reformat the outcome
  data$AKI_Category = as.character(data$AKI_Category)
  data$AKI_Category_true = ifelse((data$aki_1_icu==1 & data$AKI_Category=="no_AKI"), "stage_1_AKI", data$AKI_Category)
  data$AKI_Category = factor(data$AKI_Category, 
                             levels = c("no_AKI", "Transient_AKI", "Persistent_AKI"))
  data$AKI_Category_true = factor(data$AKI_Category_true, 
                                  levels = c("no_AKI", "stage_1_AKI", "Transient_AKI", "Persistent_AKI"))
  table(data$AKI_Category_true)
  
  data = data %>%
    filter(AKI_Category_true == "Persistent_AKI" | AKI_Category_true == "Transient_AKI") %>%
    mutate(AKI_Category_final = factor(ifelse(AKI_Category_true =="Persistent_AKI", "Yes", "No")))
  
  #factor_names = names(data)[sapply(data, is.factor)]
  #factor_names = factor_names[!(factor_names %in% c("AKI_Category", "AKI_Category_true", "AKI_Category_final"))]
  #data[, factor_names] = data.frame(lapply(data[, factor_names], function(x) as.numeric(as.character(x))))
  
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
  #completeData = subdata
  
  # find the age quartiles
  Q0 = summary(completeData$age_at_admission)[1]
  Q1 = summary(completeData$age_at_admission)[2]
  Q2 = summary(completeData$age_at_admission)[3]
  Q3 = summary(completeData$age_at_admission)[5]
  Q4 = summary(completeData$age_at_admission)[6]
  
  # save as a vector
  age_quartiles = c(Q0, Q1, Q2, Q3, Q4)
  n1 = completeData %>%
    filter(age_at_admission>=Q0 & age_at_admission<Q1) %>%
    nrow()
  
  n2 = completeData %>%
    filter(age_at_admission>=Q1 & age_at_admission<Q2) %>%
    nrow()
  
  n3 = completeData %>%
    filter(age_at_admission>=Q2 & age_at_admission<Q3) %>%
    nrow()
  
  n4 = completeData %>%
    filter(age_at_admission>=Q3 & age_at_admission<Q4) %>%
    nrow()
  
  # 90-day mortality by age quartiles
  Q1_dead = completeData %>%
    filter(age_at_admission>=Q0 & age_at_admission<Q1) %>%
    filter(dead_90==1) %>%
    nrow()
  
  Q2_dead = completeData %>%
    filter(age_at_admission>=Q1 & age_at_admission<Q2) %>%
    filter(dead_90==1) %>%
    nrow()
  
  Q3_dead = completeData %>%
    filter(age_at_admission>=Q2 & age_at_admission<Q3) %>%
    filter(dead_90==1) %>%
    nrow()
  
  Q4_dead = completeData %>%
    filter(age_at_admission>=Q3 & age_at_admission<Q4) %>%
    filter(dead_90==1) %>%
    nrow()

  Q0_list = append(Q0_list, Q0)
  Q1_list = append(Q1_list, Q1)
  Q2_list = append(Q2_list, Q2)
  Q3_list = append(Q3_list, Q3)
  Q4_list = append(Q4_list, Q4)
  
  n1_list = append(n1_list, n1)
  n2_list = append(n2_list, n2)
  n3_list = append(n3_list, n3)
  n4_list = append(n4_list, n4)
  
  Q1_dead_list = append(Q1_dead_list, Q1_dead)
  Q2_dead_list = append(Q2_dead_list, Q2_dead)
  Q3_dead_list = append(Q3_dead_list, Q3_dead)
  Q4_dead_list = append(Q4_dead_list, Q4_dead)
  
}


mean(Q0_list)
mean(Q1_list)
mean(Q2_list)
mean(Q3_list)
mean(Q4_list)

mean(pAKI_list)

mean(Q1_dead_list)
mean(Q2_dead_list)
mean(Q3_dead_list)
mean(Q4_dead_list)

mean(n1_list)
mean(n2_list)
mean(n3_list)
mean(n4_list)

mean(Q1_dead_list)/mean(n1_list)
mean(Q2_dead_list)/mean(n2_list)
mean(Q3_dead_list)/mean(n3_list)
mean(Q4_dead_list)/mean(n4_list)

mean(pAKI_list)
mean(tAKI_list)