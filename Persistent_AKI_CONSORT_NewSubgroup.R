library(tidyverse)

# setting
analysis="_CONSORT"

rawdata = read.csv("G:/Persistent_AKI/revised_pAKI_Final_Data_v2.csv")
rawdata_sub = rawdata %>%
  filter(KIDNEY_TRANSP==1)
length(unique(rawdata$patientid))

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
            'BASELINE_eGFR',
            'MIN_PULSE_PRESS_BEFORE_AKI',
            'MAX_PULSE_PRESS_BEFORE_AKI',
            'MIN_MAP_BEFORE_AKI')




ESRD_list = vector()
ECMO_list = vector()
eGFR_list = vector()
refCreat_list = vector()
CKD_4_list = vector()
CKD_5_list = vector()
qualified_list = vector()
first_encounter_only_list = vector()
pAKI_list = vector()
p3AKI_list = vector()
p2AKI_list = vector()
p1AKI_list = vector()
tAKI_list = vector()
stage1_AKI_list = vector()
nAKI_list = vector()


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
  
  
  #data = data[order(c(data$patientid, data$encounter_start_date_time_sh)),]
  #data[, names(imp)[names(imp)%in%imp_var]] = imp[, names(imp)[names(imp)%in%imp_var]]
  
  ## exclusion
  #The exclusions we need to have for the persistent AKI primary analysis are:
  #1. Kidney transplant
  #2. End stage Renal disease (by ICD-9/10 codes or ICD9/10 codes equivalent to CKD stage 4 or 5)
  #3. eGFR <15 at admission
  #4. Creatinine ≥ 4 mg/dL at admission
  #5. ECMO
  #6. Only use the first encounter
  ESRD = data %>%
    filter(ESRD==1) %>%
    nrow()
  
  ECMO = data %>%
    filter(ecmo==1) %>%
    nrow()
  
  eGFR = data %>%
    filter(BASELINE_eGFR<15) %>%
    nrow()
  
  refCreat = data %>%
    filter(refCreat>=4) %>%
    nrow()
  
  CKD_4 = data %>%
    filter(CKD_stage_4==1) %>%
    nrow()
  
  CKD_5 = data %>%
    filter(CKD_stage_5==1) %>%
    nrow()
  
  data = data %>%
    filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
  
  qualified = data %>%
    nrow()
  
  data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter id
  data = data[!duplicated(data$patientid),]
  
  #check = data %>%
  #  select(patientid, patientvisitid, KDIGO_pattern, AKI_Category_subgroup_new, dead)
  #write.csv(check, "G:/Persistent_AKI/Xinlei/sample.csv")
  
  first_encounter_only = data %>%
    nrow()
  
  # reformat the outcome
  pAKI = data %>%
    filter(AKI_Category_subgroup_new == "Persistent_Severe_AKI") %>%
    nrow()
  
  p3AKI = data %>%
    filter(AKI_Category_subgroup_new == "Persistent AKI stage 3") %>%
    nrow()
  
  p2AKI = data %>%
    filter(AKI_Category_subgroup_new == "Persistent AKI stage 2") %>%
    nrow()
  
  p1AKI = data %>%
    filter(AKI_Category_subgroup_new == "Persistent AKI stage 1") %>%
    nrow()
  
  tAKI = data %>%
    filter(AKI_Category_subgroup_new == "Transient_AKI") %>%
    nrow()
  
  nAKI = data %>%
    filter(AKI_Category_subgroup_new == "no_AKI") %>%
    nrow()
  
  
  ESRD_list[i] = ESRD
  ECMO_list[i] = ECMO
  eGFR_list[i] = eGFR
  refCreat_list[i] = refCreat
  CKD_4_list[i] = CKD_4
  CKD_5_list[i] = CKD_5
  qualified_list[i] = qualified
  285130-qualified_list
  first_encounter_only_list[i] = first_encounter_only
  pAKI_list[i] = pAKI
  p3AKI_list[i] = p3AKI
  p2AKI_list[i] = p2AKI
  p1AKI_list[i] = p1AKI
  tAKI_list[i] = tAKI
  nAKI_list[i] = nAKI
  
  
}


ESRD_excluded = mean(ESRD_list)
ECMO_excluded = mean(ECMO_list)
eGFR_excluded = mean(eGFR_list)
refCreat_excluded = mean(refCreat_list)
CKD_4_excluded = mean(CKD_4_list)
CKD_5_excluded = mean(CKD_5_list)
qualified_excluded = mean(qualified_list)
first_encounter_only_excluded = mean(first_encounter_only_list)
pAKI_list_excluded = mean(pAKI_list)
p3AKI_list_excluded = mean(p3AKI_list)
p2AKI_list_excluded = mean(p2AKI_list)
p1AKI_list_excluded = mean(p1AKI_list)
tAKI_list_excluded = mean(tAKI_list)
nAKI_excluded = mean(nAKI_list)


check = data %>%
  filter(AKI_Category_subgroup_new == "Persistent AKI stage 1") %>%
  select(KDIGO_pattern, AKI_Category_subgroup_new, AKI_Category_subgroup)






