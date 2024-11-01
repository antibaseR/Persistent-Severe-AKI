library(tidyverse)


rawdata = read.csv("G:/Persistent_AKI/Persistent_AKI_Data/revised_pAKI_Final_Data_v2.csv")
rawdata_sub = rawdata %>%
  filter(KIDNEY_TRANSP==1)
length(unique(rawdata$patientid))

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
kdigoDat = read.csv("G:/Persistent_AKI/Persistent_AKI_Data/KDIGO_Data.csv")
#check3 = kdigoDat %>%
#  filter(patientvisitid == 34593816)
#kdigoDat2 = read.csv("G:/Persistent_AKI/Persistent_AKI_Data/AKI_Subgroups.csv")
check4 = kdigoDat2 %>%
  filter(patientvisitid == 34593816)
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
rawData = rawData %>%
  select(-c("stage_2_3_less_24", "day1_kidgo", "day2_kidgo", "day3_kidgo",
            "day4_kidgo", "day5_kidgo", "day6_kidgo", "day7_kidgo", "day8_kidgo"))
rawData = merge(rawData, kdigoDat, by=c("patientid", "patientvisitid"), all.x=T)


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




# within patients classified as stage2/3, how many has stage2/3 for less than one day (using the first two days of KDIGO scores)
s23_less_24_list = vector()
s23_all_list = vector()
ptAKI_list = vector()
s23_no_less_24_list = vector()
pAKI_s23_list = vector()
psAKI_list = vector()

for (i in 1:1){
  
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
  data = data %>%
    filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
  data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter id
  data = data[!duplicated(data$patientid),]
  
  s23_all_list[i] = data %>%
    filter(aki_2_3_icu != 0) %>%
    nrow()
  
  s23_less_24_list[i] = data %>%
    filter(aki_2_3_icu != 0) %>%
    filter(day2_kidgo %in% c(1,0)) %>%
    nrow()
  
  ptAKI_list[i] = data %>%
    filter(AKI_Category_subgroup %in% c("Persistent_Severe_AKI", "No_Persistent_Severe_AKI")) %>%
    nrow()
  
  # previous definition
  s23_no_less_24_list[i] = data %>%
    filter(aki_2_3_icu != 0) %>%
    filter(day2_kidgo %in% c(2,3)) %>%
    nrow()
  
  pAKI_s23_list[i] = data %>%
    filter(aki_2_3_icu != 0) %>%
    filter(day2_kidgo %in% c(2,3)) %>%
    filter(AKI_Category_subgroup == "Persistent_Severe_AKI") %>%
    nrow()

  psAKI_list[i] = data %>%
    filter(AKI_Category_subgroup == "Persistent_Severe_AKI") %>%
    nrow()
  
  
}


mean(s23_all_list)
mean(s23_less_24_list)
mean(ptAKI_list)
mean(s23_no_less_24_list)
mean(pAKI_s23_list)
mean(psAKI_list)

#s23_no_less_24_data = data %>%
#  filter(aki_2_3_icu != 0) %>%
#  filter(day2_kidgo %in% c(2,3))

#s23_no_less_24_data2 = data %>%
#  filter(aki_2_3_icu != 0) %>%
#  filter(day2_kidgo %in% c(2,3))

# where the discrepancy comes from
# patients that are not removed in the old kdigo dataset but removed in the new kdigo dataset
check = s23_no_less_24_data2 %>%
  filter(!(patientvisitid %in% unique(s23_no_less_24_data$patientvisitid))) %>%
  select("patientvisitid", "AKI_Category_subgroup", "aki_2_3_icu", "day1_kidgo", "day2_kidgo", "day3_kidgo", "day4_kidgo", "day5_kidgo", "day6_kidgo", "day7_kidgo", "day8_kidgo")

# among those patients, the KIDGO in the new kdigo dataset
check2 = data %>%
  filter(patientvisitid %in% check$patientvisitid) %>%
  select("patientvisitid", "AKI_Category_subgroup", "aki_2_3_icu", "day1_kidgo", "day2_kidgo", "day3_kidgo", "day4_kidgo", "day5_kidgo", "day6_kidgo", "day7_kidgo", "day8_kidgo")

write.csv(check, "G:/Persistent_AKI/Xinlei/discrepancy_patients_kdigo_score1.csv")
write.csv(check2, "G:/Persistent_AKI/Xinlei/discrepancy_patients_kdigo_score2.csv")
