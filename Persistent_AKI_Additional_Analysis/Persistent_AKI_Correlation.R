library(corrplot)
library(tidyverse)

impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

renal_recovery_data = list()
for (i in 1:5){
  renal_recovery_data[[i]] = read_csv(paste("G:/Persistent_AKI/Xinlei/renal_recovery_data_imputed_", i, ".csv", sep=""))
}

# create list to save the model fit
km_fit_list = list()
logrank_fit_list = list()
cox_fit_list = list()

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

# variables we need in the analysis
all_var = c("age_at_admission",
            "gender",
            "race",
            "charlson",
            "apache3",
            "BMI",
            "TOTAL_INPUTS_BEFORE_AKI",
            "MV_DAYS",
            "SEPSIS",
            "ACTUTE_CORONARY_SYNDROME",
            "LEFT_RIGHT_BIVENTRICULAR_FAILURE",
            "CARDIOGENIC_SHOCK",
            "DISSEMINATED_INTRAVASCULAR_COAGULATION",
            "CARDIAC_ARREST",
            "ACUTE_RESPIRATORY_FAILURE",
            "CEREBRAL_VASCULAR_EVENT",
            "HEMORRHAGE",
            "PULMONARY_EMBOLISM",
            "DVT",
            "refCreat",
            "sofa_24",
            "MAX_LACTATE",
            "TOTAL_OUTPUTS_BEFORE_AKI",
            "AKI_Category_subgroup")


i=1

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
data = data %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
data = data[!duplicated(data$patientid),]
data = data[!is.na(data$renal_recov), ]

subdata = data[, c(names(data) %in% all_var)] %>%
  mutate(TOTAL_BALANCE_BEFORE_AKI = TOTAL_INPUTS_BEFORE_AKI - TOTAL_OUTPUTS_BEFORE_AKI) %>%
  select(-c(TOTAL_INPUTS_BEFORE_AKI,TOTAL_OUTPUTS_BEFORE_AKI))

subdata = subdata %>%
  mutate(AKI_Category_subgroup = ifelse(AKI_Category_subgroup == "Persistent_Severe_AKI", 1, 0)) %>%
  mutate(AKI_Category_subgroup = as.numeric(AKI_Category_subgroup))

subdata = data.frame(lapply(subdata, function(x) as.numeric(as.character(x))))

names(subdata) = c("Days of mechanical ventilation",
                   "Acute coronary syndrome",
                   "Left right biventricular failure",
                   "Cardiogenic shock",
                   "DIC",
                   "Cardiac arrest",
                   "Acute respiratory failure",
                   "CVA",
                   "Hemorrhage",
                   "Pulmonary embolism",
                   "DVT",
                   "Sepsis",
                   "Persistent severe AKI",
                   "Age",
                   "Gender",
                   "Race",
                   "Charlson",
                   "Apache 3",
                   "Max lactate",
                   "Reference creatinine",
                   "BMI",
                   "Max SOFA",
                   "Total fluid balance before AKI")

M = cor(subdata[complete.cases(subdata), ])
corrplot(M, method="color", tl.col="black") # colorful number








