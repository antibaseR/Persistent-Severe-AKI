library(tidyverse)
library(table1)
options(scipen = 999)

### summary statistics after imputation ###
# read the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# imputed variables
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

# get each imputed dataset
impDat1 = impDat[[1]] %>% mutate(imp = 1) %>% select(imp, everything()) %>% arrange(patientid, patientvisitid)
impDat2 = impDat[[2]] %>% mutate(imp = 2) %>% select(imp, everything()) %>% arrange(patientid, patientvisitid)
impDat3 = impDat[[3]] %>% mutate(imp = 3) %>% select(imp, everything()) %>% arrange(patientid, patientvisitid)
impDat4 = impDat[[4]] %>% mutate(imp = 4) %>% select(imp, everything()) %>% arrange(patientid, patientvisitid)
impDat5 = impDat[[5]] %>% mutate(imp = 5) %>% select(imp, everything()) %>% arrange(patientid, patientvisitid)

# replace the imp_var in the original dataset with the imputed ones
dat1 = rawData %>% arrange(patientid, patientvisitid)
dat1 = dat1 %>%
  select(-all_of(imp_var))
imp_sub = impDat1 %>%
  select(patientid, patientvisitid, all_of(imp_var))
dat1 = merge(dat1, imp_sub, by = c("patientid", "patientvisitid"))
dat1$imp = 1

dat2 = rawData %>% arrange(patientid, patientvisitid)
dat2 = dat2 %>%
  select(-all_of(imp_var))
imp_sub = impDat2 %>%
  select(patientid, patientvisitid, all_of(imp_var))
dat2 = merge(dat2, imp_sub, by = c("patientid", "patientvisitid"))
dat2$imp = 2

dat3 = rawData %>% arrange(patientid, patientvisitid)
dat3 = dat3 %>%
  select(-all_of(imp_var))
imp_sub = impDat3 %>%
  select(patientid, patientvisitid, all_of(imp_var))
dat3 = merge(dat3, imp_sub, by = c("patientid", "patientvisitid"))
dat3$imp = 3

dat4 = rawData %>% arrange(patientid, patientvisitid)
dat4 = dat4 %>%
  select(-all_of(imp_var))
imp_sub = impDat4 %>%
  select(patientid, patientvisitid, all_of(imp_var))
dat4 = merge(dat4, imp_sub, by = c("patientid", "patientvisitid"))
dat4$imp = 4

dat5 = rawData %>% arrange(patientid, patientvisitid)
dat5 = dat5 %>%
  select(-all_of(imp_var))
imp_sub = impDat5 %>%
  select(patientid, patientvisitid, all_of(imp_var))
dat5 = merge(dat5, imp_sub, by = c("patientid", "patientvisitid"))
dat5$imp = 5


dat1 = dat1 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0 & ICU_LOS>=3)
dat1 = dat1[order(dat1$encounter_start_date_time_sh),] #order dat1frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat1 = dat1[!duplicated(dat1$patientid),]
dat1 = dat1 %>%
  filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
  mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No"), 
                                        levels=c("No", "Yes"))) %>%
  select(imp, ICU_LOS, HOSP_LOS, AKI_Category_subgroup)

dat2 = dat2 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0 & ICU_LOS>=3)
dat2 = dat2[order(dat2$encounter_start_date_time_sh),] #order dat2frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat2 = dat2[!duplicated(dat2$patientid),]
dat2 = dat2 %>%
  filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
  mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No"), 
                                        levels=c("No", "Yes"))) %>%
  select(imp, ICU_LOS, HOSP_LOS, AKI_Category_subgroup)

dat3 = dat3 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0 & ICU_LOS>=3)
dat3 = dat3[order(dat3$encounter_start_date_time_sh),] #order dat3frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat3 = dat3[!duplicated(dat3$patientid),]
dat3 = dat3 %>%
  filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
  mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No"), 
                                        levels=c("No", "Yes"))) %>%
  select(imp, ICU_LOS, HOSP_LOS, AKI_Category_subgroup)

dat4 = dat4 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0 & ICU_LOS>=3)
dat4 = dat4[order(dat4$encounter_start_date_time_sh),] #order dat4frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat4 = dat4[!duplicated(dat4$patientid),]
dat4 = dat4 %>%
  filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
  mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No"), 
                                        levels=c("No", "Yes"))) %>%
  select(imp, ICU_LOS, HOSP_LOS, AKI_Category_subgroup)

dat5 = dat5 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0 & ICU_LOS>=3)
dat5 = dat5[order(dat5$encounter_start_date_time_sh),] #order dat5frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat5 = dat5[!duplicated(dat5$patientid),]
dat5 = dat5 %>%
  filter(!(AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI"))) %>%
  mutate(AKI_Category_subgroup = factor(ifelse(AKI_Category_subgroup =="Persistent_Severe_AKI", "Yes", "No"), 
                                        levels=c("No", "Yes"))) %>%
  select(imp, ICU_LOS, HOSP_LOS, AKI_Category_subgroup)

# combine the imputed dataset
df.list = list(dat1, dat2, dat3, dat4, dat5)
imputed_data = df.list %>% reduce(rbind)

# get the names of numeric variables
num_var = names(select_if(imputed_data, is.numeric))
num_var = num_var[ !num_var == 'imp']
num_data = cbind(imputed_data$imp, imputed_data[, names(imputed_data) %in% num_var])
colnames(num_data)[1] = 'imp'

# get AKI category index
index1 = imputed_data$AKI_Category_subgroup == 'Yes'
index0 = imputed_data$AKI_Category_subgroup == 'No'


by = 'all'
# summarize data by imp group
for (by in c('Persistent_Severe_AKI', 'No_Persistent_Severe_AKI', 'all')){
  
  if (by == 'No_Persistent_Severe_AKI') {
    num_data_sub = num_data[index0,]
  } else if (by == 'Persistent_Severe_AKI') {
    num_data_sub = num_data[index1,]
  } else {
    num_data_sub = num_data
  }
  
  # for numeric
  num_sum_by_imp = num_data_sub %>%
    group_by(imp) %>%
    summarise_all(.funs = list(min = ~ min(., na.rm=T),
                               Q1 = ~ quantile(x = ., probs = 0.25, na.rm=T),
                               median = ~ median(., na.rm=T),
                               Q3 = ~ quantile(x = ., probs = 0.75, na.rm=T),
                               max = ~ max(., na.rm=T),
                               mean = ~ mean(., na.rm=T),
                               sd = ~ sd(., na.rm=T)
    ))
  
  num_sum_by_imp = num_sum_by_imp %>% select(-imp)
  # get the mean of the statistics
  num_sum = tibble(variable = rep(num_var, times = 7),
                   statistics = rep(c('min', 'Q1', 'median', 'Q3', 'max', 'mean', 'sd'), each = length(num_var)),
                   value = round(colMeans(num_sum_by_imp), 3))
  num_sum$variable = factor(num_sum$variable, levels = num_var)
  num_sum = num_sum %>% arrange(variable)
  num_sum = spread(num_sum, key = statistics, value = value)
  num_sum = num_sum %>% select(variable, min, Q1, median, Q3, max, mean, sd)
  
  
  # save the summary tables
  write.csv(num_sum, paste('G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/Table1_numeric_', by, '.csv', sep=''))
}




num_data_sub_0 = num_data[index0,]
num_data_sub_1 = num_data[index1,]

wilcox.test(num_data_sub_0$ICU_LOS, num_data_sub_1$ICU_LOS, alternative = "two.sided")
wilcox.test(num_data_sub_0$HOSP_LOS, num_data_sub_1$HOSP_LOS, alternative = "two.sided")


