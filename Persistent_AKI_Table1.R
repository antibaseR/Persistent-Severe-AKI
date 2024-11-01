library(tidyverse)
library(table1)
options(scipen = 999)

### summary statistics after imputation ###
# read the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')
landmarkData = read.csv('G:/Persistent_AKI/Persistent_AKI_Data/Landmark_Data.csv')

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

# missingness of each variable
missing = data.frame(colSums(is.na(rawData)) / nrow(rawData))
colnames(missing) = c("missingness")

# variables that needs to be removed
## variables that have too much missingness are excluded from the variables list
highmiss_remove = missing %>%
  filter(missingness>=0.8) %>%
  rownames()

# create MV indicator
rawData = rawData %>%
  mutate(mv = factor(ifelse(MV_DAYS>0, 1, 0)))

# reformat the date variables
rawData = rawData %>%
  mutate(aki_2_3_icu_dt_sh = as.POSIXct(aki_2_3_icu_dt_sh, format = "%Y-%m-%dT%H:%M:%OSZ"),
         encounter_start_date_time_sh = as.POSIXct(encounter_start_date_time_sh, format = "%Y-%m-%dT%H:%M:%OSZ"),
         encounter_end_date_time_sh = as.POSIXct(encounter_end_date_time_sh, format = "%Y-%m-%dT%H:%M:%OSZ"),
         dt_icu_start_sh = as.POSIXct(dt_icu_start_sh, format = "%Y-%m-%dT%H:%M:%OSZ"),
         dt_icu_end_sh = as.POSIXct(dt_icu_end_sh, format = "%Y-%m-%dT%H:%M:%OSZ"),
         FIRST_AKI_DATETIME = as.POSIXct(FIRST_AKI_DATETIME, format = "%Y-%m-%dT%H:%M:%OSZ"))

# create time from hospital admission to stage2/3 AKI
rawData = rawData %>%
  mutate(encounter_to_stage2_3_time = as.numeric(difftime(aki_2_3_icu_dt_sh, 
                                                          encounter_start_date_time_sh, 
                                                          units="days")))
# create time from ICU admission to stage2/3 AKI
rawData = rawData %>%
  mutate(icu_to_stage2_3_time = as.numeric(difftime(aki_2_3_icu_dt_sh,
                                                    dt_icu_start_sh, 
                                                    units="days")))
# create time from ICU admission to ICU admission
rawData = rawData %>%
  mutate(encounter_to_icu_time = as.numeric(difftime(dt_icu_start_sh, 
                                                     encounter_start_date_time_sh,
                                                     units="days")))
# create time from hospital admission to any AKI
rawData = rawData %>%
  mutate(encounter_to_any_AKI_time = as.numeric(difftime(FIRST_AKI_DATETIME,
                                                         encounter_start_date_time_sh, 
                                                         units="days")))

# create time from ICU admission to any AKI
rawData = rawData %>%
  mutate(icu_to_any_AKI_time = as.numeric(difftime(FIRST_AKI_DATETIME,
                                                         dt_icu_start_sh, 
                                                         units="days")))

# check mld and mlsd, and create CLD
check = rawData %>%
  select(patientid, patientvisitid, mld, msld) %>%
  filter(msld==1 & mld==0)
rawData = rawData %>%
  mutate(CLD = factor(ifelse((msld==1|mld==1), 1, 0)))

# create 90-day mortality
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
rawData$dead_90 = factor(rawData$dead_90)

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
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
dat1 = dat1[order(dat1$encounter_start_date_time_sh),] #order dat1frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat1 = dat1[!duplicated(dat1$patientid),]

dat2 = dat2 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
dat2 = dat2[order(dat2$encounter_start_date_time_sh),] #order dat2frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat2 = dat2[!duplicated(dat2$patientid),]

dat3 = dat3 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
dat3 = dat3[order(dat3$encounter_start_date_time_sh),] #order dat3frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat3 = dat3[!duplicated(dat3$patientid),]

dat4 = dat4 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
dat4 = dat4[order(dat4$encounter_start_date_time_sh),] #order dat4frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat4 = dat4[!duplicated(dat4$patientid),]

dat5 = dat5 %>%
  filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
dat5 = dat5[order(dat5$encounter_start_date_time_sh),] #order dat5frame by encounter time
#subset to rows that are not duplicates of a previous encounter for that patient
dat5 = dat5[!duplicated(dat5$patientid),]

# check the weird time
check = dat5 %>% 
  filter(encounter_to_icu_time<0) %>%
  select(patientid, patientvisitid, encounter_start_date_time_sh, dt_icu_start_sh, encounter_to_icu_time)


# combine the imputed dataset
df.list = list(dat1, dat2, dat3, dat4, dat5)
imputed_data = df.list %>% reduce(rbind)

# get the names of numeric and categorical variables
num_var = names(select_if(imputed_data, is.numeric))
num_var = num_var[ !num_var == 'imp']
cat_var = names(select_if(imputed_data, is.factor))


# seperate the dataset into numeric and categorical only
num_data = cbind(imputed_data$imp, imputed_data[, names(imputed_data) %in% num_var])
colnames(num_data)[1] = 'imp'
cat_data = cbind(imputed_data$imp, imputed_data[, names(imputed_data) %in% cat_var])
colnames(cat_data)[1] = 'imp'

# get AKI category index
index_p = imputed_data$AKI_Category_subgroup == 'Persistent_Severe_AKI'
index_t = imputed_data$AKI_Category_subgroup == 'No_Persistent_Severe_AKI'
index_s1 = imputed_data$AKI_Category_subgroup == 'stage_1_AKI'
index_n = imputed_data$AKI_Category_subgroup == 'no_AKI'

by = 'all'
# summarize data by imp group
for (by in c('Persistent_Severe_AKI', 'No_Persistent_Severe_AKI', 'stage_1_AKI', 'no_AKI', 'all')){
  
  if (by == 'No_Persistent_Severe_AKI') {
    num_data_sub = num_data[index_t,]
    cat_data_sub = cat_data[index_t,]
  } else if (by == 'Persistent_Severe_AKI') {
    num_data_sub = num_data[index_p,]
    cat_data_sub = cat_data[index_p,]
  } else if (by == "stage_1_AKI"){
    num_data_sub = num_data[index_s1,]
    cat_data_sub = cat_data[index_s1,]
  } else if (by == 'no_AKI') {
    num_data_sub = num_data[index_n,]
    cat_data_sub = cat_data[index_n,]
  } else {
    num_data_sub = num_data
    cat_data_sub = cat_data
  }
  
  # for numeric
  num_sum_by_imp = num_data_sub %>%
    group_by(imp) %>%
    summarise_all(.funs = list(min = ~ min(., na.rm=T),
                               Q1 = ~ quantile(x = ., probs = 0.25, na.rm=T),
                               median = ~ median(., na.rm=T),
                               Q3 = ~ quantile(x = ., probs = 0.75, na.rm=T),
                               P90 = ~ quantile(x = ., probs = 0.90, na.rm=T),
                               max = ~ max(., na.rm=T),
                               mean = ~ mean(., na.rm=T),
                               sd = ~ sd(., na.rm=T)
    ))
  
  num_sum_by_imp = num_sum_by_imp %>% select(-imp)
  # get the mean of the statistics
  num_sum = tibble(variable = rep(num_var, times = 8),
                   statistics = rep(c('min', 'Q1', 'median', 'Q3', 'P90', 'max', 'mean', 'sd'), each = length(num_var)),
                   value = round(colMeans(num_sum_by_imp), 3))
  num_sum$variable = factor(num_sum$variable, levels = num_var)
  num_sum = num_sum %>% arrange(variable)
  num_sum = spread(num_sum, key = statistics, value = value)
  num_sum = num_sum %>% select(variable, min, Q1, median, Q3, P90, max, mean, sd)
  
  # for categorical
  summary_df <- data.frame(
    Variable = character(),
    Level = character(),
    Frequency = integer(),
    Percentage = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (col in cat_var) {
    # create a frequency table and calculate percentages
    freq_table <- table(cat_data_sub[[col]])
    percent_table <- prop.table(freq_table) * 100
    # combine the frequency and percentage tables into a data frame
    summary_table <- data.frame(
      Variable = rep(col, length(freq_table)),
      Level = names(freq_table),
      Frequency = as.vector(freq_table),
      Percentage = as.vector(percent_table),
      stringsAsFactors = FALSE
    )
    # add the summary table to the summary data frame
    summary_df <- rbind(summary_df, summary_table)
  }
  
  cat_sum = summary_df %>% mutate(Percentage = round(Percentage, 2),
                                  Frequency = Frequency/length(unique(cat_data$imp)))
  
  # save the summary tables
  write.csv(num_sum, paste('G:/Persistent_AKI/Xinlei/Summary Statistics/Table1_numeric_', by, '.csv', sep=''))
  write.csv(cat_sum, paste('G:/Persistent_AKI/Xinlei/Summary Statistics/Table1_categorical_', by, '.csv', sep=''))
}


