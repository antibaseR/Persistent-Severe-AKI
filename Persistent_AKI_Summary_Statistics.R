library(tidyverse)
library(table1)
options(scipen = 999)

### check the number of subjects in each AKI category ###
setwd("G:/Persistent_AKI")
rawdata = read.csv("revised_pAKI_Final_Data_v2.csv")
data = rawdata %>%
  filter(KIDNEY_TRANSP==0)

# the number of patientvisit in each AKI category (the true no AKI, AKI stage 1, 2, and 3. The current no AKI in the dataset has subjects with stage 1 AKI)
subdata = data %>%
  select(patientid, patientvisitid, AKI_Category, aki_1_icu, aki_2_3_icu)

## # of visits that truly has no AKI (no_AKI and stage1==NA): 128348
subdata %>%
  filter(AKI_Category=='no_AKI' & is.na(aki_1_icu)) %>%
  summarise(num_subject = n_distinct(patientvisitid))


## # of visits with stage 1 AKI: 146065
subdata %>%
  filter(aki_1_icu==1) %>%
  summarise(num_subject = n_distinct(patientvisitid))

## # of visits with stage 2 AKI: 88253
subdata %>%
  filter(aki_2_3_icu==2) %>%
  summarise(num_subject = n_distinct(patientvisitid))

## # of visits with stage 3 AKI: 21701
subdata %>%
  filter(aki_2_3_icu==3) %>%
  summarise(num_subject = n_distinct(patientvisitid))

## # of visits with stage 1 and 2: 86015
subdata %>%
  filter(aki_1_icu==1 & aki_2_3_icu==2) %>%
  summarise(num_subject = n_distinct(patientvisitid))
## # of visits with stage 1 and 3: 13904
subdata %>%
  filter(aki_1_icu==1 & aki_2_3_icu==3) %>%
  summarise(num_subject = n_distinct(patientvisitid))

## # of visits with any AKI: 156100
subdata %>%
  filter(!is.na(aki_1_icu) | !is.na(aki_2_3_icu)) %>%
  summarise(num_subject = n_distinct(patientvisitid))

## Proportion of persistent AKI
### denominator as any AKI: 8.44%
100 * sum(subdata$AKI_Category=='Persistent_AKI')/ 156100
### denominator as stage 2 and 3 AKI: 11.98%
subdata %>%
  filter(!is.na(aki_2_3_icu)) %>%
  summarise(num_subject = n_distinct(patientvisitid))
100 * sum(subdata$AKI_Category=='Persistent_AKI')/109954






















### summary statistics before imputation ##
# change the working directory to my folder
setwd("G:/Persistent_AKI/Xinlei")

# read the data
data = readRDS('./Data/Data_Revised.rds')

# get the names of numeric and categorical variables
num_var = names(select_if(data, is.numeric))
cat_var = names(select_if(data, is.factor))

# seperate the dataset into numeric and categorical only
num_data = data[, names(data) %in% num_var]
cat_data =data[, names(data) %in% cat_var]

# get AKI category index
index_p = data$AKI_Category == 'Persistent_AKI'
index_t = data$AKI_Category == 'Transient_AKI'
index_n = data$AKI_Category == 'no_AKI'

by = 'all'
# summarize data by imp group
for (by in c('Persistent_AKI', 'Transient_AKI', 'no_AKI', 'all')){
  
  if (by == 'Transient_AKI') {
    num_data_sub = num_data[index_t,]
    cat_data_sub = cat_data[index_t,]
  } else if (by == 'Persistent_AKI') {
    num_data_sub = num_data[index_p,]
    cat_data_sub = cat_data[index_p,]
  } else if (by == 'no_AKI') {
    num_data_sub = num_data[index_n,]
    cat_data_sub = cat_data[index_n,]
  } else {
    num_data_sub = num_data
    cat_data_sub = cat_data
  }
  
  # for numeric
  num_sum = num_data_sub %>%
    summarise_all(.funs = list(min = ~ min(., na.rm=T),
                               Q1 = ~ quantile(x = ., probs = 0.25, na.rm=T),
                               median = ~ median(., na.rm=T),
                               Q3 = ~ quantile(x = ., probs = 0.75, na.rm=T),
                               max = ~ max(., na.rm=T),
                               mean = ~ mean(., na.rm=T),
                               sd = ~ sd(., na.rm=T))
    )
  
  # get the mean of the statistics
  num_sum = tibble(Variable = rep(num_var, times = 7),
                   Statistics = rep(c('min', 'Q1', 'median', 'Q3', 'max', 'mean', 'sd'), each = length(num_var)),
                   Value = round(colMeans(num_sum), 3))
  num_sum$variable = factor(num_sum$Variable, levels = num_var)
  num_sum = num_sum %>% arrange(Variable)
  num_sum = spread(num_sum, key = Statistics, value = Value)
  num_sum = num_sum %>% select(Variable, min, Q1, median, Q3, max, mean, sd)
  
  num_missing = data.frame(colSums(is.na(num_data_sub)) / nrow(num_data_sub))*100
  num_missing$Missing = colSums(is.na(num_data_sub))
  num_missing$Variable = rownames(num_missing)
  colnames(num_missing) = c("Missing Pct.", "Missing", "Variable")
  num_missing = num_missing %>% select(Variable, Missing, everything())
  rownames(num_missing) = NULL
  
  num_sum = merge(num_sum, num_missing, 'Variable')
  num_sum$Variable = factor(num_sum$Variable, levels = num_missing$Variable)
  num_sum = num_sum[order(num_sum$Variable),]
  
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
  
  cat_missing = data.frame(colSums(is.na(cat_data_sub)) / nrow(cat_data_sub))*100
  cat_missing$Missing = colSums(is.na(cat_data_sub))
  cat_missing$Variable = rownames(cat_missing)
  colnames(cat_missing) = c("Missing Pct.", "Missing", "Variable")
  cat_missing = cat_missing %>% select(Variable, Missing, everything())
  rownames(cat_missing) = NULL
  
  cat_sum = merge(summary_df, cat_missing, 'Variable') %>% mutate(Percentage = round(Percentage, 2))
  cat_sum$Variable = factor(cat_sum$Variable, levels = cat_missing$Variable)
  cat_sum = cat_sum[order(cat_sum$Variable),]
  
  # save the summary tables
  write.csv(num_sum, paste('./Summary Statistics/summary_numeric_raw_', by, '.csv', sep=''))
  write.csv(cat_sum, paste('./Summary Statistics/summary_categorical_raw_', by, '.csv', sep=''))
}






















### summary statistics after imputation ###
# read the imputed data
data = readRDS('./Data/Data_Revised.rds')
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")

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

# get each imputed dataset
impDat1 = impDat[[1]] %>% mutate(imp = 1) %>% select(imp, everything())
impDat2 = impDat[[2]] %>% mutate(imp = 2) %>% select(imp, everything())
impDat3 = impDat[[3]] %>% mutate(imp = 3) %>% select(imp, everything())
impDat4 = impDat[[4]] %>% mutate(imp = 4) %>% select(imp, everything())
impDat5 = impDat[[5]] %>% mutate(imp = 5) %>% select(imp, everything())

# missingness of each variable
missing = data.frame(colSums(is.na(data)) / nrow(data))
colnames(missing) = c("missingness")

# variables that needs to be removed
## variables that have too much missingness are excluded from the variables list
highmiss_remove = missing %>%
  filter(missingness>=0.8) %>%
  rownames()

# replace the imp_var in the original dataset with the imputed ones
#dat1 = data
#dat1[, names(impDat1)%in%imp_var] = impDat1[, names(impDat1)%in%imp_var]
#dat1$imp = 1
#dat2 = data
#dat2[, names(impDat2)%in%imp_var] = impDat2[, names(impDat2)%in%imp_var]
#dat2$imp = 2
#dat3 = data
#dat3[, names(impDat3)%in%imp_var] = impDat3[, names(impDat3)%in%imp_var]
#dat3$imp = 3
#dat4 = data
#dat4[, names(impDat4)%in%imp_var] = impDat4[, names(impDat4)%in%imp_var]
#dat4$imp = 4
#dat5 = data
#dat5[, names(impDat5)%in%imp_var] = impDat5[, names(impDat5)%in%imp_var]
#dat5$imp = 5

dat1 = impDat1[, names(impDat1)%in%imp_var]
dat1$imp = 1
dat2 = impDat2[, names(impDat2)%in%imp_var]
dat2$imp = 2
dat3 = impDat3[, names(impDat3)%in%imp_var]
dat3$imp = 3
dat4 = impDat4[, names(impDat4)%in%imp_var]
dat4$imp = 4
dat5 = impDat5[, names(impDat5)%in%imp_var]
dat5$imp = 5


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
index_p = imputed_data$AKI_Category == 'Persistent_AKI'
index_t = imputed_data$AKI_Category == 'Transient_AKI'
index_n = imputed_data$AKI_Category == 'no_AKI'

by = 'all'
# summarize data by imp group
for (by in c('Persistent_AKI', 'Transient_AKI', 'no_AKI', 'all')){
  
  if (by == 'Transient_AKI') {
    num_data_sub = num_data[index_t,]
    cat_data_sub = cat_data[index_t,]
  } else if (by == 'Persistent_AKI') {
    num_data_sub = num_data[index_p,]
    cat_data_sub = cat_data[index_p,]
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
  write.csv(num_sum, paste('./Summary Statistics/summary_numeric_', by, '.csv', sep=''))
  write.csv(cat_sum, paste('./Summary Statistics/summary_categorical_', by, '.csv', sep=''))
}


