library(tidyverse)
library(mice)
library(corrplot)
library(table1)
library(missRanger)
#library(shapr)
options(scipen = 999)

# load data, exculde kidney transplant patients from the start
setwd("G:/Persistent_AKI")

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



# read the saved dataset for imputation
data = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_for_Imputation.rds")
imp_data = data[, names(data) %in% imp_var]

# read the data dictionary
library(readxl)
feature_dic = read_excel("G:/Persistent_AKI/Xinlei/Data/feature_names.xlsx")
names(feature_dic)[2] = "variable_name"

# impute the data
## check the missingness of the dataset to be imputed
missing_imp = data.frame(colSums(is.na(imp_data)) / nrow(imp_data)*100)
missing_imp$variable_name = rownames(missing_imp)
missing_df = merge(missing_imp, feature_dic, by="variable_name")
missing_df = missing_df %>%
  select(-variable_name) %>%
  select(Variable, everything())
colnames(missing_df) = c("Variable", "Missingness")
colnames(missing_df) = c("Variable", "Missingness (%)")
write.csv(missing_df, "G:/Persistent_AKI/Xinlei/Result/Missing_Table_for_Imputed_Variables.csv")







