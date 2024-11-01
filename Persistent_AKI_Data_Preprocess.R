library(tidyverse)
library(mice)
library(corrplot)
library(table1)
library(missRanger)
#library(shapr)
options(scipen = 999)


# load data, exculde kidney transplant patients from the start
setwd("G:/Persistent_AKI")
rawdata = read.csv("./Persistent_AKI_Data/revised_pAKI_Final_Data_v2.csv")
data = rawdata %>%
  filter(KIDNEY_TRANSP==0)
sgData = read.csv("G:/Persistent_AKI/Persistent_AKI_Data/AKI_Subgroups.csv")
sgData_new = sgData %>%
  select(-c("dt_icu_start_sh", "dt_icu_end_sh", "AKI_Category", "aki_2_3_icu", "aki_2_3_DATE", "aki_1_icu", "aki_1_DATE"))
data = merge(data, sgData_new, by=c("patientid", "patientvisitid"), all.x=T)

# missingness of each variable
missing = data.frame(colSums(is.na(data)) / nrow(data))
colnames(missing) = c("missingness")

# replace the NA with 0
replace_var = c('TOTAL_BLOOD_PROD_BEFORE_AKI',
                'TOTAL_ALBUMIN_BEFORE_AKI',
                'TOTAL_DOSE_HOSP_norepinephrine',
                'AVERAGE_DAILY_DOSE_HOSP_norepinephrine',
                'TOTAL_DOSE_BEFORE_AKI_norepinephrine',
                'AVERAGE_DAILY_DOSE_BEFORE_AKI_norepinephrine',
                'TOTAL_DOSE_HOSP_phenylephrine',
                'AVERAGE_DAILY_DOSE_HOSP_phenylephrine',
                'TOTAL_DOSE_BEFORE_AKI_phenylephrine',
                'AVERAGE_DAILY_DOSE_BEFORE_AKI_phenylephrine',
                'TOTAL_DOSE_HOSP_vasopressin',
                'AVERAGE_DAILY_DOSE_HOSP_vasopressin',
                'TOTAL_DOSE_BEFORE_AKI_vasopressin',
                'AVERAGE_DAILY_DOSE_BEFORE_AKI_vasopressin',
                'TOTAL_DOSE_HOSP_dopamine',
                'AVERAGE_DAILY_DOSE_HOSP_dopamine',
                'TOTAL_DOSE_BEFORE_AKI_dopamine',
                'AVERAGE_DAILY_DOSE_BEFORE_AKI_dopamine',
                'TOTAL_DOSE_HOSP_epinephrine',
                'AVERAGE_DAILY_DOSE_HOSP_epinephrine',
                'TOTAL_DOSE_BEFORE_AKI_epinephrine',
                'AVERAGE_DAILY_DOSE_BEFORE_AKI_epinephrine',
                'TOTAL_DOSE_HOSP_dobutamine',
                'AVERAGE_DAILY_DOSE_HOSP_dobutamine',
                'TOTAL_DOSE_BEFORE_AKI_dobutamine',
                'AVERAGE_DAILY_DOSE_BEFORE_AKI_dobutamine',
                'TOTAL_DOSE_HOSP_milrinone',
                'AVERAGE_DAILY_DOSE_HOSP_milrinone',
                'TOTAL_DOSE_BEFORE_AKI_milrinone',
                'AVERAGE_DAILY_DOSE_BEFORE_AKI_milrinone',
                'mv_d',
                'vp_d',
                'MAX_PEEP',
                'MAX_BG_TIDAL_VOL',
                'MAX_TIDAL_VOL_RC',
                "aki_1_icu",
                "aki_2_3_icu",
                "ecmo")


data = data %>% 
  mutate_at(replace_var, ~replace_na(.,0))

# check if all the NAs replaced by 0
rep_data = data[,replace_var]
missing_rep = data.frame(colSums(is.na(rep_data)) / nrow(rep_data))
colnames(missing_rep) = c("missingness")

# create a composite variable of all CKD/ESRD
CKD_data = data %>%
  select("CKD",
         "CKD_stage_1",
         "CKD_stage_2",
         "CKD_stage_3",
         "CKD_stage_4",
         "CKD_stage_5",
         "ESRD",
         "CKD_unspecified",
         "diabetic_CKD")
missing_ckd = data.frame(colSums(is.na(CKD_data)) / nrow(CKD_data))
colnames(missing_ckd) = c("missingness")
CKD_data$sum = rowSums(CKD_data)
CKD_data$ALL_CKD_ESRD = ifelse(CKD_data$sum==0, 0, 1)
data$ALL_CKD_ESRD = CKD_data$ALL_CKD_ESRD

# transform categorical variables into factors
unique_value = apply(data, 2, function(x) length(unique(x)))
cat_cols = unique_value<=4
cat_var = data[, cat_cols] %>% colnames()
cat_var = cat_var[cat_var!="MIN_SCVO2"]
data[, cat_var] = lapply(data[, cat_var], factor)
#str(data)

data$AKI_Category = as.character(data$AKI_Category)
data$AKI_Category_true = ifelse((data$aki_1_icu==1 & data$AKI_Category=="no_AKI"), "stage_1_AKI", data$AKI_Category)
data$AKI_Category = factor(data$AKI_Category, 
                           levels = c("no_AKI", "Transient_AKI", "Persistent_AKI"))
data$AKI_Category_true = factor(data$AKI_Category_true, 
                                levels = c("no_AKI", "stage_1_AKI", "Transient_AKI", "Persistent_AKI"))

# create a variable for the KDIGO sequence
data$KDIGO_pattern = paste0(data$day1_kidgo, data$day2_kidgo, data$day3_kidgo, data$day4_kidgo, 
                            data$day5_kidgo, data$day6_kidgo, data$day7_kidgo)

data$KDIGO_pattern_first3 = paste0(data$day1_kidgo, data$day2_kidgo, data$day3_kidgo)

# remove weird pattern subjects
#data = data %>%
#  filter(!(data$KDIGO_pattern=="NANANANANANANA" & data$AKI_Category_true=="Transient_AKI"))

data$KDIGO333 = ifelse(grepl("333", data$KDIGO_pattern), TRUE, FALSE)

# create the new labels for the AKI category subgroups
data = data %>%
  mutate(AKI_Category_subgroup = ifelse((is.na(subgroup)|AKI_Category_true == "stage_1_AKI"), as.character(AKI_Category_true), subgroup)) %>%
  mutate(AKI_Category_subgroup = ifelse((AKI_Category_true=="Transient_AKI" & KDIGO333==TRUE), 
                                        "Persistent_Severe_AKI", AKI_Category_subgroup)) %>%
  mutate(AKI_Category_subgroup = ifelse(AKI_Category_subgroup=="Persistent_AKI", "Persistent_Severe_AKI", AKI_Category_subgroup))

data$AKI_Category_subgroup = ifelse(data$AKI_Category_subgroup %in% c("no_AKI", "stage_1_AKI", "Persistent_Severe_AKI"), 
                                    data$AKI_Category_subgroup, "No_Persistent_Severe_AKI")
data$AKI_Category_subgroup = factor(data$AKI_Category_subgroup, levels = c("no_AKI", 
                                                                           "stage_1_AKI",
                                                                           "No_Persistent_Severe_AKI",
                                                                           "Persistent_Severe_AKI"))

data$AKI_Category_subgroup2 = ifelse((data$AKI_Category_subgroup_new %in% c("Persistent AKI stage 1", "Persistent AKI stage 2", "Persistent AKI stage 3")), 
                                     "Persistent_Mild_Moderate_AKI", data$AKI_Category_subgroup_new)
table(data$AKI_Category_subgroup2)
data$AKI_Category_subgroup2 = factor(data$AKI_Category_subgroup2, levels = c("no_AKI", 
                                                                             "Transient_AKI",
                                                                             "Persistent_Mild_Moderate_AKI",
                                                                             "Persistent_Severe_AKI"))



# change the working directory to my folder
setwd("G:/Persistent_AKI/Xinlei")
saveRDS(data, './Data/Data_Revised.rds')