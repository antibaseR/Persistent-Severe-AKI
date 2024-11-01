library(tidyverse)
library(survminer)
library(survival)
library(lubridate)
library(ggalluvial)

analysis = "_subgroup_trajectory_newsubgroup"

impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

# create list to save the model fit
pasta_plot_list = list()


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
surv_var = c("age_at_admission",
             "gender",
             "race",
             "charlson",
             "vp_d",
             "apache3",
             "sofa_24",
             "SEPSIS",
             "CARDIAC_SURG",
             "HEART_TRANSP",
             "LUNG_TRANSP",
             "LIVER_TRANSP",
             'MAX_LACTATE',
             "MV_DAYS",
             "rrt_discharge",
             "weight",
             "BMI",
             "TOTAL_INPUTS_BEFORE_AKI",
             "TOTAL_OUTPUTS_BEFORE_AKI",
             "dead",
             "death_date_sh",
             "encounter_start_date_time_sh",
             "encounter_end_date_time_sh",
             "dt_icu_start_sh",
             "FIRST_AKI_DATETIME",
             "AKI_Category",
             "AKI_Category_true",
             "AKI_Category_subgroup_new",
             "CKD")


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
  if (analysis == "_subgroup_trajectory_newsubgroup"){
    data = data %>%
      filter(ESRD==0 & ecmo==0 & BASELINE_eGFR>=15 & refCreat<4 & CKD_stage_4==0 & CKD_stage_5==0)
    data = data[order(data$encounter_start_date_time_sh),] #order dataframe by encounter time
    #subset to rows that are not duplicates of a previous encounter for that patient
    data = data[!duplicated(data$patientid),]
  }
  
  # make the alluvial plot
  data_n = data %>%
    select(day1_kidgo, day2_kidgo, day3_kidgo, day4_kidgo,
           day5_kidgo, day6_kidgo, day7_kidgo,
           AKI_Category_subgroup_new)
  
  ## reformat the levels of the factors
  freq_table = table(data_n) %>%
    as.data.frame() %>%
    filter(Freq!=0) %>%
    mutate(AKI_Category_subgroup_new = factor(AKI_Category_subgroup_new, levels = c("no_AKI",
                                                                                    "Transient_AKI",
                                                                                    "Persistent AKI stage 1",
                                                                                    "Persistent AKI stage 2",
                                                                                    "Persistent AKI stage 3",
                                                                                    "Persistent_Severe_AKI")))
  is_alluvia_form(as.data.frame(freq_table), axes = 1:8, silent = TRUE)
  
  
  freq_table = freq_table %>%
    mutate(day1_kidgo = factor(day1_kidgo, levels=paste(c(0,1,2,3))),
           day2_kidgo = factor(day2_kidgo, levels=paste(c(0,1,2,3))),
           day3_kidgo = factor(day3_kidgo, levels=paste(c(0,1,2,3))),
           day4_kidgo = factor(day4_kidgo, levels=paste(c(0,1,2,3))),
           day5_kidgo = factor(day5_kidgo, levels=paste(c(0,1,2,3))),
           day6_kidgo = factor(day6_kidgo, levels=paste(c(0,1,2,3))),
           day7_kidgo = factor(day7_kidgo, levels=paste(c(0,1,2,3))))
  
  freq_table %>%
    #ggplot(aes(y = Freq, 
    #           axis1 = day1_kidgo, axis2 = day2_kidgo, axis3 = day3_kidgo, axis4 = day4_kidgo,
    #           axis5 = day5_kidgo, axis6 = day6_kidgo, axis7 = day7_kidgo)) +
    ggplot(aes(y = Freq, 
               axis1 = day1_kidgo, axis2 = day2_kidgo, axis3 = day3_kidgo)) +
    geom_alluvium(aes(fill=AKI_Category_subgroup_new), width = 1/12) +
    geom_stratum(reverse = T) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), reverse = T) +
    theme(legend.position="bottom") +
    scale_fill_manual(values = c("#30c68b",
                                 "#5dcfd6",
                                 "#ffd92f",
                                 "#ff9a2f",
                                 "#c26666"),
                      labels = c("Transient AKI",
                                 "Persistent AKI stage 1",
                                 "Persistent AKI stage 2",
                                 "Persistent AKI stage 3",
                                 "Persistent Severe AKI")) +
    scale_x_discrete(limits = paste0("Day ", seq(1,3,1)), expand = c(.05, .05)) +
    labs(x = "Worse daily KDIGO stages after AKI diagnosis", y = "Frequency", fill = "AKI Category") +
    theme_bw()
  
  ggsave(filename = paste("G:/Persistent_AKI/Xinlei/Result/Pasta_Plot/pasta_plot_", i, analysis, ".png", sep=""), width = 10, height = 8)
  
}


