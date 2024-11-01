# create renal recovery data
# - this will be created using the imputed reference creatinine data, so you 
#   will need that first otherwise there will be 11,631 patientvisitids missing 
#   refCreat values
# - Recovery from AKI was assessed based on reference creatinine (SCr) at hospital
#   discharge, which was defined as a return of SCr to <= 150% of baseline 
#   without the need for dialysis. 
# - Non-recovery AKI was defined as need for dialysis within 48h prior to 
#   discharge or no return of SCr to <= 150% of baseline.

# libraries needed
library(tidyverse)
library(DBI)
options(scipen=999)

# load the imputed data
impDat = readRDS("G:/Persistent_AKI/Xinlei/Data/Data_Imputed.rds")

# get the paitentvisitid for the kidney transplant visits
kidney_transp_data = read.csv("G:/Persistent_AKI/revised_pAKI_Final_Data_v2.csv") %>%
  filter(KIDNEY_TRANSP==0)
no_kidney_transp_encounter = kidney_transp_data$patientvisitid

final_new_data6 = kidney_transp_data

for (i in 1:length(impDat)){
  
  imp = impDat[[i]]
  
  # get the patientid and patientvisitid
  patientid_data = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds') %>%
    select(patientid, patientvisitid)
  
  # combine the id to the imputed data
  imp = cbind(patientid_data, imp) %>%
    arrange(patientid, patientvisitid)
  
  # Collect data from the Server
  con <- dbConnect(odbc::odbc(), "CCM_EHR", timeout = 10)
  
  step2_Wide_File <- tbl(con, "vw_Persistent_AKI_step2_Output_Wide_File") %>% 
    collect()
  step3_Wide_File <- tbl(con, "vw_Persistent_AKI_step3_Output_Wide_File") %>% 
    collect()
  step3_Cr_Long_File <- tbl(con, "vw_Persistent_AKI_step3_Output_Cr_Long_File") %>% 
    collect()
  step3_RRT_Long_File <- tbl(con, "vw_Persistent_AKI_step3_Output_RRT_Long_File") %>% 
    collect()
  
  step2_Wide_File = step2_Wide_File %>%
    filter(patientvisitid %in% no_kidney_transp_encounter) %>%
    arrange(patientid, patientvisitid)
  step3_Wide_File = step3_Wide_File %>%
    filter(patientvisitid %in% no_kidney_transp_encounter) %>%
    arrange(patientid, patientvisitid)
  step3_Cr_Long_File = step3_Cr_Long_File %>%
    filter(patientvisitid %in% no_kidney_transp_encounter) %>%
    arrange(patientid, patientvisitid)
  step3_RRT_Long_File = step3_RRT_Long_File %>%
    filter(patientvisitid %in% no_kidney_transp_encounter) %>%
    arrange(patientid, patientvisitid)
  
  # this grabs the patientid, patientvisitid and reference creatinine, which 
  # we will use later.
  # Xinlei, this where we need to get the imputed refCreat amounts
  
  ## add the imputed refCreat to the original dataset
  imputed_creat_data = imp %>%
    select(patientid, patientvisitid, refCreat)
  colnames(imputed_creat_data)[3] = "refCreat_imp"
  step3_Wide_File = merge(step3_Wide_File, imputed_creat_data, by=c("patientid", "patientvisitid"), all.x = T)
  # mean(is.na(step3_Wide_File$refCreat_imp)) 0 missingness
  
  ## replace the unimputed one with the imputed one and drop the added one from the dataset
  step3_Wide_File = step3_Wide_File %>%
    mutate(refCreat = refCreat_imp) %>%
    select(-refCreat_imp)
  
  AKI_baseline_creatinine <- step3_Wide_File %>% 
    select(patientid, patientvisitid,
           refCreat)  
  
  AKI_baseline_creatinine %>% 
    filter(is.na(refCreat))
  
  
  # find creatinine at discharge
  # - Xinlei will run this code with the imputed data set so we can use the
  #   imputed reference creatinine values 
  ################################################################################
  # this was copied from the RRT_Script.R file 
  
  # Hospital and ICU_admission times
  ICU_admission_dt <- step2_Wide_File %>% 
    distinct(patientid, patientvisitid,
             encounter_start_date_time_sh, encounter_end_date_time_sh,
             dt_icu_start_sh, dt_icu_end_sh)
  
  # these patientvisitids have 2 or more admission dates, so we will drop them
  # because we don't know which is correct.
  dbls <- ICU_admission_dt %>% 
    count(patientvisitid) %>% 
    filter(n>1) %>% 
    pull(patientvisitid)
  
  # the cleaned hospital and ICU admission times (doubles are dropped)
  ICU_admission_dt_clean <- ICU_admission_dt %>% 
    filter(!patientvisitid %in% dbls)
  
  # clean up the daily rrt data
  RRT_long <- step3_RRT_Long_File %>% 
    relocate(patientid, .before = patientvisitid) %>% 
    arrange(patientid, patientvisitid, event_Date_sh)
  RRT_long
  
  rrt <- RRT_long %>% 
    mutate(rrt_date = as_date(event_Date_sh)) %>% 
    distinct(patientid, patientvisitid, rrt_date, kdigo, criteria) %>% 
    arrange(patientid, patientvisitid, rrt_date)
  
  check = step3_Wide_File %>%
    select(patientid, patientvisitid, refCreat)
  
  check2 = merge(step3_Cr_Long_File, check, by=c("patientid", "patientvisitid"), all.x = T)
  
  
  # clean up the daily creatinine data
  cr <- step3_Cr_Long_File %>% 
    mutate(cr_date = as_date(event_date_sh)) %>% 
    group_by(patientvisitid, cr_date) %>% 
    mutate(DAILY_MAX_CR = max(rollup_val)) %>% 
    relocate(DAILY_MAX_CR, .before = rollup_unit) %>% 
    relocate(cr_date, .before = rollup_name) %>% 
    ungroup() %>% 
    select(-c(event_date_sh, rollup_val)) %>% 
    distinct() %>% 
    arrange(patientid, patientvisitid, cr_date)
  
  
  # join the admission data, rrt, creatinine and reference creatinine data
  ICU_rrt_cr_data <- ICU_admission_dt_clean %>% 
    left_join(rrt, by = c("patientid", "patientvisitid"), multiple="all") %>% 
    left_join(cr, by = c("patientid", "patientvisitid"), multiple="all") %>% 
    left_join(AKI_baseline_creatinine, by = c("patientid", "patientvisitid")) %>% 
    select(patientid, patientvisitid,
           encounter_end_date_time_sh, dt_icu_start_sh, 
           rrt_date,
           cr_date, 
           DAILY_MAX_CR,
           refCreat) %>% 
    distinct()
  
  
  # creates the rrt_discharge data... from the RRT_Script.R file
  # this is how we defined RRT at discharge
  rrt_discharge <- ICU_rrt_cr_data %>% 
    select(-dt_icu_start_sh) %>%
    mutate(minus7 = as_date(encounter_end_date_time_sh) - days(7),
           plus7 = as_date(encounter_end_date_time_sh) + days(7)) %>% 
    mutate(rrt_discharge = if_else(rrt_date >= minus7 & rrt_date <= plus7, 1, 0)) %>% 
    mutate(rrt_discharge = if_else(is.na(rrt_discharge), 0, rrt_discharge)) %>% 
    mutate(cr_pm7 = if_else(cr_date >= minus7 & cr_date <= plus7, 1, 0)) %>% 
    mutate(CR_over200 = if_else(rrt_discharge == 0 & cr_pm7 == 1 & DAILY_MAX_CR >= 2*refCreat, 1, 0)) %>% 
    mutate(CR_over200 = if_else(is.na(CR_over200), 0, CR_over200)) %>% 
    mutate(rrt_discharge = if_else(rrt_discharge == 0 & CR_over200 == 1, 1, rrt_discharge))
  
  rrt_discharge_data <- rrt_discharge %>% 
    select(patientid, patientvisitid, 
           CR_over200, rrt_discharge) %>% 
    group_by(patientvisitid) %>% 
    mutate(CR_over200 = max(CR_over200),
           rrt_discharge = max(rrt_discharge)) %>% 
    ungroup() %>% 
    distinct()
  
  
  rrt_discharge_data # looks good!
  
  ################################################################################
  ICU_rrt_cr_data # this comes from the RRT_script.R file (see above if needed)
  
  # we want either the max creatinine level at discharge or the last creatinine 
  # level taken and then to see if the received RRT at discharge.
  # If the last creatinine level is less than or equal to 150% of the reference 
  # creatinine an did not have RRT at discharge, then they had renal recovery
  cr_discharge_data <- ICU_rrt_cr_data %>% 
    mutate(discharge_date = as_date(encounter_end_date_time_sh)) %>% 
    relocate(discharge_date, .before = rrt_date) %>% 
    select(-rrt_date) %>% 
    #select(-c(encounter_end_date_time_sh, dt_icu_start_sh, rrt_date)) %>% 
    #mutate(CR_at_discharge = if_else(cr_date == discharge_date, 1, 0)) %>% 
    mutate(CR_lt_150 = if_else(DAILY_MAX_CR <= 1.5*refCreat, 1, 0))   # identifies if creatinine is less than 150% of the reference creatinine
  cr_discharge_data
  
  cr_discharge_last_data <- cr_discharge_data %>% 
    group_by(patientvisitid) %>% 
    slice_tail(n=1) %>%  # this looks good but gives us some NA values because patients don't have DAILY_MAX_CR or refCreat... what do we want to do with these people???
    ungroup()
  
  cr_discharge_last_data
  rrt_discharge_data  # comes from the RRT_Script.R file
  
  renal_recovery_data <- cr_discharge_last_data %>% 
    left_join(rrt_discharge_data) %>% 
    select(-CR_over200) %>% 
    mutate(renal_recov = if_else(CR_lt_150 == 1 & rrt_discharge == 0, 1, 0))
  
  renal_recovery_data  # this looks good
  
  
  # clean up the data and only keep keys and what we want to join
  renal_recovery_data_clean <- renal_recovery_data %>% 
    select(patientid, patientvisitid, refCreat, rrt_discharge, renal_recov) %>%
    mutate(patientid = as.numeric(patientid),
           patientvisitid = as.numeric(patientvisitid))
  
  
  # join the renal recovery data with the data to create the final dataset
  final_new_data7 = merge(final_new_data6, renal_recovery_data_clean, by=c("patientid", "patientvisitid"), all.x=T)
  final_new_data7 = final_new_data7 %>%
    select(patientid, patientvisitid, renal_recov)
  
  # save the final_new_data7 set as the revised_pAKI_Final_Data_v2.csv file
  write_csv(final_new_data7, paste("renal_recovery_data_imputed_", i, ".csv", sep=""))
  
}




