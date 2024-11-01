library(tidyverse)

# read the unimputed data
rawData = readRDS('G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')

#Persistent AKI stage 1: 1-1-1s + those ending in 1 (i.e., 2-2-1, or 2-1-1)
#Persistent AKI stage 2: 2-2-2s + those ending in 2 (i.e., 3-3-2, or 2-1-2)
#Persistent AKI stage 3: Any stage combination lasting 3 days, with stage 3 on day 3 (i.e., 1-2-3, or 1-1-3)
#Transient AKI = Return to no KDIGO stage within 72 hours of first AKI (i.e., 2-2-0, 2-0-0, 1-0-0, 3-3-0)

stage1_pattern = c("NANA1", "NA11", "NA21", "NA31",
                   "1NA1", "111", "121", "131", 
                   "2NA1", "211", "221", "231",
                   "3NA1", "311", "321", "331")

stage2_pattern = c("NANA2", "NA12", "NA22", "NA32",
                   "1NA2", "112", "122", "132",
                   "2NA2", "212", "222", "232",
                   "3NA2", "312", "322", "332")

stage3_pattern = c("NANA3", "NA13", "NA23", "NA33",
                   "1NA3", "113", "123", "133", 
                   "2NA3", "213", "223", "233",
                   "3NA3", "313", "323", "333")

tAKI_pattern = c("10NA", "100", "101", "102", "103", "1NA0", "110", "120", "130",
                 "20NA", "200", "201", "202", "203", "2NA0", "210", "220", "230",
                 "30NA", "300", "301", "302", "303", "3NA0", "310", "320", "330")

pattern_all = c(stage1_pattern, stage2_pattern, stage3_pattern, tAKI_pattern)

# Function to compare patterns in a given string
compare_patterns <- function(pattern, text) {
  matches = gregexpr(pattern, text)
  first_position <- min(matches[[1]])
  return(first_position)
}

result_df = tibble(int = c(1:nrow(rawData))) %>% as.data.frame()
for (pattern in pattern_all){
  result = mapply(compare_patterns, pattern = pattern, text = rawData$KDIGO_pattern)
  result = ifelse(result == -1, 999, result)
  result_df_new = tibble(position = result) %>% as.data.frame()
  colnames(result_df_new)[1] = paste(pattern)
  result_df = cbind(result_df, result_df_new)
}

result_df = result_df %>% select(-int)

get_min_position = function(row) {
  return(which.min(row))
}

get_min_value = function(row) {
  return(min(row))
}

# Apply the function to each row of the data frame
min_position_no_tAKI = apply(result_df[,1:48], 1, get_min_position)
min_position_tAKI = apply(result_df[,48:75], 1, get_min_position) + 27
min_value_no_tAKI = apply(result_df[,1:48], 1, get_min_value)
min_value_tAKI = apply(result_df[,48:75], 1, get_min_value)
min_position_df = tibble(min_position_no_tAKI=min_position_no_tAKI, min_position_tAKI=min_position_tAKI, 
                         min_value_no_tAKI = min_value_no_tAKI, min_value_tAKI = min_value_tAKI)
min_position_df = min_position_df %>%
  mutate(min_value_tAKI_true = ifelse(min_value_tAKI>1, 999, 1)) %>%
  select(-(min_value_tAKI))

min_position_df = min_position_df %>%
  mutate(min_position = ifelse(min_value_tAKI_true < min_value_no_tAKI, min_position_tAKI,
                               ifelse(min_value_no_tAKI < min_value_tAKI_true, min_position_no_tAKI, 999)))

min_position_df = min_position_df %>%
  mutate(subgroup_new = ifelse((min_position<=16), "Persistent AKI stage 1", 
                               ifelse((min_position>16 & min_position<=32), "Persistent AKI stage 2",
                                      ifelse((min_position>32 & min_position<=48), "Persistent AKI stage 3",
                                             ifelse((min_position>48 & min_position<999), "Transient_AKI", NA)))))

rawData$subgroup_new = min_position_df$subgroup_new

rawData = rawData %>%
  mutate(AKI_Category_subgroup_new = ifelse((AKI_Category_subgroup == "Persistent_Severe_AKI" | AKI_Category_subgroup == "no_AKI"), 
                                            as.character(AKI_Category_subgroup),
                                            subgroup_new))

check = rawData %>%
  select(patientid, patientvisitid, KDIGO_pattern, subgroup_new, AKI_Category_subgroup_new, AKI_Category_subgroup) %>%
  filter(is.na(AKI_Category_subgroup_new))

saveRDS(rawData, 'G:/Persistent_AKI/Xinlei/Data/Data_Revised.rds')


