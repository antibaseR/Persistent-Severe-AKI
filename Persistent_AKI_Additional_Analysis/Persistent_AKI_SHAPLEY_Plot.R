library(tidyverse)
library(caret)
library(xgboost)
library(caretEnsemble)
library(pROC)
library(caTools)
library(randomForest)
library(nnet)
library(DALEX)
library(iml)

for (analysis in c("_BaseCreat", "_BaseCreat_sensitivity",
                   "_sepsis", "_sepsis_sensitivity",
                   "_no_sepsis_no_cardiac_surgery", "_no_sepsis_no_cardiac_surgery_sensitivity",
                   "_cardiac_surgery", "_cardiac_surgery_sensitivity"
                   )){
  
  # create list to save the results
  test_sum_list = list()
  shap_val_list = list()
  
  # run only after having the SHAPLEY values
  shap_val_list = list()
  shap_val_list[[1]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_1", analysis, ".rds", sep=""))
  shap_val_list[[2]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_2", analysis, ".rds", sep=""))
  shap_val_list[[3]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_3", analysis, ".rds", sep=""))
  shap_val_list[[4]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_4", analysis, ".rds", sep=""))
  shap_val_list[[5]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_5", analysis, ".rds", sep=""))
  saveRDS(shap_val_list, paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/shap_val_list", analysis, ".rds", sep=""))
  
  
  
  # read the result
  library(readxl)
  library(tidyverse)
  
  test_sum_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/test_sum_list", analysis, ".rds", sep=""))
  shap_val_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/shap_val_list", analysis, ".rds", sep=""))
  
  # the testing AUC
  test_sum_df = data.frame(glmnet=as.numeric(),
                           xgb = as.numeric(),
                           ensemble = as.numeric())
  for (i in 1:length(test_sum_list)){
    test_sum_df[i,] = round(test_sum_list[[i]], 3)
  }
  
  
  # SHAPLEY values
  shap_val_df = data.frame()
  for (i in 1:length(shap_val_list)){
    shap_val_df_new = shap_val_list[[i]]
    shap_val_df_new$imp = i
    shap_val_df = rbind(shap_val_df, shap_val_df_new)
  }
  shap_val_df = shap_val_df %>%
    mutate(value = as.numeric(sub(".*=", "", feature.value))) %>%
    filter(class=="Yes")
  
  ## get the top features
  shap_top_features = shap_val_df %>%
    group_by(feature) %>%
    summarise(shap_feature = mean(abs(phi))) %>%
    arrange(desc(shap_feature)) %>%
    filter(shap_feature!=0)
  top_features_order_candidate = shap_top_features$feature
  
  shap_selected_features = shap_val_df %>%
    filter(feature %in% top_features_order_candidate) %>%
    group_by(feature) %>%
    summarise(shap_feature_max = max(abs(phi)),
              shap_feature = mean(abs(phi))) %>%
    arrange(desc(shap_feature)) %>%
    filter(shap_feature_max>0.2)
  
  top_features_order = shap_selected_features$feature
  shap_val_df$feature = factor(shap_val_df$feature, levels=rev(top_features_order))
  
  ## take a subset of the data
  shap_val_df_sub = shap_val_df %>%
    filter(feature %in% top_features_order) %>%
    group_by(feature) %>%
    mutate(rescaled_value = (value - min(value)) / (max(value) - min(value)))
  
  ## calculate the global shap value
  shap_val_df_sub_global = shap_val_df_sub %>%
    group_by(feature) %>%
    summarize(mean_ab_value_phi = mean(abs(phi)))
  
  ## get the label of each variable
  feature_dic = read_excel("G:/Persistent_AKI/Xinlei/Data/feature_names.xlsx")
  feature_label = feature_dic[feature_dic$`Name of variable in the data set` %in% top_features_order,] %>%
    mutate(`Name of variable in the data set` = factor(`Name of variable in the data set`, level = rev(top_features_order))) %>%
    arrange(`Name of variable in the data set`) %>%
    select(Variable)
  
  ## make the SHAPLEY plot
  shap_val_df_sub %>%
    ggplot(aes(x = phi, y = feature, color = rescaled_value)) +
    geom_point(size=2, alpha=0.5) +
    labs(x = "SHAP Values", y = "Feature Names", color = "Feature Values") +
    scale_color_gradient(low="#ffcc31", high="#5c1ca4", breaks=c(0,1), label = c("Low", "High")) +
    scale_x_continuous(breaks = seq(-1,1,0.2)) +
    scale_y_discrete(labels = feature_label$Variable) +
    theme_classic()
  
  ggsave(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/shap_plot", analysis, ".png", sep=""), height=5, width=10)
  
  ## make the SHAPLEY bar plot
  shap_val_df_sub_global %>%
    ggplot(aes(x = mean_ab_value_phi, y = feature)) +
    geom_bar(stat="identity", fill = "blue") +
    labs(x = "SHAP Values", y = "Feature Names", color = "Feature Values") +
    scale_x_continuous(breaks = seq(0,0.1,0.02)) +
    scale_y_discrete(labels = feature_label$Variable) +
    theme_classic()
  
  ggsave(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/shap_bar_plot", analysis, ".png", sep=""), height=5, width=10)
  
  
  # AUC
  AUC_df = data.frame()
  for (i in 1:length(test_sum_list)){
    AUC_df_new = data.frame(test_sum_list[[i]])
    AUC_df = rbind(AUC_df, AUC_df_new)
  }
  colnames(AUC_df) = c("LASSO Regression", "XGBoost", "Super Learner")
  write.csv(AUC_df, paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/AUC_df", analysis, ".csv", sep=""))
}




for (analysis in c("_ageq", "_ageq_sensitivity")){
  
  # create list to save the results
  test_sum_list = list()
  shap_val_list = list()
  
  # create list for the outer loop
  test_sum_list_imp = list()
  shap_val_list_imp = list()
  
  # read the result
  library(readxl)
  library(tidyverse)
  
  test_sum_list = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/test_sum_list", analysis, ".rds", sep=""))
  
  for (j in 1:4){
    
    shap_val_list = list()
    shap_val_list[[1]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_1_", j, ".rds", sep=""))
    shap_val_list[[2]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_2_", j, ".rds", sep=""))
    shap_val_list[[3]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_3_", j, ".rds", sep=""))
    shap_val_list[[4]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_4_", j, ".rds", sep=""))
    shap_val_list[[5]] = readRDS(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/SHAP_df_5_", j, ".rds", sep=""))
    saveRDS(shap_val_list, paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Data/shap_val_list", analysis, "_", j, ".rds", sep=""))
    
    
    # the testing AUC
    test_sum_df = data.frame(glmnet=as.numeric(),
                             xgb = as.numeric(),
                             ensemble = as.numeric())
    for (i in 1:length(test_sum_list)){
      test_sum_df[i,] = round(test_sum_list[[i]][[j]], 3)
    }
    
    
    # SHAPLEY values
    shap_val_df = data.frame()
    for (i in 1:length(shap_val_list)){
      shap_val_df_new = shap_val_list[[i]]
      shap_val_df_new$imp = i
      shap_val_df = rbind(shap_val_df, shap_val_df_new)
    }
    shap_val_df = shap_val_df %>%
      mutate(value = as.numeric(sub(".*=", "", feature.value))) %>%
      filter(class=="Yes")
    
    ## get the top features
    shap_top_features = shap_val_df %>%
      group_by(feature) %>%
      summarise(shap_feature = mean(abs(phi))) %>%
      arrange(desc(shap_feature)) %>%
      filter(shap_feature!=0)
    top_features_order_candidate = shap_top_features$feature
    
    shap_selected_features = shap_val_df %>%
      filter(feature %in% top_features_order_candidate) %>%
      group_by(feature) %>%
      summarise(shap_feature_max = max(abs(phi)),
                shap_feature = mean(abs(phi))) %>%
      arrange(desc(shap_feature)) %>%
      filter(shap_feature_max>0.2)
    
    top_features_order = shap_selected_features$feature
    shap_val_df$feature = factor(shap_val_df$feature, levels=rev(top_features_order))
    
    ## take a subset of the data
    shap_val_df_sub = shap_val_df %>%
      filter(feature %in% top_features_order) %>%
      group_by(feature) %>%
      mutate(rescaled_value = (value - min(value)) / (max(value) - min(value)))
    
    ## calculate the global shap value
    shap_val_df_sub_global = shap_val_df_sub %>%
      group_by(feature) %>%
      summarize(mean_ab_value_phi = mean(abs(phi)))
    
    ## get the label of each variable
    feature_dic = read_excel("G:/Persistent_AKI/Xinlei/Data/feature_names.xlsx")
    feature_label = feature_dic[feature_dic$`Name of variable in the data set` %in% top_features_order,] %>%
      mutate(`Name of variable in the data set` = factor(`Name of variable in the data set`, level = rev(top_features_order))) %>%
      arrange(`Name of variable in the data set`) %>%
      select(Variable)
    
    ## make the SHAPLEY plot
    shap_val_df_sub %>%
      ggplot(aes(x = phi, y = feature, color = rescaled_value)) +
      geom_point(size=2, alpha=0.5) +
      labs(x = "SHAP Values", y = "Feature Names", color = "Feature Values") +
      scale_color_gradient(low="#ffcc31", high="#5c1ca4", breaks=c(0,1), label = c("Low", "High")) +
      scale_x_continuous(breaks=seq(-1,1,0.2)) +
      scale_y_discrete(labels = feature_label$Variable) +
      theme_classic()
    
    ggsave(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/shap_plot", analysis, "_", j, ".png", sep=""), height=5, width=10)
    
    ## make the SHAPLEY bar plot
    shap_val_df_sub_global %>%
      ggplot(aes(x = mean_ab_value_phi, y = feature)) +
      geom_bar(stat="identity", fill = "blue") +
      labs(x = "SHAP Values", y = "Feature Names", color = "Feature Values") +
      scale_x_continuous(breaks = seq(0,0.1,0.02)) +
      scale_y_discrete(labels = feature_label$Variable) +
      theme_classic()
    
    ggsave(paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/shap_bar_plot", analysis, "_", j, ".png", sep=""), height=5, width=10)
    
    # AUC
    AUC_df = data.frame()
    for (i in 1:length(test_sum_list)){
      AUC_df_new = data.frame(test_sum_list[[i]][[j]])
      AUC_df = rbind(AUC_df, AUC_df_new)
    }
    colnames(AUC_df) = c("LASSO Regression", "XGBoost", "Super Learner")
    write.csv(AUC_df, paste("G:/Persistent_AKI/Xinlei/Persistent_AKI_Additional_Analysis/Result/AUC_df", analysis, "_", j, ".csv", sep=""))
  }
}

