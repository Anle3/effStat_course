library(ggplot2)
library(dplyr)
library(randomForestSRC)
library(randomForest)
library(gridExtra) 

num_simulations <- 200  
n <- 200
p <- 20
set.seed(123)

# This function used in all the scenarios in common. 
individual_trt_effect <- function(data, treat , time, status) {
  
  data_trt0 <- base::subset(data, treat == 0)
  data_trt1 <- base::subset(data, treat == 1)

   # Fitting the RandomForests
  rf_model_trt0 <- randomForestSRC::rfsrc(Surv(time,status) ~., 
                                          data_trt0,
                                          ntree = 1000
  )
  
  rf_model_trt1 <- randomForestSRC::rfsrc(Surv(time,status) ~., 
                                          data_trt1,
                                          ntree = 1000
  )
  
  # THIS MIGHT NOT BE TRUE: in survival case, direct substraction might not be the case. 
  twin0 <- predict.rfsrc(rf_model_trt0, data)$predicted
  twin1 <- predict.rfsrc(rf_model_trt1, data)$predicted
  
  # Individual treatment effect 
  z <- twin1 - twin0
  return(z = z)
}






################################ Scenario  A-1 ############################## 

Xprog <- c("bm01","bm02","bm03","bm04")        # Set of predictive biomarkers
Xpred <- c("bm04","bm05","bm06")               # Set of prognostic biomarkers
q <- length(Xpred)                             # number of the predictive biomarkers
XsolelyProg <- dplyr::setdiff(Xprog, Xpred)    # set differentiation 
Xprogpred <- dplyr::intersect(Xpred,Xprog)

# Function to run the simulation
run_VTsimulation <- function(num_simulations) {
  # Storage for results
  VT_simulation_results <- list()
  
  # Set seed for reproducibility
  set.seed(123)
  
  # Loop through simulations
  for (i in 1:num_simulations) {
    # Generate data
    
    dataBM <- biospear::simdata(
      n = n,
      p = p, 
      q.main = 4 ,                          # the number of true prognostic biomarkers.
      q.inter = 3,                          # the number of true biomarkers interacting with the treatement.
      prob.tt = 0.5,                        # the treatment assignment probability.
      m0 = 5,                               # baseline median survival time. 
      alpha.tt = -1,                        # the effect of the treatment (in log-scale).
      beta.main = c(0.2, 0.3, 0.4, 0.1),
      beta.inter = rep(-1, 3),
      b.corr = 0, 
      b.corr.by = 1, 
      wei.shape = 1, 
      recr = 1, 
      fu = 10, 
      timefactor = 1, 
      active.main = Xprog, 
      active.inter = Xpred 
    )
    
    # adjust the simulated data frame columns to start time, status, treat, bm01... and turns the trt assignment +-0.5 to 0 and 1
    dataBM <- dataBM %>%
      dplyr::relocate(time, status, treat) %>%
      dplyr::mutate(treat = ifelse(treat == "-0.5", 1, ifelse(treat == "0.5", 0, treat))) 
    
    
    # Compute individual treatment effect
    z <- individual_trt_effect(
      data = dataBM,
      treat = "treat",
      time = "time",
      status = "status"
    )
    
    # Combine data and treatment effect
    combineddata <- base::cbind(z, dataBM) %>% 
      dplyr::select(-c(time, status, treat))
    
    # Set seed for random forest reproducibility
    set.seed(123)
    # Fit random forest and compute variable importance: meaning try to understand how biomarkers important to predict differenctial trt effect z
    important_features <- randomForest(z ~ .,
                                       data = combineddata,
                                       ntree = 500,
                                       importance = TRUE)
    
    # Extract variable importance
    importance_features_ordered <- important_features %>%
                                   randomForest::importance(type = 1, scale = TRUE) %>%
                                   as.data.frame() %>%
                                   dplyr::arrange(desc(`%IncMSE`))
    
    
    
    # Store results
    VT_simulation_results[[i]] <- importance_features_ordered
  }
  
  return(VT_simulation_results)
}



################################################  ANALYSIS    ##########################################################

results_VT_A1 <- run_VTsimulation(num_simulations)
#saveRDS(results_VT_A1, "results_VT_A1.rds")



# Evaluation measures I created accordingly described in the paper Sechidis et al. 

# Extract variable names from each simulation trial. (this will give you Xpred^q)
top_q_variable_names_A1 <- lapply(results_VT_A1, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_VT_A1[[1]])


# True Positive Rate: Captures the how accurately alygorithm captures the algorithm correctly identifies the predictive biomarkers. 
VT_TPR_A1 <- lapply(top_q_variable_names_A1, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})



# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
VT_FNR_prog_A1 <- lapply(top_q_variable_names_A1, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})


# set of irrelevant features: obtained by subtracting prognostic and predictive biomarkers from all variable names
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))


VT_FNR_irrlevant_A1 <- lapply(top_q_variable_names_A1, function(fnr_irr_calculation) {
  intersection <- intersect(Xirr, fnr_irr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})


# This actually the same as 1-TPR
VT_FNR_A1 <- base::Map("+", VT_FNR_prog_A1, VT_FNR_irrlevant_A1) 


# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
VT_TPR_A1 <- base::mean(base::unlist(VT_TPR_A1))                 
VT_FNR_prog_A1 <- base::mean(base::unlist(VT_FNR_prog_A1))





# Create a list to store variable importance scores for each simulation
importance_scores_list <- lapply(all_variable_names, function(var_name) {
  sapply(results_VT_A1, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "%IncMSE"])
    } else {
      return(NA)
    }
  })
})


#### Preparation for Visualizations ####
# Combine the importance scores into a data frame
importance_scores_df <- setNames(as.data.frame(do.call(cbind, importance_scores_list)), all_variable_names)

# This helps me to format data create plots. 
# PLOT1: Variable Importance Box plot for each variable obtained from N run Simulation 

importance_score_melted_VT_A1 <- reshape2::melt(importance_scores_df) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xprogpred ~ "Xprogpred",
      variable %in% Xpred ~ "Xpred",
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  )
                                                                


importance_score_melted_VT_A1$variable <- base::sub("^bm", "", importance_score_melted_VT_A1$variable)
VT_importanceplot_A1 <- ggplot2::ggplot(importance_score_melted_VT_A1, 
                                        aes(x = reorder(variable, -value), 
                                            y = value, 
                                            fill = biomarker_types)) +
                                 ggplot2::geom_boxplot(alpha = 0.5) +
                                 ggplot2::scale_fill_manual(values = c("Xprogpred" = "darkolivegreen3",
                                                                       "Xpred" = "slateblue", 
                                                                       "Xprog" = "#EE8844",  
                                                                       "Xirr" = "grey")) +
                                 ggplot2::labs(title = sprintf("Scenario A1", 
                                                               num_simulations),
                                                               x = "Biomarker-ID",
                                                               y = "Importance",
                                                               fill = "Biomarker Type") + 
                                 ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_VT_A1$value), by = 10)) + 
                                 ggplot2::theme_bw() 
VT_importanceplot_A1



################################ Scenario  A-2 ############################## 

Xprog <- c("bm01","bm02","bm03")               # Set of predictive biomarkers
Xpred <- c("bm04","bm05","bm06")               # Set of prognostic biomarkers
q <- length(Xpred)                             # number of the predictive biomarkers
XsolelyProg <- dplyr::setdiff(Xprog, Xpred)    # set differentiation 



# Function to run the simulation
run_VTsimulation_A2 <- function(num_simulations) {
  # Storage for results
  VT_simulation_results <- list()
  
  set.seed(123)
  # Loop through simulations
  for (i in 1:num_simulations) {
  
    # Generate data
    dataBM <- biospear::simdata(
      n = n,
      p = p, 
      q.main = 3 ,                          # the number of true prognostic biomarkers.
      q.inter = 3,                          # the number of true biomarkers interacting with the treatement.
      prob.tt = 0.5,                        # the treatment assignment probability.
      m0 = 5,                               # baseline median survival time. 
      alpha.tt = -1,                        # the effect of the treatment (in log-scale).
      beta.main = c(0.2, 0.3, 0.4),
      beta.inter = rep(-1, 3),
      b.corr = 0, 
      b.corr.by = 1, 
      wei.shape = 1, 
      recr = 1, 
      fu = 10, 
      timefactor = 1, 
      active.main = Xprog, 
      active.inter = Xpred 
    )
    # adjust the simulated data frame columns to start time, status, treat, bm01... and turns the trt assignment +-0.5 to 0 and 1
    dataBM <- dataBM %>%
      dplyr::relocate(time, status, treat) %>%
      dplyr::mutate(treat = ifelse(treat == "-0.5", 1, ifelse(treat == "0.5", 0, treat))) 
    
    # Compute individual treatment effect
    z <- individual_trt_effect(
      data = dataBM,
      treat = "treat",
      time = "time",
      status = "status"
    )
    
    # Combine data and treatment effect
    combineddata <- base::cbind(z, dataBM) %>% 
      dplyr::select(-c(time, status, treat))
    
    # Set seed for random forest reproducibility
    set.seed(123)
    # Fit random forest and compute variable importance: meaning try to understand how biomarkers important to predict differenctial trt effect z
    important_features <- randomForest(z ~ .,
                                       data = combineddata,
                                       ntree = 500,
                                       importance = TRUE)
    
    # Extract variable importance
    importance_features_ordered <- important_features %>%
      randomForest::importance(type = 1, scale = TRUE) %>%
      as.data.frame() %>%
      dplyr::arrange(desc(`%IncMSE`))
    
    # Store results
    VT_simulation_results[[i]] <- importance_features_ordered
  }
  
  return(VT_simulation_results)
}


################################################  ANALYSIS    ##########################################################

results_VT_A2 <- run_VTsimulation_A2(num_simulations)
#saveRDS(results_VT_A2, "results_VT_A2.rds")


# Evaluation measures I created accordingly described in the paper Sechidis et al. 
# Extract variable names from each simulation trial. (this will give you Xpred^q)
top_q_variable_names_A2 <- lapply(results_VT_A2, function(df) rownames(df)[1:q])

# set of irrelevant features: obtained by subtracting prognostic and predictive biomarkers from all variable names
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))

# True Positive Rate: Captures the how accurately alygorithm captures the algorithm correctly identifies the predictive biomarkers. 
VT_TPR_A2 <- lapply(top_q_variable_names_A2, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})



# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
VT_FNR_prog_A2 <- lapply(top_q_variable_names_A2, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})


VT_FNR_irrlevant_A2 <- lapply(top_q_variable_names_A2, function(fnr_irr_calculation) {
  intersection <- intersect(Xirr, fnr_irr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})


# This actually the same as 1-TPR
VT_FNR_A2 <- base::Map("+", VT_FNR_prog_A2, VT_FNR_irrlevant_A2) 


# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
VT_TPR_A2 <- base::mean(base::unlist(VT_TPR_A2))                 
VT_FNR_prog_A2 <- base::mean(base::unlist(VT_FNR_prog_A2))



# Create a list to store variable importance scores for each simulation
importance_scores_list_A2 <- lapply(all_variable_names, function(var_name) {
  sapply(results_VT_A2, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "%IncMSE"])
    } else {
      return(NA)
    }
  })
})


#### Preparation for Visualizations ####
# Combine the importance scores into a data frame
importance_scores_df_A2 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_A2)), all_variable_names)

# This helps me to format data create plots. 
# PLOT1: Variable Importance Box plot for each variable obtained from N run Simulation 

importance_score_melted_VT_A2 <- reshape2::melt(importance_scores_df_A2) %>%
  dplyr::mutate(biomarker_types = dplyr::case_when(variable %in% Xpred ~ "Xpred",
                                                   variable %in% Xprog ~ "Xprog",
                                                   variable %in% Xirr ~ "Xirr"))


importance_score_melted_VT_A2$variable <- base::sub("^bm", "", importance_score_melted_VT_A2$variable)
VT_importanceplot_A2 <- ggplot2::ggplot(importance_score_melted_VT_A2, 
                                        aes(x = reorder(variable, -value), 
                                            y = value, 
                                            fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xprogpred" = "darkolivegreen3",
                                        "Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = sprintf("Scenario A2", 
                                num_simulations),
                x = "Biomarker-ID",
                y = "Importance",
                fill = "Biomarker Type") + 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_VT_A2$value), by = 10)) + 
  ggplot2::theme_bw() 
VT_importanceplot_A2



######### Setting-A: Scenario-3: 20 BM: 3-Prog, 3-Pred, 14 Irrelevant + 0.7 correlation with 4 blocks###########

Xprog <- c("bm01","bm02","bm03")               # Set of predictive biomarkers
Xpred <- c("bm04","bm05","bm06")               # Set of prognostic biomarkers
q <- length(Xpred)                             # number of the predictive biomarkers
XsolelyProg <- dplyr::setdiff(Xprog, Xpred)    # set differentiation 



# Define the simulation function
run_VTsimulation_A3 <- function(num_simulations) {
  # Define empty list to store your results 
  VT_simulation_results <- list()
  
  set.seed(123)
  # Loop through simulations
  for (i in 1:num_simulations) {
    
    # Generate data
    dataBM <- biospear::simdata(
      n = n,
      p = p, 
      q.main = 3 ,                          # the number of true prognostic biomarkers.
      q.inter = 3,                          # the number of true biomarkers interacting with the treatement.
      prob.tt = 0.5,                        # the treatment assignment probability.
      m0 = 5,                               # baseline median survival time. 
      alpha.tt = -1,                        # the effect of the treatment (in log-scale).
      beta.main = c(0.2, 0.3, 0.4),
      beta.inter = rep(-1, 3),
      b.corr = 0.7, 
      b.corr.by = 4, 
      wei.shape = 1, 
      recr = 1, 
      fu = 10, 
      timefactor = 1, 
      active.main = Xprog, 
      active.inter = Xpred 
    )
    # adjust the simulated data frame columns to start time, status, treat, bm01... and turns the trt assignment +-0.5 to 0 and 1
    dataBM <- dataBM %>%
      dplyr::relocate(time, status, treat) %>%
      dplyr::mutate(treat = ifelse(treat == "-0.5", 1, ifelse(treat == "0.5", 0, treat))) 
    
    # Compute individual treatment effect
    z <- individual_trt_effect(
      data = dataBM,
      treat = "treat",
      time = "time",
      status = "status"
    )
    
    # Combine data and treatment effect
    combineddata <- base::cbind(z, dataBM) %>% 
      dplyr::select(-c(time, status, treat))
    
    # Set seed for random forest reproducibility
    set.seed(123)
    # Fit random forest and compute variable importance: meaning try to understand how biomarkers important to predict differenctial trt effect z
    important_features <- randomForest(z ~ .,
                                       data = combineddata,
                                       ntree = 500,
                                       importance = TRUE)
    
    # Extract variable importance
    importance_features_ordered <- important_features %>%
      randomForest::importance(type = 1, scale = TRUE) %>%
      as.data.frame() %>%
      dplyr::arrange(desc(`%IncMSE`))
    
    # Store results
    VT_simulation_results[[i]] <- importance_features_ordered
  }
  
  return(VT_simulation_results)
}


################################################  ANALYSIS    ##########################################################

results_VT_A3 <- run_VTsimulation_A3(num_simulations)

# Evaluation measures I created accordingly described in the paper Sechidis et al. 
# Extract variable names from each simulation trial. (this will give you Xpred^q)
top_q_variable_names_A3 <- lapply(results_VT_A3, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_VT_A[[1]])

# set of irrelevant features: obtained by subtracting prognostic and predictive biomarkers from all variable names
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))

# True Positive Rate: Captures the how accurately alygorithm captures the algorithm correctly identifies the predictive biomarkers. 
VT_TPR_A3 <- lapply(top_q_variable_names_A3, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})



# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
VT_FNR_prog_A3 <- lapply(top_q_variable_names_A3, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})


VT_FNR_irrlevant_A3 <- lapply(top_q_variable_names_A3, function(fnr_irr_calculation) {
  intersection <- intersect(Xirr, fnr_irr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})


# This actually the same as 1-TPR
VT_FNR_A3 <- base::Map("+", VT_FNR_prog_A3, VT_FNR_irrlevant_A3) 


# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
VT_TPR_A3 <- base::mean(base::unlist(VT_TPR_A3))                 
VT_FNR_prog_A3 <- base::mean(base::unlist(VT_FNR_prog_A3))



# Create a list to store variable importance scores for each simulation
importance_scores_list_A3 <- lapply(all_variable_names, function(var_name) {
  sapply(results_VT_A3, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "%IncMSE"])
    } else {
      return(NA)
    }
  })
})


#### Preparation for Visualizations ####
# Combine the importance scores into a data frame
importance_scores_df_A3 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_A3)), all_variable_names)

# This helps me to format data create plots. 
# PLOT1: Variable Importance Box plot for each variable obtained from N run Simulation 

importance_score_melted_VT_A3 <- reshape2::melt(importance_scores_df_A3) %>%
  dplyr::mutate(biomarker_types = dplyr::case_when(variable %in% Xpred ~ "Xpred",
                                                   variable %in% Xprog ~ "Xprog",
                                                   variable %in% Xirr ~ "Xirr"))


importance_score_melted_VT_A3$variable <- base::sub("^bm", "", importance_score_melted_VT_A3$variable)
VT_importanceplot_A3 <- ggplot2::ggplot(importance_score_melted_VT_A3, 
                                        aes(x = reorder(variable, -value), 
                                            y = value, 
                                            fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xprogpred" = "darkolivegreen3",
                                        "Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = sprintf("Scenario A3", 
                                num_simulations),
                x = "Biomarker-ID",
                y = "Importance",
                fill = "Biomarker Type") + 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_VT_A3$value), by = 10)) + 
  ggplot2::theme_bw() 
VT_importanceplot_A3



######### Setting-A: Scenario-4: 20 BM: 6 prognostic 14 Irrelevant ###########

Xprog <- c("bm01","bm02","bm03","bm04","bm05","bm06")     # Set of predictive biomarkers
                                                          # Set of prognostic biomarkers
# q <- length(Xpred)                                        # number of the predictive biomarkers
# XsolelyProg <- dplyr::setdiff(Xprog, Xpred)               # set differentiation 



# Function to run the simulation
run_VTsimulation_A4 <- function(num_simulations) {
  # Storage for results
  VT_simulation_results <- list()
  
  set.seed(123)
  # Loop through simulations
  for (i in 1:num_simulations) {
    
    # Generate data
    dataBM <- biospear::simdata(
      n = n,
      p = p, 
      q.main = 6 ,                          # the number of true prognostic biomarkers.
      q.inter = 0,                          # the number of true biomarkers interacting with the treatement.
      prob.tt = 0.5,                        # the treatment assignment probability.
      m0 = 5,                               # baseline median survival time. 
      alpha.tt = 0,                         # the effect of the treatment (in log-scale).
      beta.main = c(0.2, 0.3, 0.4, 0.2, 0.8, 0.3),
      b.corr = 0, 
      b.corr.by = 1, 
      wei.shape = 1, 
      recr = 1, 
      fu = 10, 
      timefactor = 1,                       # shape of the weibull dist,=1 exp. 
      active.main = Xprog 
    )
    # adjust the simulated data frame columns to start time, status, treat, bm01... and turns the trt assignment +-0.5 to 0 and 1
    dataBM <- dataBM %>%
      dplyr::relocate(time, status, treat) %>%
      dplyr::mutate(treat = ifelse(treat == "-0.5", 1, ifelse(treat == "0.5", 0, treat))) 
    
    # Compute individual treatment effect
    z <- individual_trt_effect(
      data = dataBM,
      treat = "treat",
      time = "time",
      status = "status"
    )
    
    # Combine data and treatment effect
    combineddata <- base::cbind(z, dataBM) %>% 
      dplyr::select(-c(time, status, treat))
    
    # Set seed for random forest reproducibility
    set.seed(123)
    # Fit random forest and compute variable importance: meaning try to understand how biomarkers important to predict differenctial trt effect z
    important_features <- randomForest(z ~ .,
                                       data = combineddata,
                                       ntree = 500,
                                       importance = TRUE)
    
    # Extract variable importance
    importance_features_ordered <- important_features %>%
      randomForest::importance(type = 1, scale = TRUE) %>%
      as.data.frame() %>%
      dplyr::arrange(desc(`%IncMSE`))
    
    # Store results
    VT_simulation_results[[i]] <- importance_features_ordered
  }
  
  return(VT_simulation_results)
}


################################################  ANALYSIS    ##########################################################

results_VT_A4 <- run_VTsimulation_A4(num_simulations)
#saveRDS(results_VT_A4, "results_VT_A4.rds")





# Create a list to store variable importance scores for each simulation
importance_scores_list_A4 <- lapply(all_variable_names, function(var_name) {
  sapply(results_VT_A4, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "%IncMSE"])
    } else {
      return(NA)
    }
  })
})


#### Preparation for Visualizations ####
# Combine the importance scores into a data frame
importance_scores_df_A4 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_A4)), all_variable_names)

# This helps me to format data create plots. 
# PLOT1: Variable Importance Box plot for each variable obtained from N run Simulation 


importance_score_melted_VT_A4 <- reshape2::melt(importance_scores_df_A4) %>%
  dplyr::mutate(biomarker_types = dplyr::case_when(variable %in% Xprog ~ "Xprog"))


importance_score_melted_VT_A4$variable <- base::sub("^bm", "", importance_score_melted_VT_A4$variable)
VT_importanceplot_A4 <- ggplot2::ggplot(importance_score_melted_VT_A4, 
                                        aes(x = reorder(variable, -value), 
                                            y = value, 
                                            fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = sprintf("Scenario A4", 
                                num_simulations),
                x = "Biomarker-ID",
                y = "Importance",
                fill = "Biomarker Type") + 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_VT_A3$value), by = 10)) + 
  ggplot2::theme_bw() 
VT_importanceplot_A4





######### Setting-A: Scenario-0: 20 BM: 20 Irrelevant ###########


# Function to run the simulation
run_VTsimulation_A0 <- function(num_simulations) {
  # Storage for results
  VT_simulation_results <- list()
  
  set.seed(123)
  # Loop through simulations
  for (i in 1:num_simulations) {
    
    # Generate data
    dataBM <- biospear::simdata(n = n,
                                p = p, 
                                q.main = 0 ,                          # the number of true prognostic biomarkers.
                                q.inter = 0,                          # the number of true biomarkers interacting with the treatement.
                                prob.tt = 0.5,                        # the treatment assignment probability.
                                m0 = 5,                               # baseline median survival time. 
                                alpha.tt = 0,                         # the effect of the treatment (in log-scale).
                                b.corr = 0, 
                                b.corr.by = 1, 
                                wei.shape = 1, 
                                recr = 1, 
                                fu = 10, 
                                timefactor = 1,                       # shape of the weibull dist,=1 exp. 
    )
    # adjust the simulated data frame columns to start time, status, treat, bm01... and turns the trt assignment +-0.5 to 0 and 1
    dataBM <- dataBM %>%
      dplyr::relocate(time, status, treat) %>%
      dplyr::mutate(treat = ifelse(treat == "-0.5", 1, ifelse(treat == "0.5", 0, treat))) 
    
    # Compute individual treatment effect
    z <- individual_trt_effect(
      data = dataBM,
      treat = "treat",
      time = "time",
      status = "status"
    )
    
    # Combine data and treatment effect
    combineddata <- base::cbind(z, dataBM) %>% 
      dplyr::select(-c(time, status, treat))
    
    # Set seed for random forest reproducibility
    set.seed(123)
    # Fit random forest and compute variable importance: meaning try to understand how biomarkers important to predict differenctial trt effect z
    important_features <- randomForest(z ~ .,
                                       data = combineddata,
                                       ntree = 500,
                                       importance = TRUE)
    
    # Extract variable importance
    importance_features_ordered <- important_features %>%
      randomForest::importance(type = 1, scale = TRUE) %>%
      as.data.frame() %>%
      dplyr::arrange(desc(`%IncMSE`))
    
    # Store results
    VT_simulation_results[[i]] <- importance_features_ordered
  }
  
  return(VT_simulation_results)
}


################################################  ANALYSIS    ##########################################################

results_VT_A0 <- run_VTsimulation_A0(num_simulations)
#saveRDS(results_VT_A0, "results_VT_A0.rds")


# Create a list to store variable importance scores for each simulation
importance_scores_list_A0 <- lapply(all_variable_names, function(var_name) {
  sapply(results_VT_A0, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "%IncMSE"])
    } else {
      return(NA)
    }
  })
})


#### Preparation for Visualizations ####
# Combine the importance scores into a data frame
importance_scores_df_A0 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_A0)), all_variable_names)

# This helps me to format data create plots. 
# PLOT1: Variable Importance Box plot for each variable obtained from N run Simulation 
# set of irrelevant biomarkers
Xirr <- all_variable_names

importance_score_melted_VT_A0 <- reshape2::melt(importance_scores_df_A0) %>%
  dplyr::mutate(biomarker_types = dplyr::case_when(variable %in% Xirr ~ "Xirr"))


importance_score_melted_VT_A0$variable <- base::sub("^bm", "", importance_score_melted_VT_A0$variable)
VT_importanceplot_A0 <- ggplot2::ggplot(importance_score_melted_VT_A0, 
                                        aes(x = reorder(variable, -value), 
                                            y = value, 
                                            fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xirr" = "grey")) +
  ggplot2::labs(title = sprintf("Scenario A0", 
                                num_simulations),
                x = "Biomarker-ID",
                y = "Importance",
                fill = "Biomarker Type") + 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_VT_A3$value), by = 10)) + 
  ggplot2::theme_bw() 
VT_importanceplot_A0








grid.arrange(VT_importanceplot_A0,
             VT_importanceplot_A1,
             VT_importanceplot_A2,
             VT_importanceplot_A3,
             VT_importanceplot_A4,
             ncol = 1)
             
             
             
             
             
             
             