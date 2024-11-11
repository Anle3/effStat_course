library(Rcpp)
library(xml2)
library(biospear)
library(dplyr)
library(ggplot2)
library(gridExtra)
dyn.load("sides64.dll")
source("CSIDES.r") 

num_simulations = 200 
n = 200
p = 20

Xprog <- c("bm01","bm02","bm03","bm04")
Xpred <- c("bm04","bm05","bm06")
Xprogpred <- dplyr::intersect(Xpred,Xprog)
XsolelyProg <- dplyr::setdiff(Xprog, Xpred)    # set differentiation 
q = length(Xpred)


run_SIDESsimulation_A1 <- function(num_simulations) {

  # Define empty list to store the simulation results 
  SIDES_simulation_results <- list()   

  set.seed(123)
  for (i in 1:num_simulations){
    
dataBM <- biospear::simdata(n = n,
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
                            timefactor = 1,                       # shape of the weibull dist,=1 exp. 
                            active.main = Xprog,
                            active.inter = Xpred
                            )

dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                        "control", 
                                                        "active"))                # preparing the data in suited form. -0.5 was 
  

biomarker_names <- colnames(dataBM)[2:21]
biomarker_names

# MIKA to SANDRA: Maybe below is the correct form I should use not the one above (which I used in my simulations already):
# Problem with above set of parameters when I parse I get the following error: <simpleError in get_standardizedSIDESparam(parameters): Error: missing censoring variable for survival outcome>
#<simpleError in create_xml(datafname, parameters): exiting..> 
# data_set_parameters = list(
#                            data_set = dataBM,
#                            outcome_variable_type = "survival",
#                            outcome_variable_name = "time", 
#                            outcome_cencor_name = "status",
#                            outcome_variable_direction = -1,
#                            treatment_variable_name = "treat",
#                            treatment_variable_control_value = "control",
#                            covariate_names = biomarker_names,
#                            covariate_types = rep("numeric", length(biomarker_names))
# )


  data_set_parameters = list(
  data_set = dataBM,
  outcome_variable_type = "binary",
  outcome_variable_name = "status",
  outcome_variable_direction = -1,
  treatment_variable_name = "treat",
  treatment_variable_control_value = "control",
  covariate_names = biomarker_names,
  covariate_types = rep("numeric", length(biomarker_names)
                         )
  )

  
  algorithm_parameters = list(min_subgroup_size = 30,
                             criterion_type = 1,                   # 1 means Differential effect splitting criterion.
                             depth = 3,
                             width = 5,
                             gamma = c(NA,NA,NA),                  # turning off the complexity control  
                             local_mult_adj = 1,
                             n_perms_mult_adjust = 250,
                             subgroup_search_algorithm = "Adaptive SIDEScreen procedure",
                             multiplier = 1,                       # the VI trhesholding rule is > E0+multipler*S0, E0 and S0 are means and SD from null distribution of max VI
                             n_perms_vi_score = 100                # number of permutations to compute max VI threshold
                             )


  parameters = list(data_set_parameters = data_set_parameters, algorithm_parameters = algorithm_parameters)

  CSIDES_results <- SIDES(parameters)
  VI_CIDES <- CSIDES_results[[2]] %>% dplyr::select(vi)
  
  
                             
  # Convert row names to numeric
  ordinal_biomarker <- as.numeric(rownames(VI_CIDES %>% as.data.frame()))
  
  # Format numbers with leading zeros and concatenate with "bm"
  nominal_biomarker <- paste0("bm", 
                          sprintf("%02d", ordinal_biomarker))
  
  # Assign new row names to the data frame
  rownames(VI_CIDES) <- nominal_biomarker
  
  SIDES_simulation_results[[i]] <- VI_CIDES
  }
return(SIDES_simulation_results)
}


results_CIDES_A1 <- run_SIDESsimulation_A1(num_simulations)


top_q_variable_names <- lapply(results_CIDES_A1, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_CIDES_A1[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_A1 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_A1, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_A1 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_A1)), 
                                    all_variable_names)


# set of irrelevant biomarkers
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_A1 <- reshape2::melt(importance_scores_df_CIDES_A1) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xprogpred ~ "Xprogpred",
      variable %in% Xpred ~ "Xpred",
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  )

importance_score_melted_CIDES_A1$variable <- base::sub("^bm", "", importance_score_melted_CIDES_A1$variable)
CIDES_importanceplot_A1 <- ggplot2::ggplot(importance_score_melted_CIDES_A1, 
                                        aes(x = reorder(variable, -value), 
                                            y = value, 
                                            fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xprogpred" = "darkolivegreen3",
                                        "Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario A1",
                x = "Biomarker",
                y = "Importance", 
                fill = "Biomarker Type")+ 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_A1$value), by = 1)) + 
  ggplot2::theme_bw()
CIDES_importanceplot_A1


# TPR 
CIDES_TPR_A1 <- lapply(top_q_variable_names, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})
CIDES_TPR_A1

# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
CIDES_meanTPR_A1 <- base::mean(base::unlist(CIDES_TPR_A1)) 


# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
CIDES_FNR_prog_A1 <- lapply(top_q_variable_names, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})

CIDES_meanFNR_prog_A1 <- base::mean(base::unlist(CIDES_FNR_prog_A1))






########################## CIDES SCENARIO A2 ###########################################

Xprog <- c("bm01","bm02","bm03")               # Set of predictive biomarkers
Xpred <- c("bm04","bm05","bm06")               # Set of prognostic biomarkers
Xprogpred <- dplyr::intersect(Xpred,Xprog)
XsolelyProg <- dplyr::setdiff(Xprog, Xpred)    # set differentiation 
q = length(Xpred)


run_SIDESsimulation_A2 <- function(num_simulations) {
  
  # Define empty list to store the simulation results 
  SIDES_simulation_results <- list()   
  
  set.seed(123)
  for (i in 1:num_simulations){
    
    dataBM <- biospear::simdata(n = n,
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
    
    dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                      "control", 
                                                      "active"))                # preparing the data in suited form. -0.5 was 
    
    
    biomarker_names <- colnames(dataBM)[2:21]
    biomarker_names
    
    
    data_set_parameters = list(
      data_set = dataBM,
      outcome_variable_type = "binary",
      outcome_variable_name = "status",
      outcome_variable_direction = -1,
      treatment_variable_name = "treat",
      treatment_variable_control_value = "control",
      covariate_names = biomarker_names,
      covariate_types = rep("numeric", length(biomarker_names)
      )
    )
    
    
    algorithm_parameters = list(min_subgroup_size = 30,
                                criterion_type = 1,                   # 1 means Differential effect splitting criterion.
                                depth = 3,
                                width = 5,
                                gamma = c(NA,NA,NA),                  # turning off the complexity control  
                                local_mult_adj = 1,
                                n_perms_mult_adjust = 250,
                                subgroup_search_algorithm = "Adaptive SIDEScreen procedure",
                                multiplier = 1,                       # the VI trhesholding rule is > E0+multipler*S0, E0 and S0 are means and SD from null distribution of max VI
                                n_perms_vi_score = 100                # number of permutations to compute max VI threshold
    )
    
    
    parameters = list(data_set_parameters = data_set_parameters, algorithm_parameters = algorithm_parameters)
    
    CSIDES_results <- SIDES(parameters)
    VI_CIDES <- CSIDES_results[[2]] %>% dplyr::select(vi)
    
    
    
    # Convert row names to numeric
    ordinal_biomarker <- as.numeric(rownames(VI_CIDES %>% as.data.frame()))
    
    # Format numbers with leading zeros and concatenate with "bm"
    nominal_biomarker <- paste0("bm", 
                                sprintf("%02d", ordinal_biomarker))
    
    # Assign new row names to the data frame
    rownames(VI_CIDES) <- nominal_biomarker
    
    SIDES_simulation_results[[i]] <- VI_CIDES
  }
  return(SIDES_simulation_results)
}


results_CIDES_A2 <- run_SIDESsimulation_A2(num_simulations)


top_q_variable_names_A2 <- lapply(results_CIDES_A2, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_CIDES_A2[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_A2 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_A2, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_A2 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_A2)), 
                                          all_variable_names)


# set of irrelevant biomarkers
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_A2 <- reshape2::melt(importance_scores_df_CIDES_A2) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xpred ~ "Xpred",
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  )

importance_score_melted_CIDES_A2$variable <- base::sub("^bm", "", importance_score_melted_CIDES_A2$variable)
CIDES_importanceplot_A2 <- ggplot2::ggplot(importance_score_melted_CIDES_A2, 
                                           aes(x = reorder(variable, -value), 
                                               y = value, 
                                               fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xprogpred" = "darkolivegreen3",
                                        "Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario A2",
                x = "Biomarker",
                y = "Importance", 
                fill = "Biomarker Type") + 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_A2$value), by = 1)) + 
  ggplot2::theme_bw() 
CIDES_importanceplot_A2


# TPR 
CIDES_TPR_A2 <- lapply(top_q_variable_names_A2, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})
CIDES_TPR_A2

# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
CIDES_meanTPR_A2 <- base::mean(base::unlist(CIDES_TPR_A2)) 


# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
CIDES_FNR_prog_A2 <- lapply(top_q_variable_names_A2, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})

CIDES_meanFNR_prog_A2 <- base::mean(base::unlist(CIDES_FNR_prog_A2))






########################## CIDES SCENARIO A3: A2 + correlation ###########################################

Xprog <- c("bm01","bm02","bm03")               # Set of predictive biomarkers
Xpred <- c("bm04","bm05","bm06")               # Set of prognostic biomarkers


run_SIDESsimulation_A3 <- function(num_simulations) {
  
  # Define empty list to store the simulation results 
  SIDES_simulation_results <- list()   
  
  set.seed(123)
  for (i in 1:num_simulations){
    
    dataBM <- biospear::simdata(n = n,
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
    
    dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                      "control", 
                                                      "active"))                # preparing the data in suited form. -0.5 was 
    
    
    biomarker_names <- colnames(dataBM)[2:21]
    biomarker_names
    
    
    data_set_parameters = list(
      data_set = dataBM,
      outcome_variable_type = "binary",
      outcome_variable_name = "status",
      outcome_variable_direction = -1,
      treatment_variable_name = "treat",
      treatment_variable_control_value = "control",
      covariate_names = biomarker_names,
      covariate_types = rep("numeric", length(biomarker_names)
      )
    )
    
    
    algorithm_parameters = list(min_subgroup_size = 30,
                                criterion_type = 1,                   # 1 means Differential effect splitting criterion.
                                depth = 3,
                                width = 5,
                                gamma = c(NA,NA,NA),                  # turning off the complexity control  
                                local_mult_adj = 1,
                                n_perms_mult_adjust = 250,
                                subgroup_search_algorithm = "Adaptive SIDEScreen procedure",
                                multiplier = 1,                       # the VI trhesholding rule is > E0+multipler*S0, E0 and S0 are means and SD from null distribution of max VI
                                n_perms_vi_score = 100                # number of permutations to compute max VI threshold
    )
    
    
    parameters = list(data_set_parameters = data_set_parameters, algorithm_parameters = algorithm_parameters)
    
    CSIDES_results <- SIDES(parameters)
    VI_CIDES <- CSIDES_results[[2]] %>% dplyr::select(vi)
    
    
    
    # Convert row names to numeric
    ordinal_biomarker <- as.numeric(rownames(VI_CIDES %>% as.data.frame()))
    
    # Format numbers with leading zeros and concatenate with "bm"
    nominal_biomarker <- paste0("bm", 
                                sprintf("%02d", ordinal_biomarker))
    
    # Assign new row names to the data frame
    rownames(VI_CIDES) <- nominal_biomarker
    
    SIDES_simulation_results[[i]] <- VI_CIDES
  }
  return(SIDES_simulation_results)
}


results_CIDES_A3 <- run_SIDESsimulation_A3(num_simulations)


top_q_variable_names_A3 <- lapply(results_CIDES_A3, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_CIDES_A3[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_A3 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_A3, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_A3 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_A3)), 
                                          all_variable_names)


# set of irrelevant biomarkers
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_A3 <- reshape2::melt(importance_scores_df_CIDES_A3) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xpred ~ "Xpred",
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  )

importance_score_melted_CIDES_A3$variable <- base::sub("^bm", "", importance_score_melted_CIDES_A3$variable)
CIDES_importanceplot_A3 <- ggplot2::ggplot(importance_score_melted_CIDES_A3, 
                                           aes(x = reorder(variable, -value), 
                                               y = value, 
                                               fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xprogpred" = "darkolivegreen3",
                                        "Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario A3",
                x = "Biomarker",
                y = "Importance", 
                fill = "Biomarker Type") + 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_A3$value), by = 1)) + 
  ggplot2::theme_bw() 
CIDES_importanceplot_A3


# TPR 
CIDES_TPR_A3 <- lapply(top_q_variable_names_A3, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})
CIDES_TPR_A3

# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
CIDES_meanTPR_A3 <- base::mean(base::unlist(CIDES_TPR_A3)) 


# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
CIDES_FNR_prog_A3 <- lapply(top_q_variable_names_A3, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})

CIDES_meanFNR_prog_A3 <- base::mean(base::unlist(CIDES_FNR_prog_A3))





########################## CIDES SCENARIO A4: ###########################################

Xprog <- c("bm01","bm02","bm03","bm04","bm05","bm06")     # Set of predictive biomarkers


run_SIDESsimulation_A4 <- function(num_simulations) {
  
  # Define empty list to store the simulation results 
  SIDES_simulation_results <- list()   
  
  set.seed(123)
  for (i in 1:num_simulations){
    
    dataBM <- biospear::simdata(n = n,
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
    
    dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                      "control", 
                                                      "active"))                # preparing the data in suited form. -0.5 was 
    
    
    biomarker_names <- colnames(dataBM)[2:21]
    biomarker_names
    
    
    data_set_parameters = list(
      data_set = dataBM,
      outcome_variable_type = "binary",
      outcome_variable_name = "status",
      outcome_variable_direction = -1,
      treatment_variable_name = "treat",
      treatment_variable_control_value = "control",
      covariate_names = biomarker_names,
      covariate_types = rep("numeric", length(biomarker_names)
      )
    )
    
    
    algorithm_parameters = list(min_subgroup_size = 30,
                                criterion_type = 1,                   # 1 means Differential effect splitting criterion.
                                depth = 3,
                                width = 5,
                                gamma = c(NA,NA,NA),                  # turning off the complexity control  
                                local_mult_adj = 1,
                                n_perms_mult_adjust = 250,
                                subgroup_search_algorithm = "Adaptive SIDEScreen procedure",
                                multiplier = 1,                       # the VI trhesholding rule is > E0+multipler*S0, E0 and S0 are means and SD from null distribution of max VI
                                n_top_biomarkers = 6,
                                n_perms_vi_score = 100                # number of permutations to compute max VI threshold
    )
    
    
    parameters = list(data_set_parameters = data_set_parameters, algorithm_parameters = algorithm_parameters)
    
    CSIDES_results <- SIDES(parameters)
    VI_CIDES <- CSIDES_results[[2]] %>% dplyr::select(vi)
    
    
    
    # Convert row names to numeric
    ordinal_biomarker <- as.numeric(rownames(VI_CIDES %>% as.data.frame()))
    
    # Format numbers with leading zeros and concatenate with "bm"
    nominal_biomarker <- paste0("bm", 
                                sprintf("%02d", ordinal_biomarker))
    
    # Assign new row names to the data frame
    rownames(VI_CIDES) <- nominal_biomarker
    
    SIDES_simulation_results[[i]] <- VI_CIDES
  }
  return(SIDES_simulation_results)
}


results_CIDES_A4 <- run_SIDESsimulation_A4(num_simulations)


#top_q_variable_names_A4 <- lapply(results_CIDES_A4, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_CIDES_A4[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_A4 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_A4, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_A4 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_A4)), 
                                          all_variable_names)


# set of irrelevant biomarkers
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_A4 <- reshape2::melt(importance_scores_df_CIDES_A4) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  )

importance_score_melted_CIDES_A4$variable <- base::sub("^bm", "", importance_score_melted_CIDES_A4$variable)
CIDES_importanceplot_A4 <- ggplot2::ggplot(importance_score_melted_CIDES_A4, 
                                           aes(x = reorder(variable, -value), 
                                               y = value, 
                                               fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario A4" ,
                x = "Biomarker",
                y = "Importance", 
                fill = "Biomarker Type") + 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_A4$value), by = 1)) + 
  ggplot2::theme_bw() 
CIDES_importanceplot_A4





########################## CIDES SCENARIO A0: NULL ###########################################



run_SIDESsimulation_A0 <- function(num_simulations) {
  
  # Define empty list to store the simulation results 
  SIDES_simulation_results <- list()   
  
  set.seed(123)
  for (i in 1:num_simulations){
    
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
    
    dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                      "control", 
                                                      "active"))                # preparing the data in suited form. -0.5 was 
    
    
    biomarker_names <- colnames(dataBM)[2:21]
    biomarker_names
    
    
    data_set_parameters = list(
      data_set = dataBM,
      outcome_variable_type = "binary",
      outcome_variable_name = "status",
      outcome_variable_direction = -1,
      treatment_variable_name = "treat",
      treatment_variable_control_value = "control",
      covariate_names = biomarker_names,
      covariate_types = rep("numeric", length(biomarker_names)
      )
    )
    
    
    algorithm_parameters = list(min_subgroup_size = 30,
                                criterion_type = 1,                   # 1 means Differential effect splitting criterion.
                                depth = 3,
                                width = 5,
                                gamma = c(NA,NA,NA),                  # turning off the complexity control  
                                local_mult_adj = 1,
                                n_perms_mult_adjust = 250,
                                subgroup_search_algorithm = "Adaptive SIDEScreen procedure",
                                multiplier = 1,                       # the VI trhesholding rule is > E0+multipler*S0, E0 and S0 are means and SD from null distribution of max VI
                                n_perms_vi_score = 100                # number of permutations to compute max VI threshold
    )
    
    
    parameters = list(data_set_parameters = data_set_parameters, algorithm_parameters = algorithm_parameters)
    
    CSIDES_results <- SIDES(parameters)
    VI_CIDES <- CSIDES_results[[2]] %>% dplyr::select(vi)
    
    
    
    # Convert row names to numeric
    ordinal_biomarker <- as.numeric(rownames(VI_CIDES %>% as.data.frame()))
    
    # Format numbers with leading zeros and concatenate with "bm"
    nominal_biomarker <- paste0("bm", 
                                sprintf("%02d", ordinal_biomarker))
    
    # Assign new row names to the data frame
    rownames(VI_CIDES) <- nominal_biomarker
    
    SIDES_simulation_results[[i]] <- VI_CIDES
  }
  return(SIDES_simulation_results)
}


results_CIDES_A0 <- run_SIDESsimulation_A0(num_simulations)


all_variable_names <- rownames(results_CIDES_A0[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_A0 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_A0, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_A0 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_A0)), 
                                          all_variable_names)


# set of irrelevant biomarkers
Xirr <- all_variable_names



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_A0 <- reshape2::melt(importance_scores_df_CIDES_A0) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xirr ~ "Xirr"
    )
  )


importance_score_melted_CIDES_A0$variable <- base::sub("^bm", "", importance_score_melted_CIDES_A0$variable)
CIDES_importanceplot_A0 <- ggplot2::ggplot(importance_score_melted_CIDES_A0, 
                                           aes(x = reorder(variable, -value), 
                                               y = value, 
                                               fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario A0",
                x = "Biomarker",
                y = "Importance", 
                fill = "Biomarker Type") + 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_A0$value), by = 1)) + 
  ggplot2::theme_bw() 
CIDES_importanceplot_A0


grid.arrange(CIDES_importanceplot_A0,
             CIDES_importanceplot_A1,
             CIDES_importanceplot_A2,
             CIDES_importanceplot_A3,
             CIDES_importanceplot_A4,
             ncol = 1)



#saveRDS(results_CIDES_A4, "results_CIDES_A4.rds")
#saveRDS(results_CIDES_A3, "results_CIDES_A3.rds")
# saveRDS(results_CIDES_A2, "results_CIDES_A2.rds")
# saveRDS(results_CIDES_A1, "results_CIDES_A1.rds")
# saveRDS(results_CIDES_A0, "results_CIDES_A0.rds")
# 
results_CIDES_A0 <- readRDS("results_CIDES_A0.rds")
results_CIDES_A1 <- readRDS("results_CIDES_A1.rds")
results_CIDES_A2 <- readRDS("results_CIDES_A2.rds")
results_CIDES_A3 <- readRDS("results_CIDES_A3.rds")
results_CIDES_A4 <- readRDS("results_CIDES_A4.rds")


