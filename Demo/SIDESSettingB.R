library(Rcpp)
library(xml2)
library(biospear)
library(dplyr)
library(ggplot2)
library(gridExtra)
dyn.load("sides64.dll")
source("CSIDES.r") 

# User defined simulation settings 
num_simulations <- 200
n = 70
p = 250



run_SIDESsimulation_B1 <- function(num_simulations) {

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
                            active.main = c("bm001","bm002","bm003","bm004"),
                            active.inter = c("bm004","bm005","bm006")
                            )


dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                  "control", 
                                                  "active"))                # preparing the data in suited form. -0.5 was 


biomarker_names <- colnames(dataBM)[2:251]
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


results_CIDES_B1 <- run_SIDESsimulation_B1(num_simulations)


Xprog <- c("bm01","bm02","bm03","bm04")
Xpred <- c("bm04","bm05","bm06")
Xprogpred <- dplyr::intersect(Xpred,Xprog)
XsolelyProg <- dplyr::setdiff(Xprog, Xpred)    # set differentiation 
q = length(Xpred)


top_q_variable_names_B1 <- lapply(results_CIDES_B1, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_CIDES_B1[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_B1 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_B1, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_B1 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_B1)), 
                                    all_variable_names)

importance_scores_df_CIDES_B1_reordered <- importance_scores_df_CIDES_B1[, 
                                                               order(-colSums(importance_scores_df_CIDES_B1))][,1:20]


# set of irrelevant biomarkers
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_B1 <- reshape2::melt(importance_scores_df_CIDES_B1_reordered) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xprogpred ~ "Xprogpred",
      variable %in% Xpred ~ "Xpred",
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  ) 


importance_score_melted_CIDES_B1$variable <- base::sub("^bm", "", importance_score_melted_CIDES_B1$variable)
CIDES_importanceplot_B1 <- ggplot2::ggplot(importance_score_melted_CIDES_B1, 
                                        aes(x = reorder(variable, -value), 
                                            y = value, 
                                            fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xprogpred" = "darkolivegreen3",
                                        "Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario B1",
               x = "Biomarker-ID",
                y = "Importance", 
                fill = "Biomarker Type")+ 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_B1$value), by = 0.5)) + 
  ggplot2::theme_bw()
CIDES_importanceplot_B1


# TPR 
CIDES_TPR_B1 <- lapply(top_q_variable_names_B1, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})
CIDES_TPR_B1

# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
CIDES_meanTPR_B1 <- base::mean(base::unlist(CIDES_TPR_B1)) 


# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
CIDES_FNR_prog_B1 <- lapply(top_q_variable_names_B1, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})

CIDES_meanFNR_prog_B1 <- base::mean(base::unlist(CIDES_FNR_prog_B1))






########################## CIDES Setting-B: Scenario-2: 250 BM: 3-Prog, 3-Pred, 244 Irrelevant ###########


run_SIDESsimulation_B2 <- function(num_simulations) {
  
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
                                timefactor = 1,                       # shape of the weibull dist,=1 exp. 
                                active.main = c("bm001","bm002","bm003"),
                                active.inter = c("bm004","bm005","bm006")
    )
    
    dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                      "control", 
                                                      "active"))                # preparing the data in suited form. -0.5 was 
    
    
    biomarker_names <- colnames(dataBM)[2:251]
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


results_CIDES_B2 <- run_SIDESsimulation_B2(num_simulations)


Xprog <- c("bm01","bm02","bm03")
Xpred <- c("bm04","bm05","bm06")
Xprogpred <- dplyr::intersect(Xpred,Xprog)
XsolelyProg <- dplyr::setdiff(Xprog, Xpred)    # set differentiation 
q = length(Xpred)


top_q_variable_names_B2 <- lapply(results_CIDES_B2, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_CIDES_B2[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_B2 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_B2, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_B2 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_B2)), 
                                          all_variable_names)

importance_scores_df_CIDES_B2_reordered <- importance_scores_df_CIDES_B2[, 
                                                                         order(-colSums(importance_scores_df_CIDES_B2))][,1:20]


# set of irrelevant biomarkers
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_B2 <- reshape2::melt(importance_scores_df_CIDES_B2_reordered) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xprogpred ~ "Xprogpred",
      variable %in% Xpred ~ "Xpred",
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  ) 


importance_score_melted_CIDES_B2$variable <- base::sub("^bm", "", importance_score_melted_CIDES_B2$variable)
CIDES_importanceplot_B2 <- ggplot2::ggplot(importance_score_melted_CIDES_B2, 
                                           aes(x = reorder(variable, -value), 
                                               y = value, 
                                               fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario B2",
               x = "Biomarker-ID",
                y = "Importance", 
                fill = "Biomarker Type")+ 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_B2$value), by = 0.5)) + 
  ggplot2::theme_bw()
CIDES_importanceplot_B2


# TPR 
CIDES_TPR_B2 <- lapply(top_q_variable_names_B2, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})
CIDES_TPR_B2

# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
CIDES_meanTPR_B2 <- base::mean(base::unlist(CIDES_TPR_B2)) 


# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
CIDES_FNR_prog_B2 <- lapply(top_q_variable_names_B2, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})

CIDES_meanFNR_prog_B2 <- base::mean(base::unlist(CIDES_FNR_prog_B2))






########################## CIDES SCENARIO B3: B2 + correlation ###########################################

run_SIDESsimulation_B3 <- function(num_simulations) {
  
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
                                b.corr.by = 25, 
                                wei.shape = 1, 
                                recr = 1, 
                                fu = 10, 
                                timefactor = 1,                       # shape of the weibull dist,=1 exp. 
                                active.main = c("bm001","bm002","bm003"),
                                active.inter = c("bm004","bm005","bm006"))
    
    
    dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                      "control", 
                                                      "active"))                # preparing the data in suited form. -0.5 was 

    
    biomarker_names <- colnames(dataBM)[2:251]
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


results_CIDES_B3 <- run_SIDESsimulation_B3(num_simulations)


top_q_variable_names_B3 <- lapply(results_CIDES_B3, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_CIDES_B3[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_B3 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_B3, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_B3 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_B3)), 
                                          all_variable_names)

importance_scores_df_CIDES_B3_reordered <- importance_scores_df_CIDES_B3[, 
                                                                         order(-colSums(importance_scores_df_CIDES_B3))][,1:20]


# set of irrelevant biomarkers
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_B3 <- reshape2::melt(importance_scores_df_CIDES_B3_reordered) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xpred ~ "Xpred",
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  ) 



importance_score_melted_CIDES_B3$variable <- base::sub("^bm", "", importance_score_melted_CIDES_B3$variable)
CIDES_importanceplot_B3 <- ggplot2::ggplot(importance_score_melted_CIDES_B3, 
                                           aes(x = reorder(variable, -value), 
                                               y = value, 
                                               fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario B3",
               x = "Biomarker-ID",
                y = "Importance", 
                fill = "Biomarker Type")+ 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_B3$value), by = 0.5)) + 
  ggplot2::theme_bw()
CIDES_importanceplot_B3


# TPR 
CIDES_TPR_B3 <- lapply(top_q_variable_names_B3, function(tpr_calculation) {
  intersection <- intersect(Xpred, tpr_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})
CIDES_TPR_B3

# Takes the TPR produced in each simulation and get the means to receive single TPR score for the scenario.
CIDES_meanTPR_B3 <- base::mean(base::unlist(CIDES_TPR_B3)) 


# FNR_prognostic: captures how often an algorithm selects as predictive covariates those that are solely prognostic.
CIDES_FNR_prog_B3 <- lapply(top_q_variable_names_B3, function(fnr_prog_calculation) {
  intersection <- intersect(XsolelyProg, fnr_prog_calculation)
  cardinality <- length(intersection)
  ratio <- cardinality / length(Xpred)
  ratio
})

CIDES_meanFNR_prog_B3 <- base::mean(base::unlist(CIDES_FNR_prog_B3))






######### Setting-B: Scenario-0: 250 BM: 250 Irrelevant ###########


run_SIDESsimulation_B0 <- function(num_simulations) {
  
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
    
    
    biomarker_names <- colnames(dataBM)[2:251]
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


results_CIDES_B0 <- run_SIDESsimulation_B0(num_simulations)

all_variable_names <- rownames(results_CIDES_B0[[1]])
Xirr <- all_variable_names



#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_B0 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_B0, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_B0 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_B0)), 
                                          all_variable_names)

importance_scores_df_CIDES_B0_reordered <- importance_scores_df_CIDES_B0[, 
                                                                         order(-colSums(importance_scores_df_CIDES_B0))][,1:20]

# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_B0 <- reshape2::melt(importance_scores_df_CIDES_B0_reordered) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xirr ~ "Xirr"
    )
  ) 


importance_score_melted_CIDES_B0$variable <- base::sub("^bm", "", importance_score_melted_CIDES_B0$variable)
CIDES_importanceplot_B0 <- ggplot2::ggplot(importance_score_melted_CIDES_B0, 
                                           aes(x = reorder(variable, -value), 
                                               y = value, 
                                               fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario B0",
               x = "Biomarker-ID",
                y = "Importance", 
                fill = "Biomarker Type")+ 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_B0$value), by = 0.5)) + 
  ggplot2::theme_bw()
CIDES_importanceplot_B0












######### Setting-B: Scenario-4: 250 BM: 6 Prognostic ###########

run_SIDESsimulation_B4 <- function(num_simulations) {
  
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
                                active.main = c("bm001","bm002","bm003","bm004","bm005","bm006"))
    
    
    
    dataBM <- dataBM %>% dplyr::mutate(treat = ifelse(treat == "-0.5",          # before: A total of n patients is generated and randomly assigned to the experimental (coded as +0.5, with probability prob.tt) and control treatment (coded as -0.5
                                                      "control", 
                                                      "active"))                # preparing the data in suited form. -0.5 was 
    
    
    biomarker_names <- colnames(dataBM)[2:251]
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


results_CIDES_B4 <- run_SIDESsimulation_B4(num_simulations)

Xprog <- c("bm01","bm02","bm03","bm04","bm05","bm06") 
top_q_variable_names_B4 <- lapply(results_CIDES_B4, function(df) rownames(df)[1:q])
all_variable_names <- rownames(results_CIDES_B4[[1]])


#### Preparation for Visualizations ####
# Collect vi score for each biomarker from each simulation run.
importance_scores_list_CIDES_B4 <- lapply(all_variable_names, function(var_name) {
  sapply(results_CIDES_B4, function(sim_result) {
    if (var_name %in% rownames(sim_result)) {
      return(sim_result[var_name, "vi"])
    } else {
      return(NA)
    }
  })
})


# Turn the previous step into a nice data frame so you can use. 
importance_scores_df_CIDES_B4 <- setNames(as.data.frame(do.call(cbind, importance_scores_list_CIDES_B4)), 
                                          all_variable_names)

importance_scores_df_CIDES_B4_reordered <- importance_scores_df_CIDES_B4[, 
                                                                         order(-colSums(importance_scores_df_CIDES_B4))][,1:20]


# set of irrelevant biomarkers
Xirr <- dplyr::setdiff(all_variable_names,                        
                       dplyr::union(Xpred, Xprog))



# label your biomarkers in their types so you can argue better. 
importance_score_melted_CIDES_B4 <- reshape2::melt(importance_scores_df_CIDES_B4_reordered) %>%
  dplyr::mutate(
    biomarker_types = dplyr::case_when(
      variable %in% Xpred ~ "Xpred",
      variable %in% Xprog ~ "Xprog",
      variable %in% Xirr ~ "Xirr"
    )
  ) 



importance_score_melted_CIDES_B4$variable <- base::sub("^bm", "", importance_score_melted_CIDES_B4$variable)
CIDES_importanceplot_B4 <- ggplot2::ggplot(importance_score_melted_CIDES_B4, 
                                           aes(x = reorder(variable, -value), 
                                               y = value, 
                                               fill = biomarker_types)) +
  ggplot2::geom_boxplot(alpha = 0.5) +
  ggplot2::scale_fill_manual(values = c("Xpred" = "slateblue", 
                                        "Xprog" = "#EE8844",  
                                        "Xirr" = "grey")) +
  ggplot2::labs(title = "Scenario B4",
                x = "Biomarker-ID",
                y = "Importance", 
                fill = "Biomarker Type")+ 
  ggplot2::scale_y_continuous(breaks = seq(0, max(importance_score_melted_CIDES_B4$value), by = 0.5)) + 
  ggplot2::theme_bw()
CIDES_importanceplot_B4



grid.arrange(CIDES_importanceplot_B0,
             CIDES_importanceplot_B1,
             CIDES_importanceplot_B2,
             CIDES_importanceplot_B3,
             CIDES_importanceplot_B4,
             ncol = 1)


# saveRDS(results_CIDES_B4, "results_CIDES_B4.rds")
# saveRDS(results_CIDES_B3, "results_CIDES_B3.rds")
# saveRDS(results_CIDES_B2, "results_CIDES_B2.rds")
# saveRDS(results_CIDES_B1, "results_CIDES_B1.rds")
# saveRDS(results_CIDES_B0, "results_CIDES_B0.rds")
# 
# results_CIDES_B0 <- readRDS("results_CIDES_B0.rds")
# results_CIDES_B1 <- readRDS("results_CIDES_B1.rds")
# results_CIDES_B2 <- readRDS("results_CIDES_B2.rds")
# results_CIDES_B3 <- readRDS("results_CIDES_B3.rds")
# results_CIDES_B4 <- readRDS("results_CIDES_B4.rds")




