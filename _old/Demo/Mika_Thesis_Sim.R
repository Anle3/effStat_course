# Load the simstudy package
library(simstudy)
library(ggplot2)
library(survminer)
library(survival)
library(tidyr)
library(glmnet)
library(tidymodels)

# Set the seed for reproducibility
set.seed(123)

# Define the number of observations
n <- 20

# Define binary Treatment variable (equally sized)
def <- defData(varname = "trt", formula = "0.5", dist = "binary")

#Define biomarker variables as characteristic
for (i in 1:30) {
  varname <- paste("X", i, sep = "")
  def <- defData(def, varname = varname, formula = "0", variance = "1", dist = "normal", link = "log")
}


# Interaction-Term for predictive biomarkers (trt*biomarker)
def <- defData(def, varname = "X28TimesTreatment", formula = "X28 * trt") 
def <- defData(def, varname = "X29TimesTreatment", formula = "X29 * trt") 
def <- defData(def, varname = "X30TimesTreatment", formula = "X30 * trt")

# FIGUREOUT THIS PART WHEN YOU INCREASE PREDICTIVE MARKER!!! Worse case loop it.
# for (i in c(28, 29, 30)) {
#   varname <- paste("X", i, "TimesTreatment", sep = "")
#   formula <- paste("X", i, " * trt", sep = "")
#   def <- defData(def, varname = varname, formula = formula)
# }


# Define the formula for generating the survival time
survival_formula <- "exp(0.5*trt + 0.2*X1 + 0.3*X2 + 0.4*X3 + 0.2*X4 + 0.3*X5 + 0.4*X6 + 0.5*X7 + 0.4*X8 + 0.2*X9 + 0.3*X10 
                             + 0.2*X11 + 0.3*X12 + 0.4*X13 + 0.2*X14 + 0.3*X15 + 0.4*X16 + 0.5*X17 + 0.4*X18 + 0.2*X19 + 0.3*X20 
                             + 0.2*X21 + 0.3*X22 + 0.4*X23 + 0.2*X24 + 0.3*X25 + 0.4*X26 + 0.5*X27 + 0.4*X28 + 0.2*X29 + 0.3*X30 
                             - 1*X28*trt- 1*X29*trt- 1*X30*trt)"

# Survival information Time and Cencoring Time
# scale: Scales lambda paramtr for the Weibull dist: smaller leads short survtimes vice versa
# shape: Shapes of the survtime distribution: Shape = 1 for an exp dist. Shape < 1 the hazard decreases/event more likely to occur. 
sdef <- defSurv(varname = "surv_time", formula = survival_formula, scale = 100, shape = 1)
sdef <- defSurv(sdef, varname = "censorTime", scale = 100, shape = 1)


# Generate survival time with censoring (Weibull Distribution used to generate survival times)
set.seed(1234)
base_data <- genData(n, def)
base_data <- genSurv(base_data, sdef, timeName = "obsTime", censorName = "censorTime",
                     eventName = "status", keepEvents = TRUE)
base_data


############################# Visualization #############################


# Histogram of the variables
par(mfrow=c(5, 6))  # To arrange plots in a grid
for (i in 1:30) {
  hist(base_data[[paste("X", i, sep = "")]], main = paste("X", i), xlab = "")
}


# Box plot of the variables
base_data_long_box <- tidyr::gather(base_data, biomarker, value, starts_with("X"))

ggplot(base_data_long_box, aes(x = factor(trt), y = value, fill = factor(trt))) +
  geom_boxplot() +
  facet_wrap(~reorder(biomarker, as.numeric(gsub("X", "", biomarker))), scales = "fixed") +
  labs(
    x = "Treatment Group ",
    y = "Biomarker Value",
    fill = "Group",
    title = "Box Plots of Biomarkers by Treatment Group"
    ) + theme_bw()

# In case in ERROR: try to adjust your Plots part of R studio bigger scale it by holding and dragging with the cursor.
# In case in ERROR: try to adjust your Plots part of R studio bigger scale it by holding and dragging with the cursor.
# In case in ERROR: try to adjust your Plots part of R studio bigger scale it by holding and dragging with the cursor.


# Kaplan-Maier Survival cure for trt groups
fit1 <- survival::survfit(Surv(base_data$obsTime, base_data$status) ~ base_data$trt, data = base_data)
survminer::ggsurvplot(fit1, data = base_data,
                      size = .5)


# Hazard? 
fit2 <- survival::coxph(Surv(base_data$obsTime, base_data$status) ~ base_data$X2, data = base_data)
gtsummary::tbl_regression(fit2) # optional to read statistics


# Create a ggplot scatter plot
base_data_long <- tidyr::pivot_longer(base_data, cols = c(X26,X27,X28, X29, X30), names_to = "Biomarker", values_to = "Value")

ggplot(base_data_long, aes(x = Value, y = surv_time, color = factor(trt))) +
  geom_point(shape = 19) +
  labs(
    x = "Biomarker Value",
    y = "Survival Time",
    title = "Scatter Plot of Biomarkers vs. Survival Time"
  ) +
  facet_wrap(~ Biomarker, scales = "free_x") +
  scale_color_discrete(name = "Group") +
  theme_minimal()


# Extract the subset of columns from X1 to X30
subset_data <- base_data[, c(paste0("X", 1:30), "X28TimesTreatment", "X29TimesTreatment", "X30TimesTreatment")]

# Calculate the correlation matrix
correlation_matrix <- cor(subset_data)

# Create a larger correlation heatmap using the heatmap function
heatmap(correlation_matrix,
        col = colorRampPalette(c("blue", "white", "red"))(100),
        main = "Correlation Heatmap",
        cexRow = 1,  # Increase the row label size
        cexCol = 1,   # Increase the column label size
        lhei = c(0.5, 2.5),  # Control height of row labels (0.5 for diagonal)
        lwid = c(0.5, 2.5)  # Control width of column labels (0.5 for diagonal)
)


####### LASSO ####### 

library(tidymodels)
set.seed(123)


# Splitting the data into test (%30) and training (%70)
splitted_data <- initial_split(base_data, prop = 0.7)
training_data <- rsample::training(splitted_data)
testing_data <- rsample::testing(splitted_data)


# Defining the outcome and predictors
outcome <- "surv_time"
predictors <- c(paste0("X", 1:30))#, "X28TimesTreatment", "X29TimesTreatment", "X30TimesTreatment") # interactionterm excluded from set of predictors


# Defining the Recipe: preprocessing + formula for the model 
recipe_data <- recipes::recipe(formula = as.formula(paste(outcome, "~", paste(predictors, collapse = " + "))),
                               data = training_data) %>%
  step_scale(all_predictors()) %>% # scales and centeres the predictors might not be needed in simulation setting since covariates created as desired.
  step_center(all_predictors())


# Speficiy your model: LASSO
specification_lasso <- parsnip::linear_reg(penalty = 1, mixture = 1) %>% 
  parsnip::set_engine("glm")


# Create the Workflow : Combine the recipe and the model(LASSO)
workflow_lasso <- workflow() %>%
  add_recipe(recipe_data) %>%
  add_model(specification_lasso)


# Fit the model (if neccessary tune it later)
fitted_lasso <- fit(workflow_lasso, data = training_data)


# Prediction
predictions_lasso <- predict(fitted_lasso, new_data = testing_data) %>%
  as_tibble() %>%
  dplyr::bind_cols(testing_data)

# Problem!! Check the Correlation btween predictors! Obvious trt*biomarker but I have not 
# included trt*biomarker inside of my predictors (line 145) They are in the training data.

cor_matrix_lasso_predictors <- cor(training_data[, c(paste0("X", 1:30), "X28TimesTreatment", 
                                                     "X29TimesTreatment", "X30TimesTreatment")])

highly_correlated_pairs <- which(cor_matrix_lasso_predictors > 0.7 & cor_matrix_lasso_predictors < 1, arr.ind = TRUE)


####### LASSO FAILED CORRELATED FEATURES #######

###### ELASTICNET #######

spesification_elasticnet <- linear_reg(penalty = 1, mixture = 0.5) %>%
  set_engine("glmnet")


workflow_elasticnet <- workflow() %>%
  add_recipe(recipe_data) %>%
  add_model(spesification_elasticnet)


fitted_elasticnet <- fit(workflow_elasticnet, data = training_data)

#I am not sure how to interpret What I want from this model. 
predictions_elasticnet <- predict(fitted_elasticnet, new_data = testing_data) %>%
  bind_cols(testing_data) %>%
  select(.pred, id, status)


###### I will fix continue with below models #####
####### PP-LASSO (Priority Penalized LASSO) ####### Prognostic Predictive Lasso for Biomarker Selection

#You assign the spesific weights for spesific variables. There is existed R package for this I will 
#in my modelling scheme as selected penalized models. https://pubmed.ncbi.nlm.nih.gov/36690931/

####### IPP-LASSO (Integrative LASSO with penalthy factors) #######
#https://cran.r-project.org/web/packages/ipflasso/ipflasso.pdf
