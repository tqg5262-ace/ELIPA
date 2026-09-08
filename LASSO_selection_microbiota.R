# =============================================================================
# Title:   Longitudinal Feature Velocity Pipeline & LASSO Feature Selection
# Purpose: 1) Compute rate-of-change ("velocity") features across study
#             time points.
#          2) Regress out selected covariates and apply an inverse-normal
#             transformation to the resulting residuals.
#          3) Filter highly correlated features within each study period.
#          4) Apply LASSO regression with leave-one-out cross-validation
#             (LOOCV) to identify predictors associated with the continuous
#             outcome.
#
# NOTE ON ANONYMIZATION: All local file paths, working directories, source
# filenames, and data-type-specific labels have been replaced with generic
# placeholders. Update the "User-defined inputs" section below before running.
# =============================================================================


# ---- 1. Required packages -----------------------------------------------------

library(caret)
library(dplyr)
library(glmnet)
library(readxl)
library(tidyr)


# ---- 2. User-defined inputs ---------------------------------------------------

# Update these paths/parameters for your local environment before running.

longitudinal_feature_file <- "path/to/longitudinal_feature_data.xlsx"
participant_data_file     <- "path/to/participant_level_data.csv"

# Column names expected in the anonymised input data:
id_col        <- "ID"
time_col      <- "Time"
age_col       <- "age"
sex_col       <- "sex"
group_col     <- "group"
outcome_col   <- "outcome"

# Correlation threshold used to remove highly correlated features.
correlation_cutoff <- 0.80

# Seed used for LASSO cross-validation.
lasso_seed <- 123

longitudinal_data <- read_excel(longitudinal_feature_file)
participant_data  <- read.csv(participant_data_file)

covariates <- participant_data %>% 
  select(all_of(c(id_col,age_col,sex_col,group_col)))

# Standardise ID name for the anonymised analysis.
names(covariates)[names(covariates) == id_col] <- "ID"
feature_covariate_data <- left_join(feature_velocity,covariates,by = "ID")

# Rank-based inverse-normal transformation, scaled to approximately [-1, 1].
inverse_normal_scaled <- function(x) {
  non_missing <- !is.na(x)
  z <- rep(NA_real_,length(x))
  if (sum(non_missing) < 2) {return(z)}
  ranked_values <- rank(x[non_missing],ties.method = "average")
  transformed_values <- qnorm((ranked_values - 0.5) /sum(non_missing))
  max_abs <- max(abs(transformed_values),na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) {z[non_missing] <- 0} else {z[non_missing] <- transformed_values / max_abs}
  z}

# Identify velocity features while retaining ID and covariates unchanged.
exclude_cols <- c("ID",age_col,sex_col,group_col)
velocity_cols <- setdiff(names(feature_covariate_data),exclude_cols)

# For each velocity feature:
#   1) regress out age, sex, and study group;
#   2) extract residuals;
#   3) apply the inverse-normal transformation.
feature_normed <- lapply(velocity_cols,function(var) {
  model_formula <- reformulate(c(age_col,sex_col,group_col),response = var)
    fit <- lm(
      model_formula,
      data = feature_covariate_data,
      na.action = na.exclude
    )
    residual_values <- residuals(fit)
    inverse_normal_scaled(
      residual_values)}
)

feature_normed <- as.data.frame(feature_normed,check.names = FALSE)
names(feature_normed) <- velocity_cols

# Add participant ID back after transformation.
feature_normed <- cbind(ID = feature_covariate_data$ID,feature_normed)

# Split transformed features into WL and WM periods.
feature_normed_WL <- feature_normed %>%select(ends_with("_WL"))
feature_normed_WM <- feature_normed %>%select(ends_with("_WM"))

filter_correlated_features <- function(df,label,cutoff = 0.80) {
  # Remove zero-variance features before calculating correlations.
  non_zero_variance <- vapply(
    df,
    function(x) {feature_sd <- sd(x,na.rm = TRUE)
      is.finite(feature_sd) &&
        feature_sd > 0},
    logical(1))
  df <- df[ ,non_zero_variance,drop = FALSE]
  if (ncol(df) <= 1) {message(
      sprintf(
        "%s: correlation filtering skipped because <=1 feature remained.",
        label
      )
    )
    return(df)
  }
  cor_matrix <- cor(df,use = "pairwise.complete.obs")
  if (anyNA(cor_matrix)) {
    stop(sprintf(
        "%s: undefined correlations detected after removing zero-variance features.",
        label
      )
    )
  }
 high_cor_index <- findCorrelation(cor_matrix,cutoff = cutoff)
  message(
    sprintf(
      "%s: removing %d highly correlated feature(s)",
      label,
      length(high_cor_index)
    )
  )
  if (length(high_cor_index) > 0) {df <- df[ , -high_cor_index,drop = FALSE]
  }
  df
}


df_filtered_WL <- filter_correlated_features(feature_normed_WL,label = "WL",cutoff = correlation_cutoff)
df_filtered_WM <- filter_correlated_features(feature_normed_WM,label = "WM",cutoff = correlation_cutoff)

selected_features <- cbind(df_filtered_WL,df_filtered_WM)
selected_features <- cbind(ID = feature_normed$ID,selected_features)

outcome_data <- participant_data %>%
  select(all_of(c(id_col,outcome_col)))

names(outcome_data)[names(outcome_data) == id_col]      <- "ID"
names(outcome_data)[names(outcome_data) == outcome_col] <- "Outcome"

lasso_data <- left_join(selected_features,outcome_data,by = "ID")

y <- lasso_data$Outcome
x <- lasso_data %>%
  select(-ID,-Outcome) %>%
  as.matrix()

# alpha = 1 specifies LASSO regression.
# Setting nfolds equal to the number of observations implements LOOCV.

set.seed(lasso_seed)
cv_lasso <- cv.glmnet(
  x = x,
  y = y,
  alpha = 1,
  nfolds = nrow(x)
)

lambda_min_lasso <- cv_lasso$lambda.min

final_coef_lasso <- as.matrix(coef(cv_lasso,s = "lambda.min"))

selected_predictors_lasso <- data.frame(
  Predictor = rownames(final_coef_lasso),
  Coefficient = final_coef_lasso[, 1],
  row.names = NULL
) %>%
  filter(
    Predictor != "(Intercept)",
    Coefficient != 0
  )

lambda_index <- which.min(abs(cv_lasso$lambda -lambda_min_lasso))

mse_min_lasso <- cv_lasso$cvm[lambda_index]
rsq_min_lasso <- 1 -(mse_min_lasso /var(y,na.rm = TRUE))

selected_predictor_names <- selected_predictors_lasso$Predictor

lasso_selected_data <- lasso_data %>%
  select(ID,all_of( selected_predictor_names))


write.csv(selected_predictors_lasso,"lasso_selected_predictor_coefficients.csv",row.names = FALSE)
write.csv(lasso_selected_data,"lasso_selected_predictors.csv",row.names = FALSE)

