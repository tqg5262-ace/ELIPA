# =============================================================================
# Title:   Metabolomics/Anthropometric Velocity Pipeline & Random Forest
#          Classification of BMI Outcome Trajectory
# Purpose: 1) Compute rate-of-change ("velocity") features for metabolomic
#             and anthropometric variables across study time points.
#          2) Impute missing data, regress out covariates, and apply an
#             inverse-normal transformation to the resulting features.
#          3) Filter highly correlated features and fit a tuned random
#             forest classifier to predict BMI outcome trajectory
#             (gain / loss / stable).
#          4) Summarize and visualize variable importance and group
#             differences for top predictors.
#
# NOTE ON ANONYMIZATION: All local file paths, working directories, and
# source filenames have been replaced with generic placeholders. Update the
# "User-defined inputs" section below to point to your own files before
# running.
# =============================================================================


# Update these paths/parameters for your local environment before running.

covariate_file           <- "path/to/covariate_data.RData_or_csv"  # provides `data`: ID, age, sex, group, k2_0gluc, k4_0gluc, outcome, etc.
body_summary_file        <- "path/to/body_summary.xlsx"
metabolites_stage1_file  <- "path/to/metabolites_stage1.xlsx"
metabolites_stage2_file  <- "path/to/metabolites_stage2.xlsx"
metabolites_stage4_file  <- "path/to/metabolites_stage4.xlsx"
metabolites_stage1_sheet <- 1
metabolites_stage2_sheet <- 1


# NOTE: The object `data` (containing at minimum ID, age, sex, group,
# k2_0gluc, k4_0gluc, and the outcome variable) is referenced below but its
# creation step was not included in the code provided. Load/derive it here,
# e.g.:
#   data <- read_excel(covariate_file)
# before running the rest of the pipeline.


k1_metabolites <- read_excel(metabolites_stage1_file, sheet = metabolites_stage1_sheet)  # Stage 1 features
k2_metabolites <- read_excel(metabolites_stage2_file, sheet = metabolites_stage2_sheet)  # Stage 2 features
k4_metabolites <- read_excel(metabolites_stage4_file)                                    # Stage 4 features

ID <- k4_metabolites$ID

# Drop non-numeric ID/metadata column, coerce remaining columns to numeric
k1_metabolites <- as.data.frame(lapply(k1_metabolites[, -1], as.numeric))
k2_metabolites <- as.data.frame(lapply(k2_metabolites[, -1], as.numeric))
k4_metabolites <- as.data.frame(lapply(k4_metabolites[, -1], as.numeric))

# Harmonize column names across stages (strip stage-specific prefixes)
names(k1_metabolites) <- sub("^k1_", "", names(k1_metabolites))
names(k2_metabolites) <- sub("^k2_", "", names(k2_metabolites))


velocity_WL <- (k2_metabolites - k1_metabolites) / time_diff_WL  # Stage 1 -> Stage 2
velocity_WM <- (k4_metabolites - k2_metabolites) / time_diff_WM  # Stage 2 -> Stage 4

colnames(velocity_WL) <- paste0(colnames(velocity_WL), "_WL")
colnames(velocity_WM) <- paste0(colnames(velocity_WM), "_WM")

metabolites_velocity <- cbind(ID, velocity_WL, velocity_WM)


variable    <- data[, c(1, 3, 5, 7)]  # NOTE: column positions assumed stable; consider selecting by name for robustness
metabolites <- left_join(variable, metabolites_velocity, by = "ID")
ma          <- left_join(metabolites, body, by = "ID")

write.csv(ma, "combined_analysis_data.csv", row.names = FALSE)


set.seed(500)

# NOTE: columns 184:185 are dropped by position prior to imputation (assumed
# to be non-numeric/identifier columns unsuitable for missForest). Consider
# replacing with an explicit column-name-based selection for clarity and to
# guard against column-order changes upstream.
ma_for_imputation <- ma[, -c(184, 185)]

mf_ma <- missForest(ma_for_imputation, maxiter = 10, ntree = 100)
df_ma <- mf_ma$ximp

# For each metabolite/velocity feature, regress out age, sex, and group;
# apply an inverse-normal transformation to the residuals.
#
# NOTE: `inverse_normal_scaled()` is referenced but was not defined in the
# code provided. Define it before running, e.g.:
#   inverse_normal_scaled <- function(x) {
#     qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
#   }

exclude_cols <- c("ID", "age", "sex", "group")
ma_cols      <- setdiff(names(df_ma), exclude_cols)

ma_normed <- lapply(ma_cols, function(var) {
  model_formula <- as.formula(paste(var, "~ age + sex + group"))
  fit <- lm(model_formula, data = df_ma)
  res <- residuals(fit)
  inverse_normal_scaled(res)
})
ma_normed <- as.data.frame(ma_normed)
names(ma_normed) <- ma_cols

# NOTE: `unchanged` is intended to carry an identifier/outcome column through
# the transformation (originally "ID and Outcome" per the source comment),
# but only a single column (body[, 3]) is selected. Confirm this is the
# intended column before proceeding.
unchanged    <- body[, 3]
df_ma_normed <- cbind(unchanged, ma_normed)

write.csv(ma_normed, "metabolite_velocity_INT.csv", row.names = FALSE)

# Split into "_WL" and "_WM" feature sets, then drop highly correlated
# features (|r| > 0.8) within each set.

filter_correlated_features <- function(df, label, cutoff = 0.8) {
  cor_matrix <- cor(df, use = "complete.obs")
  corrplot(cor_matrix, method = "color", tl.cex = 0.3)
  high_cor_index <- findCorrelation(cor_matrix, cutoff = cutoff)
  message(sprintf("%s: removing %d highly correlated feature(s)", label, length(high_cor_index)))
  print(high_cor_index)df[, -high_cor_index]
}

ma_normed_WL <- ma_normed[, grep("_WL$", names(ma_normed))]
ma_normed_WM <- ma_normed[, grep("_WM$", names(ma_normed))]

df_ma_filtered_WL <- filter_correlated_features(ma_normed_WL, "WL")
df_ma_filtered_WM <- filter_correlated_features(ma_normed_WM, "WM")

write.csv(df_ma_filtered_WL, "metabolite_velocity_WL_filtered.csv", row.names = FALSE)
write.csv(df_ma_filtered_WM, "metabolite_velocity_WM_filtered.csv", row.names = FALSE)

rf_ma    <- cbind(df_ma_filtered_WL, df_ma_filtered_WM)
df_rf_ma <- cbind(unchanged, rf_ma)

names(df_rf_ma) <- sub("^velocity_", "", names(df_rf_ma))

df_rf_ma$pct_bmi_outcome_WM <- factor(
  df_rf_ma$pct_bmi_outcome_WM,
  levels = c("gain", "loss", "stable")
)

# Split into training (80%) and test (20%) sets
set.seed(500)
train_index <- createDataPartition(df_rf_ma$pct_bmi_outcome_WM, p = 0.8, list = FALSE)
train_set   <- df_rf_ma[train_index, ]
test_set    <- df_rf_ma[-train_index, ]

train_control <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
tune_grid     <- expand.grid(mtry = seq(2, ncol(df_rf_ma) - 1, by = 3))

set.seed(500)
rf_model <- train(
  pct_bmi_outcome_WM ~ .,
  data       = train_set,
  method     = "rf",
  trControl  = train_control,
  tuneGrid   = tune_grid,
  ntree      = 200,
  importance = TRUE
)

print(rf_model)
plot(rf_model)

predictions <- predict(rf_model, newdata = test_set)
conf_mat    <- confusionMatrix(predictions, test_set$pct_bmi_outcome_WM)
print(conf_mat)

overall_importance <- varImp(rf_model)
print(overall_importance)

imp_df <- overall_importance$importance %>%
  rownames_to_column(var = "Predictor")

plot_lollipop <- function(data, outcome, title, color) {
  top20 <- data %>%
    arrange(desc(.data[[outcome]])) %>%
    slice(1:20) %>%
    mutate(ImportanceLabel = round(.data[[outcome]], 2))
  
  ggplot(top20, aes(x = .data[[outcome]], y = reorder(Predictor, .data[[outcome]]))) +
    geom_segment(aes(x = 0, xend = .data[[outcome]], y = Predictor, yend = Predictor),
                 color = color) +
    geom_point(color = color, size = 3) +
    geom_text(aes(label = ImportanceLabel), hjust = -0.2, vjust = 0.5, color = "black", size = 2.5) +
    labs(title = title, x = "Importance", y = "Predictor") +
    theme_minimal() +
    theme(
      plot.margin = ggplot2::margin(1, 1, 1, 1, unit = "cm"),
      axis.text.y = element_text(size = 7)
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)))
}

plot_gain   <- plot_lollipop(imp_df, "gain",   "Top 20 Predictors for Gain",   "#D55E00")
plot_loss   <- plot_lollipop(imp_df, "loss",   "Top 20 Predictors for Loss",   "#009E73")
plot_stable <- plot_lollipop(imp_df, "stable", "Top 20 Predictors for Stable", "#56B4E9")

combined_plot <- arrangeGrob(plot_gain, plot_loss, plot_stable, ncol = 3)

ggsave(
  filename    = "RF_directional_top20_predictors.tiff",
  plot        = combined_plot,
  width       = 24,
  height      = 8,
  units       = "in",
  dpi         = 600,
  device      = "tiff",
  compression = "lzw"
)


importance_data <- varImp(rf_model)$importance

top_loss_vars <- importance_data %>%
  rownames_to_column(var = "Predictor") %>%
  arrange(desc(loss)) %>%
  slice(1:20) %>%
  pull(Predictor)

print(top_loss_vars)

vars_WL <- top_loss_vars[str_detect(top_loss_vars, "_WL$")]
vars_WM <- top_loss_vars[str_detect(top_loss_vars, "_WM$")]

df_plot <- df_rf_ma %>%
  rename(Class = pct_bmi_outcome_WM) %>%
  filter(Class != "stable") %>%
  select(Class, all_of(top_loss_vars)) %>%
  pivot_longer(cols = -Class, names_to = "Predictor", values_to = "Value") %>%
  mutate(Group = if_else(str_detect(Predictor, "_WL$"), "WL", "WM"))

plot_predictor_group <- function(df_plot, group_label, title) {
  df_plot %>%
    filter(Group == group_label) %>%
    ggplot(aes(x = Class, y = Value, fill = Class, color = Class)) +
    geom_boxplot(outlier.shape = NA, color = "black") +
    geom_jitter(width = 0.1, alpha = 0.5, shape = 21, color = "black") +
    facet_wrap(~ Predictor, scales = "free_y") +
    scale_fill_manual(values = c(gain = "#D55E00", loss = "#009E73")) +
    labs(title = title, x = "BMI Outcome Trajectory", y = "Scaled Value") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, hjust = 0.5),
      axis.title = element_text(size = 7.5),
      axis.text  = element_text(size = 10)
    )
}
