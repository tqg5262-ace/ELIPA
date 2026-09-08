# =============================================================================
# Title:   Calculation of Feature Velocity Between Two Time Points
# Purpose: Computes the rate of change ("velocity") for a set of numeric
#          features between two measurement time points (Stage 1 and Stage 2),
#          using a fixed time interval.
# =============================================================================


library(readxl)

# NOTE: Update these paths/parameters for your local environment before running.
# File paths, sheet indices, and directory structure have been anonymized/
# generalized for code sharing and peer review.

stage1_file  <- "path/to/stage1_data.xlsx"   # Stage 1 (baseline) feature file
stage2_file  <- "path/to/stage2_data.xlsx"   # Stage 2 (follow-up) feature file
stage1_sheet <- 1                            # Sheet index/name for Stage 1 data
stage2_sheet <- 1                            # Sheet index/name for Stage 2 data

time_diff_weeks <- stage_duration / total_duration   # Time interval between Stage 1 and Stage 2, in weeks

output_file <- "velocity_data.csv"

stage1_data <- read_excel(stage1_file, sheet = stage1_sheet)
stage2_data <- read_excel(stage2_file, sheet = stage2_sheet)

# Drop the first row if it contains non-numeric metadata (e.g., units, labels),
# then coerce all remaining columns to numeric.

clean_numeric_df <- function(df) {
  df_clean <- df[-1, ]
  df_clean <- as.data.frame(lapply(df_clean, as.numeric))
  return(df_clean)
}

stage1_clean <- clean_numeric_df(stage1_data)
stage2_clean <- clean_numeric_df(stage2_data)

# Verify that both datasets share identical column structure before proceeding
if (!identical(colnames(stage1_clean), colnames(stage2_clean))) {
  stop("Column names in Stage 1 and Stage 2 data do not match. Please check the input files.")
}

# Velocity = (Stage 2 value - Stage 1 value) / time interval
velocity <- (stage2_clean - stage1_clean) / time_diff_weeks

rownames(velocity) <- paste0("Velocity_", seq_len(nrow(velocity)))
colnames(velocity) <- colnames(stage1_clean)

write.csv(velocity, output_file, row.names = TRUE)
