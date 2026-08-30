# =====================================================================
# Script 01: Data Loading and Preprocessing
# Rationale: Cleans raw clinical data, calculates derived pathophysiological 
# indices, and prepares the dataset for downstream statistical and ML analysis.
# =====================================================================
library(dplyr)

# 1. Load data using relative paths (ensures reproducibility across different machines)
file_path <- "../data/markers_results_comma.csv"
df <- read.csv2(file_path, stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

# Helper function to safely convert comma-separated decimals to numeric
clean_numeric <- function(x) {
  as.numeric(gsub(",", ".", as.character(x)))
}

# 2. Exclude specific outlier (ID 678528) and somatic controls (ISS = 0)
df <- df %>% filter(ID != 678528 & clean_numeric(ISS) > 0)

# 3. Format clinical and demographic variables
df$Age <- clean_numeric(df$Вік)
df$ISS <- clean_numeric(df$ISS)
df$NISS <- clean_numeric(df$NISS)

df$Date <- as.Date(df$Дата, format = "%d.%m.%Y")
df$Date_of_injury <- as.Date(df$Дата.травми, format = "%d.%m.%Y") 
df$Time_since_injury_months <- as.numeric(difftime(df$Date, df$Date_of_injury, units = "days")) / 30.44

# 4. Format raw biomarker concentrations
markers_raw <- c("HSPG2", "GPC1", "SDC1", "TIMP4", "MMP9", "MMP2", "eNOS", "Fibronectin")
df[markers_raw] <- lapply(df[markers_raw], clean_numeric)

# 5. Calculate pathophysiological indices (ratios) with a small constant to prevent division by zero
df$MMP9_TIMP4_ratio <- df$MMP9 / (df$TIMP4 + 1e-6) 
df$SDC1_eNOS_ratio <- df$SDC1 / (df$eNOS + 1e-6)

markers <- c(markers_raw, "MMP9_TIMP4_ratio", "SDC1_eNOS_ratio")
predictors <- c(markers, "Age", "Time_since_injury_months")

# 6. Remove missing values to ensure matrix integrity for ML algorithms
df <- df %>% na.omit()

# Save preprocessed environment for next scripts
saveRDS(df, "../data/preprocessed_data.rds")