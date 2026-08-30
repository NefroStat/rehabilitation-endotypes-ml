# =====================================================================
# Script 02: Statistical Analysis by Anatomical Severity (NISS)
# Rationale: Stratifies the cohort based on physical injury volume and 
# evaluates baseline biomarker differences using non-parametric testing.
# =====================================================================
library(dplyr)
library(gtsummary)
library(pwr)

df <- readRDS("../data/preprocessed_data.rds")
predictors <- c("HSPG2", "GPC1", "SDC1", "TIMP4", "MMP9", "MMP2", "eNOS", "Fibronectin", "MMP9_TIMP4_ratio", "SDC1_eNOS_ratio", "Age", "Time_since_injury_months")

# 1. Stratify by NISS Score
df$NISS_Category <- ifelse(df$NISS >= 25, "Severe polytrauma (NISS >= 25)", "Moderate trauma (NISS < 25)")
df$NISS_Category <- factor(df$NISS_Category, levels = c("Moderate trauma (NISS < 25)", "Severe polytrauma (NISS >= 25)"))

# 2. Generate descriptive statistics table
theme_gtsummary_language(language = "en", decimal.mark = ".", big.mark = " ")

table_niss <- df %>%
  dplyr::select(all_of(predictors), NISS_Category) %>%
  tbl_summary(
    by = NISS_Category,
    statistic = all_continuous() ~ "{median} [{p25}; {p75}]",
    missing = "no"
  ) %>%
  add_p(test = all_continuous() ~ "wilcox.test") %>%
  modify_header(label = "**Parameter**")

res_table_niss <- as_tibble(table_niss)

# 3. Calculate Effect Sizes (Cohen's d) and Statistical Power
results_list_niss <- list()
for (m in predictors) {
  w_test <- wilcox.test(df[[m]] ~ df$NISS_Category, exact = FALSE)
  p_mw <- w_test$p.value
  p_calc <- ifelse(p_mw == 0, 1e-16, p_mw) # Prevent infinite Z values
  
  n1 <- sum(df$NISS_Category == "Moderate trauma (NISS < 25)")
  n2 <- sum(df$NISS_Category == "Severe polytrauma (NISS >= 25)")
  
  Z <- qnorm(p_calc / 2, lower.tail = FALSE)
  r_eff <- Z / sqrt(n1 + n2)
  d_eff <- 2 * r_eff / sqrt(1 - r_eff^2)
  d_eff <- ifelse(is.na(d_eff) | is.infinite(d_eff), 0.01, d_eff)
  
  # Calculate 95% Confidence Intervals for Cohen's d
  var_d <- ((n1 + n2) / (n1 * n2)) + (d_eff^2 / (2 * (n1 + n2)))
  se_d <- sqrt(var_d)
  ci_lower <- d_eff - 1.96 * se_d
  ci_upper <- d_eff + 1.96 * se_d
  
  power_mw <- pwr.t2n.test(n1 = n1, n2 = n2, d = d_eff, sig.level = 0.05)$power
  
  results_list_niss[[m]] <- data.frame(Marker = m, P_Value = p_mw, Cohen_d = d_eff, 
                                       CI_95_Lower = ci_lower, CI_95_Upper = ci_upper, Power = power_mw)
}

res_niss_stats <- bind_rows(results_list_niss)
# Apply Benjamini-Hochberg correction for multiple testing
res_niss_stats$FDR <- p.adjust(res_niss_stats$P_Value, method = "BH")
res_niss_stats <- res_niss_stats %>% 
  dplyr::select(Marker, P_Value, FDR, Cohen_d, CI_95_Lower, CI_95_Upper, Power) %>% 
  mutate_if(is.numeric, round, 4)

saveRDS(list(df = df, table = res_table_niss, stats = res_niss_stats), "../data/niss_stats.rds")