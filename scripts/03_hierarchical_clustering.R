# =====================================================================
# Script 03: Molecular Endotyping (Hierarchical Clustering)
# Rationale: Utilizes unsupervised agglomerative clustering to identify 
# distinct pathophysiological trajectories independent of physical injury.
# =====================================================================
library(cluster)
library(dplyr)
library(gtsummary)
library(pwr)

data_list <- readRDS("../data/niss_stats.rds")
df <- data_list$df
markers <- c("HSPG2", "GPC1", "SDC1", "TIMP4", "MMP9", "MMP2", "eNOS", "Fibronectin", "MMP9_TIMP4_ratio", "SDC1_eNOS_ratio")
predictors <- c(markers, "Age", "Time_since_injury_months")

# 1. Z-score normalization and Ward's clustering
df_scaled <- scale(df[, markers]) 
hc <- hclust(dist(df_scaled), method = "ward.D2")
df$Cluster_Raw <- cutree(hc, k = 2)

# Dynamically assign the hyper-remodeling label to the cluster with higher Fibronectin
med_fibro <- tapply(df$Fibronectin, df$Cluster_Raw, median)
hyper_cluster_id <- as.numeric(names(which.max(med_fibro)))

df$Group <- ifelse(df$Cluster_Raw == hyper_cluster_id, 
                   "Endotype 2 (Hyper-remodeling)", 
                   "Endotype 1 (Stabilized)")
df$Group <- factor(df$Group, levels = c("Endotype 1 (Stabilized)", "Endotype 2 (Hyper-remodeling)"))

# 2. Endotype Descriptive Statistics
table_endotypes <- df %>%
  dplyr::select(all_of(c(predictors, "NISS")), Group) %>%
  tbl_summary(
    by = Group,
    statistic = all_continuous() ~ "{median} [{p25}; {p75}]",
    missing = "no"
  ) %>%
  add_p(test = all_continuous() ~ "wilcox.test") %>%
  modify_header(label = "**Parameter**")

res_table_endotypes <- as_tibble(table_endotypes)

# 3. Calculate Effect Sizes and Statistical Power for Endotypes
results_list_endo <- list()
for (m in c(predictors, "NISS")) {
  w_test <- wilcox.test(df[[m]] ~ df$Group, exact = FALSE)
  p_mw <- w_test$p.value
  p_calc <- ifelse(p_mw == 0, 1e-16, p_mw)
  
  n1 <- sum(df$Group == "Endotype 1 (Stabilized)")
  n2 <- sum(df$Group == "Endotype 2 (Hyper-remodeling)")
  
  Z <- qnorm(p_calc / 2, lower.tail = FALSE)
  r_eff <- Z / sqrt(n1 + n2)
  d_eff <- 2 * r_eff / sqrt(1 - r_eff^2)
  d_eff <- ifelse(is.na(d_eff) | is.infinite(d_eff), 0.01, d_eff)
  
  var_d <- ((n1 + n2) / (n1 * n2)) + (d_eff^2 / (2 * (n1 + n2)))
  se_d <- sqrt(var_d)
  ci_lower <- d_eff - 1.96 * se_d
  ci_upper <- d_eff + 1.96 * se_d
  
  power_mw <- pwr.t2n.test(n1 = n1, n2 = n2, d = d_eff, sig.level = 0.05)$power
  
  results_list_endo[[m]] <- data.frame(Marker = m, P_Value = p_mw, Cohen_d = d_eff, 
                                       CI_95_Lower = ci_lower, CI_95_Upper = ci_upper, Power = power_mw)
}

res_endo_stats <- bind_rows(results_list_endo)
res_endo_stats$FDR <- p.adjust(res_endo_stats$P_Value, method = "BH")
res_endo_stats <- res_endo_stats %>% 
  dplyr::select(Marker, P_Value, FDR, Cohen_d, CI_95_Lower, CI_95_Upper, Power) %>% 
  mutate_if(is.numeric, round, 4)

# 4. Generate Volcano df for Visualization (Log2 Fold Change)
medians_endo1 <- apply(df[df$Group == "Endotype 1 (Stabilized)", predictors], 2, median, na.rm=TRUE)
medians_endo2 <- apply(df[df$Group == "Endotype 2 (Hyper-remodeling)", predictors], 2, median, na.rm=TRUE)
fold_change <- medians_endo2 / medians_endo1

df_volcano <- res_endo_stats %>%
  filter(Marker %in% predictors) %>%
  mutate(Log2FC = log2(fold_change[Marker]),
         Significance = case_when(
           FDR < 0.05 & Log2FC > 0 ~ "Upregulated in Endotype 2",
           FDR < 0.05 & Log2FC < 0 ~ "Downregulated in Endotype 2",
           TRUE ~ "Non-significant (FDR > 0.05)"
         ))

data_list$df <- df
data_list$table_endo <- res_table_endotypes
data_list$stats_endo <- res_endo_stats
data_list$df_volcano <- df_volcano
saveRDS(data_list, "../data/clustered_data.rds")