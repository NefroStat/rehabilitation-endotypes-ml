# =====================================================================
# Script 04: Principal Component Analysis (PCA) & Linear Modeling
# Rationale: Evaluates multidimensional biomarker topology and demonstrates 
# the limitations of standard linear models (LASSO, Robust Regression) 
# in capturing post-traumatic complexity.
# =====================================================================
library(psych)
library(factoextra)
library(robustbase)
library(glmnet)
library(dplyr)
library(openxlsx)

data_list <- readRDS("../data/clustered_data.rds")
df <- data_list$df
predictors <- c("HSPG2", "GPC1", "SDC1", "TIMP4", "MMP9", "MMP2", "eNOS", "Fibronectin", "MMP9_TIMP4_ratio", "SDC1_eNOS_ratio", "Age", "Time_since_injury_months")

X_matrix <- as.matrix(df[, predictors])
y_vector <- df$NISS

# 1. PCA and KMO Diagnostics
cor_matrix <- cor(X_matrix)
kmo_res <- KMO(cor_matrix)
bartlett_res <- cortest.bartlett(cor_matrix, n = nrow(X_matrix))

df_pca_diag <- data.frame(
  Metric = c("KMO (Overall MSA)", "Bartlett Test (Chi-Square)", "Bartlett Test (p-value)"),
  Value = c(kmo_res$MSA, bartlett_res$chisq, bartlett_res$p.value)
)

pca_model <- prcomp(X_matrix, center = TRUE, scale. = TRUE)
pca_var <- get_eigenvalue(pca_model)
df_pca_var <- data.frame(
  Component = rownames(pca_var),
  Eigenvalue = pca_var$eigenvalue,
  Variance_Percent = pca_var$variance.percent,
  Cumulative_Variance_Percent = pca_var$cumulative.variance.percent
)

df_pca_loadings <- as.data.frame(pca_model$rotation)
df_pca_loadings$Predictor <- rownames(df_pca_loadings)
df_pca_loadings <- df_pca_loadings %>% dplyr::select(Predictor, everything())

# 2. Linear Modeling (Robust Regression)
formula_robust <- as.formula(paste("NISS ~", paste(predictors, collapse = " + ")))
model_robust <- lmrob(formula_robust, data = df)
summary_robust <- summary(model_robust)
df_robust_coef <- data.frame(
  Predictor = rownames(summary_robust$coefficients),
  Estimate = summary_robust$coefficients[, "Estimate"],
  Std_Error = summary_robust$coefficients[, "Std. Error"],
  t_value = summary_robust$coefficients[, "t value"],
  p_value = summary_robust$coefficients[, "Pr(>|t|)"]
)
df_robust_diag <- data.frame(Metric = c("Robust R-squared", "Adjusted Robust R-squared", "Residual Standard Error"),
                             Value = c(summary_robust$r.squared, summary_robust$adj.r.squared, summary_robust$sigma))

# 3. L1-Regularization (LASSO)
set.seed(42)
cv_lasso <- cv.glmnet(X_matrix, y_vector, alpha = 1, nfolds = 5)
model_lasso <- glmnet(X_matrix, y_vector, alpha = 1, lambda = cv_lasso$lambda.min)
lasso_coefs <- as.matrix(coef(model_lasso))
df_lasso_coef <- data.frame(Predictor = rownames(lasso_coefs), Coefficient = as.numeric(lasso_coefs[, 1])) %>% filter(Coefficient != 0)
df_lasso_diag <- data.frame(Metric = c("Optimal Lambda (min)", "Lambda (1 standard error)", "Deviance Ratio (Explained Variance)"),
                            Value = c(cv_lasso$lambda.min, cv_lasso$lambda.1se, model_lasso$dev.ratio))

data_list$pca_model <- pca_model
data_list$linear_exports <- list(pca_diag = df_pca_diag, pca_var = df_pca_var, pca_loadings = df_pca_loadings, 
                                 robust_coef = df_robust_coef, robust_diag = df_robust_diag, 
                                 lasso_coef = df_lasso_coef, lasso_diag = df_lasso_diag)
saveRDS(data_list, "../data/pca_linear_data.rds")