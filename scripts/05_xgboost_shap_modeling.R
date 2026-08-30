# =====================================================================
# Script 05: Non-Linear Machine Learning & SHAP Explanation
# Rationale: Applies Random Forest and XGBoost to decode non-linear 
# interactions, utilizing SHAP values to prioritize biomarker importance.
# =====================================================================
library(randomForest)
library(xgboost)
library(caret)
library(dplyr)
library(openxlsx)

data_list <- readRDS("../data/pca_linear_data.rds")
df <- data_list$df
predictors <- c("HSPG2", "GPC1", "SDC1", "TIMP4", "MMP9", "MMP2", "eNOS", "Fibronectin", "MMP9_TIMP4_ratio", "SDC1_eNOS_ratio", "Age", "Time_since_injury_months")
X_matrix <- as.matrix(df[, predictors])
y_vector <- df$NISS

# 1. Random Forest (with LOOCV validation)
set.seed(42)
model_rf <- randomForest(x = X_matrix, y = y_vector, ntree = 500, importance = TRUE)
train_control <- trainControl(method = "LOOCV")
rf_loocv <- train(X_matrix, y_vector, method = "rf", trControl = train_control, tuneLength = 1)

df_rf_importance <- data.frame(Predictor = rownames(importance(model_rf)),
                               IncMSE = as.numeric(importance(model_rf)[, "%IncMSE"]),
                               IncNodePurity = as.numeric(importance(model_rf)[, "IncNodePurity"])) %>% arrange(desc(IncMSE))
df_rf_diag <- data.frame(Metric = c("Mean of Squared Residuals (OOB)", "Variance Explained (OOB, %)", "RMSE (LOOCV)", "MAE (LOOCV)"),
                         Value = c(model_rf$mse[length(model_rf$mse)], model_rf$rsq[length(model_rf$rsq)] * 100, 
                                   rf_loocv$results$RMSE, rf_loocv$results$MAE))

# 2. Extreme Gradient Boosting (XGBoost) and SHAP Extraction
set.seed(42)
dtrain <- xgb.DMatrix(data = X_matrix, label = y_vector)
params <- list(booster = "gbtree", objective = "reg:squarederror", eta = 0.1, max_depth = 4)
model_xgb <- xgb.train(params = params, data = dtrain, nrounds = 100, verbose = 0)
xgb_pred <- predict(model_xgb, dtrain)
xgb_rmse <- sqrt(mean((y_vector - xgb_pred)^2))

shap_matrix <- predict(model_xgb, newdata = dtrain, predcontrib = TRUE)
if ("BIAS" %in% colnames(shap_matrix)) {
  shap_features <- shap_matrix[, colnames(shap_matrix) != "BIAS"]
} else {
  shap_features <- shap_matrix[, -ncol(shap_matrix)]
}
shap_means <- colMeans(abs(shap_features))
df_xgb_shap <- data.frame(Predictor = names(shap_means), Mean_Absolute_SHAP = as.numeric(shap_means)) %>% arrange(desc(Mean_Absolute_SHAP))
df_xgb_diag <- data.frame(Metric = c("Train RMSE"), Value = c(xgb_rmse))

# 3. Export full pipeline results to Excel Table
wb <- createWorkbook()
sheets <- c("Desc_NISS_Scale", "Stats_FDR_NISS_Scale", "Desc_Endotypes", "Stats_FDR_Endotypes",
            "PCA_Variance_Diag", "PCA_Loadings", "Robust_Regression_Coefs", "Robust_Diag",
            "LASSO_Coefs", "LASSO_Diag", "Random_Forest_Importance", "Random_Forest_Diag",
            "XGBoost_SHAP", "XGBoost_Diag")

for (s in sheets) addWorksheet(wb, s)

# Write data from previous steps
writeData(wb, "Desc_NISS_Scale", data_list$table)
writeData(wb, "Stats_FDR_NISS_Scale", data_list$stats)
writeData(wb, "Desc_Endotypes", data_list$table_endo)
writeData(wb, "Stats_FDR_Endotypes", data_list$stats_endo)
writeData(wb, "PCA_Variance_Diag", bind_rows(data_list$linear_exports$pca_diag, data.frame(Metric="", Value=NA), 
                                             data_list$linear_exports$pca_var %>% rename(Metric=Component, Value=Eigenvalue)))
writeData(wb, "PCA_Loadings", data_list$linear_exports$pca_loadings)
writeData(wb, "Robust_Regression_Coefs", data_list$linear_exports$robust_coef)
writeData(wb, "Robust_Diag", data_list$linear_exports$robust_diag)
writeData(wb, "LASSO_Coefs", data_list$linear_exports$lasso_coef)
writeData(wb, "LASSO_Diag", data_list$linear_exports$lasso_diag)

# Write ML data
writeData(wb, "Random_Forest_Importance", df_rf_importance)
writeData(wb, "Random_Forest_Diag", df_rf_diag)
writeData(wb, "XGBoost_SHAP", df_xgb_shap)
writeData(wb, "XGBoost_Diag", df_xgb_diag)

saveWorkbook(wb, "../output/tables/Final_Analysis_Complete.xlsx", overwrite = TRUE)

data_list$df_xgb_shap <- df_xgb_shap
saveRDS(data_list, "../data/final_workspace.rds")