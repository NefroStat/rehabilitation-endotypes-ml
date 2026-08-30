# =====================================================================
# Script 06: Data Visualization (Publication Quality Formatting)
# Rationale: Generates Figures 1-4 formatted for Q1 journal submission,
# ensuring English labels, robust color palettes, and strict 300 DPI limits.
# =====================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales) 
library(factoextra)
library(patchwork)

# Define universal output directory
out_dir <- "../output/figures/"
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

data_list <- readRDS("../data/final_workspace.rds")
df <- data_list$df
df_volcano <- data_list$df_volcano
df_xgb_shap <- data_list$df_xgb_shap

# ---------------------------------------------------------------------
# Figure 1: Diverging Lollipop Chart
# ---------------------------------------------------------------------
# Використовуємо $Marker замість $Predictor
df_volcano$Predictor_Clean <- recode(df_volcano$Marker,
                                     "SDC1_eNOS_ratio" = "SDC1/eNOS ratio",
                                     "Time_since_injury_months" = "Time since injury, months",
                                     "Age" = "Age, years",
                                     "MMP9_TIMP4_ratio" = "MMP9/TIMP4 ratio")

df_volcano <- df_volcano %>% arrange(Log2FC)
df_volcano$Predictor_Clean <- factor(df_volcano$Predictor_Clean, levels = df_volcano$Predictor_Clean)

fig1 <- ggplot(df_volcano, aes(x = Predictor_Clean, y = Log2FC, color = Significance)) +
  geom_segment(aes(x = Predictor_Clean, xend = Predictor_Clean, y = 0, yend = Log2FC), size = 1.2) +
  geom_point(size = 5) +
  coord_flip() +
  scale_y_continuous(breaks = c(-2, -1, 0, 1, 2, 3), labels = c("0.25", "0.5", "1", "2", "4", "8")) +
  scale_color_manual(values = c("Downregulated in Endotype 2" = "#2b83ba", 
                                "Non-significant (FDR > 0.05)" = "gray70", 
                                "Upregulated in Endotype 2" = "#d7191c")) +
  theme_classic() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_blank(),
        axis.text.y = element_text(size = 11, color = "black"), axis.text.x = element_text(size = 11, color = "black")) +
  labs(y = "Fold Change (Endotype 2 vs Endotype 1)")

ggsave(paste0(out_dir, "Figure_1_Lollipop.tiff"), plot = fig1, width = 8, height = 6, dpi = 300, compression = "lzw")
# ---------------------------------------------------------------------
# Figure 2: Boxplots with Jitter (Log Scale)
# ---------------------------------------------------------------------
key_markers <- c("SDC1", "TIMP4", "Fibronectin", "SDC1_eNOS_ratio")
df_long <- df %>% dplyr::select(all_of(key_markers), Group) %>%
  pivot_longer(cols = all_of(key_markers), names_to = "Marker", values_to = "Value")

df_long$Group_Eng <- recode(df_long$Group,
                            "Endotype 1 (Stabilized)" = "Endotype 1\n(Stabilized)",
                            "Endotype 2 (Hyper-remodeling)" = "Endotype 2\n(Hyper-remodeling)")
df_long$Marker_Clean <- factor(recode(df_long$Marker, "SDC1_eNOS_ratio" = "SDC1/eNOS ratio"), levels = c("SDC1", "TIMP4", "Fibronectin", "SDC1/eNOS ratio"))
endo_colors <- c("Endotype 1\n(Stabilized)" = "#2b83ba", "Endotype 2\n(Hyper-remodeling)" = "#d7191c")

fig2 <- ggplot(df_long, aes(x = Group_Eng, y = Value, fill = Group_Eng)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.6) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, color = "black") +
  facet_wrap(~ Marker_Clean, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = endo_colors) +
  scale_y_log10(labels = label_number(drop0trailing = TRUE, big.mark = "")) +
  theme_classic() +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(), axis.text = element_text(size = 10, color = "black"), axis.title.y = element_text(size = 12, face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5)) +
  labs(y = "Concentration / Ratio")

ggsave(paste0(out_dir, "Figure_2_Boxplots.tiff"), plot = fig2, width = 8, height = 7, dpi = 300, compression = "lzw")

# ---------------------------------------------------------------------
# Figure 3: PCA Biplot with Manual Annotation
# ---------------------------------------------------------------------
pca_custom <- data_list$pca_model
rownames(pca_custom$rotation) <- recode(rownames(pca_custom$rotation), "SDC1_eNOS_ratio" = "SDC1/eNOS ratio", "MMP9_TIMP4_ratio" = "MMP9/TIMP4 ratio", "Age" = "Age", "Time_since_injury_months" = "Time since injury")

# Mask specific overlapping labels using spaces (leaves the arrow intact)
rownames(pca_custom$rotation)[rownames(pca_custom$rotation) == "SDC1/eNOS ratio"] <- " "
rownames(pca_custom$rotation)[rownames(pca_custom$rotation) == "Fibronectin"] <- "  "
rownames(pca_custom$rotation)[rownames(pca_custom$rotation) == "MMP2"] <- "   "
rownames(pca_custom$rotation)[rownames(pca_custom$rotation) == "eNOS"] <- "    "
rownames(pca_custom$rotation)[rownames(pca_custom$rotation) == "MMP9"] <- "     "

fig3_base <- fviz_pca_biplot(pca_custom, geom.ind = "point", col.ind = df$Group, palette = c("#2b83ba", "#d7191c"),
                             addEllipses = TRUE, ellipse.type = "confidence", col.var = "black", repel = TRUE, 
                             labelsize = 4.5, arrowsize = 0.6, legend.title = "Endotype", title = "") +           
  theme_classic() + theme(legend.position = "top", legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = "bold"))

fig3 <- fig3_base +
  annotate("text", x = 1.0, y = 3.7, label = "SDC1/eNOS ratio", size = 4.5) +
  annotate("text", x = 2.5, y = 2.7, label = "Fibronectin", size = 4.5) +
  annotate("text", x = 1.7, y = 1.1, label = "MMP2", size = 4.5) +
  annotate("text", x = -0.45, y = -1.1, label = "eNOS", size = 4.5) + 
  annotate("text", x = -0.2, y = -3.4, label = "MMP9", size = 4.5)

ggsave(paste0(out_dir, "Figure_3_PCA.tiff"), plot = fig3, width = 8, height = 7, dpi = 300, compression = "lzw")

# ---------------------------------------------------------------------
# Figure 4: SHAP + NISS Boxplot (Patchwork)
# ---------------------------------------------------------------------
df_xgb_shap$Predictor_Clean <- factor(recode(df_xgb_shap$Predictor, "SDC1_eNOS_ratio" = "SDC1/eNOS ratio", "MMP9_TIMP4_ratio" = "MMP9/TIMP4 ratio", "Time_since_injury_months" = "Time since injury, months", "Age" = "Age, years"), levels = rev(recode(df_xgb_shap$Predictor, "SDC1_eNOS_ratio" = "SDC1/eNOS ratio", "MMP9_TIMP4_ratio" = "MMP9/TIMP4 ratio", "Time_since_injury_months" = "Time since injury, months", "Age" = "Age, years")))

plot_shap <- ggplot(df_xgb_shap, aes(x = Predictor_Clean, y = Mean_Absolute_SHAP)) +
  geom_bar(stat = "identity", fill = "#7fc97f", color = "black", width = 0.7) +
  coord_flip() + theme_classic() + labs(x = "Predictors", y = "Mean Absolute SHAP Value (Impact on NISS)")

df_niss <- df %>% mutate(Group_Eng = factor(recode(Group, "Endotype 1 (Stabilized)" = "Endotype 1\n(Stabilized)", "Endotype 2 (Hyper-remodeling)" = "Endotype 2\n(Hyper-remodeling)"), levels = c("Endotype 1\n(Stabilized)", "Endotype 2\n(Hyper-remodeling)")))

plot_niss <- ggplot(df_niss, aes(x = Group_Eng, y = NISS, fill = Group_Eng)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black", width = 0.5) +
  geom_jitter(width = 0.15, size = 1.5, alpha = 0.5, color = "black") +
  scale_fill_manual(values = endo_colors) + theme_classic() + theme(legend.position = "none", axis.title.x = element_blank())

fig4 <- plot_shap + plot_niss + plot_annotation(tag_levels = 'A')
ggsave(paste0(out_dir, "Figure_4_SHAP_NISS.tiff"), plot = fig4, width = 11, height = 5, dpi = 300, compression = "lzw")