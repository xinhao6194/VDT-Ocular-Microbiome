# 清除工作空间中的所有对象
rm(list = ls())
# 设置工作路径
setwd("/Users/xinhao/Desktop/")

# 加载包
library(tidyverse)
library(pheatmap)
library(corrplot)
library(Hmisc)
library(ggpubr)
library(DESeq2)
library(readxl)
library(vegan)

# 1) 安装 Bioconductor（推荐 R ≥ 4.2）
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")   # 2025‑最新 Bioconductor 版本

# 2) 安装 MaAsLin2（当前最新 1.23.0，2025‑09‑11）
BiocManager::install("Maaslin2") # [[1]]



# 3) 加载包
library(Maaslin2)

# 设定文件路径
feat_file <- "2.csv"
meta_file <- "1.csv"
out_dir   <- "maaslin2_output"
colnames(clinical_data)
# 读取文件（根据实际分隔符修改 sep）
feat <- read.csv("3.csv", header = TRUE, sep = ",", row.names = 1)
meta <- read.csv("1.csv", header = TRUE, sep = ",", row.names = 1)

# 调用函数
fit <- Maaslin2(
  input_data      = feat,
  input_metadata  = meta,
  output          = "output_dir_tax",
  
  # ---------- 数据预处理 ----------
  min_prevalence  = 0.1,          # 过滤出现频率 <10% 的通路
  normalization   = "TSS",        # 总计数归一化
  transform       = "LOG",        # log2(x+1) 变换
  standardize     = TRUE,        # 连续协变量 Z‑score
  
  # ---------- 模型设定 ----------
  analysis_method = "LM",         # 线性模型（对数化后适用）
  random_effects  = "subject_id", # 同一受试者的重复测量作为随机效应
  
  # --------- fixed_effects（多个临床变量） ----------
  # 这里假设在 metadata.tsv 中已经把列名中的空格改为下划线
  fixed_effects = c("TBUT", "OSDI", "CFS", "Meibomian_Gland_Scoring"),
  
  # ---------- 多重检验 ----------
  correction      = "BH",        # Benjamini‑Hochberg FDR 校正
  
  # ---------- 可视化 ----------
  plot_heatmap    = FALSE,
  heatmap_first_n = 50,
  plot_scatter    = FALSE,
  
  # ---------- 计算资源 ----------
  cores           = 4           # 根据机器核数自行调整
)

vars <- c("TBUT", "OSDI", "CFS", "Meibomian_Gland_Scoring")

all_res <- list()

for (v in vars) {
  fit <- Maaslin2(
    input_data = feat,
    input_metadata = meta,
    output = paste0("output_dir_tax_", v),
    min_prevalence = 0.1,
    normalization = "TSS",
    transform = "LOG",
    standardize = TRUE,
    analysis_method = "LM",
    random_effects = "subject_id",
    fixed_effects = v,
    correction = "BH",
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    cores = 4
  )
  
  res <- fit$results
  res$clinical_metric <- v
  all_res[[v]] <- res
}

final_res <- bind_rows(all_res)

# 全局 BH 校正，基于所有模型的 pval
final_res$qval_global <- p.adjust(final_res$pval, method = "BH")

# 可选：按 qval_global 排序
final_res <- final_res %>% arrange(qval_global)

# 查看显著结果
sig_res <- final_res %>% filter(qval_global < 0.25)

# 保存
write.csv(final_res, "output_dir_tax/maaslin2_single_predictor_globalBH.csv", row.names = FALSE)
library(pheatmap)
library(viridis)

# 准备热图数据（前30个显著特征）
sig_results <- fit$results %>%
  filter(qval < 0.05) %>%
  arrange(qval) %>%
  head(30)


# 提取每个协变量的效应量
# 注意：这里假设每个特征对多个协变量都有结果
# 您需要根据实际数据结构调整
heatmap_data <- fit$results %>%
  filter(feature %in% sig_results$feature) %>%
  select(feature, metadata, coef) %>%
  pivot_wider(names_from = metadata, values_from = coef) %>%
  column_to_rownames("feature")

# 绘制热图
pheatmap(
  fit$results,
  color = colorRampPalette(c("#2E74B5", "white", "#B33D3D"))(100),
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Effect Size Heatmap",
  filename = "output_dir_tax/effect_heatmap.pdf",
  width = 8,
  height = 10
)

# 1. 加载所需的包
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

# 2. 读取数据 (请确保文件路径正确)
df <- read.csv("maaslin2_single_predictor_globalBH.csv")

# 3. 数据预处理
# - 过滤掉 "Others" (如果有的话)
# - 修改临床指标的名称，使其在图上更好看
# - 根据 pval 或 qval 生成显著性星号 (这里演示基于 pval < 0.05 标星号)
#   如果您想基于 qval < 0.25 标星号，请修改下面 sig_label 的判断条件
df_clean <- df %>%
  filter(feature != "Others") %>%
  mutate(
    clinical_metric = recode(clinical_metric, 
                             "OSDI" = "OSDI Score",
                             "Meibomian_Gland_Scoring" = "Meibomian Gland Score",
                             "CFS" = "CFS",
                             "TBUT" = "TBUT"),
    # 设定显著性标签
    sig_label = case_when(
      qval < 0.05 ~ "**",
      qval < 0.25 ~ "*",
      TRUE ~ ""
    ),
    # 将 coef 格式化为两位小数，用于格子内显示
    coef_text = sprintf("%.2f", coef)
  )
# 4. 补全缺失的网格 (让热图呈完美的矩形，没有空白)
# 创建一个所有 feature 和 clinical_metric 的组合
all_features <- unique(df_clean$feature)
all_metrics <- unique(df_clean$clinical_metric)
grid <- expand.grid(feature = all_features, clinical_metric = all_metrics)

df_plot <- grid %>%
  left_join(df_clean, by = c("feature", "clinical_metric")) %>%
  # 将 NA 替换为 0（或者你想留白也可以把 0 改为 NA）
  mutate(coef = ifelse(is.na(coef), 0, coef),
         sig_label = ifelse(is.na(sig_label), "", sig_label),
         coef_text = ifelse(is.na(coef_text), "", coef_text))

# 5. 排序：让特征和指标有序排列，提升美观度
# 可以自定义顺序，这里默认按字母排序，指标按指定顺序
df_plot$clinical_metric <- factor(df_plot$clinical_metric, 
                                  levels = c("OSDI Score", "TBUT", "CFS", "Meibomian Gland Score"))

# 6. 绘制热图
p <- ggplot(df_plot, aes(x = clinical_metric, y = feature, fill = coef)) +
  # 【修改这里】：width 和 height 控制格子大小（1 表示撑满，0.9 表示留白）
  geom_tile(color = "white", linewidth = 0.5, width = 0.97, height = 0.95) + 
  geom_text(aes(label = sig_label), size = 5, color = "black", vjust = 0.75) +
  scale_fill_gradient2(low = "#4575b4", mid = "white", high = "#d73027", 
                       midpoint = 0, name = "Coefficient") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11, 
                               face = "bold",color = "black", family = "serif"),
    axis.text.y = element_text(size = 10,color = "black", 
                               family = "serif"),
    axis.title = element_blank(), 
    panel.grid = element_blank(), 
    legend.position = "right"
  ) +
  coord_fixed(ratio = 0.6) 

print(p)
# 7. 显示并保存图片
print(p)
ggsave("Microbiome_Clinical_Heatmap.pdf", plot = p, width = 8, height = 7, dpi = 300)
ggsave("Microbiome_Clinical_Heatmap.png", plot = p, width = 8, height = 7, dpi = 300)
