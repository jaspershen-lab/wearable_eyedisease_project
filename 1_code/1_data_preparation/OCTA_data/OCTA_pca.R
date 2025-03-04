######OCTA pca
# 准备OCTA变量
retina_vars <- c(
  # SVD related
  "SVD_NerveFiber_0_6_T0",
  "SVD_Superficial_0_6_T0",

  # PA related
  "PA_Avascular_0_6_T0",
  "PA_Choriocapillaris_0_6_T0",
  "PA_Choroid_0_6_T0",
  "PA_DCP_0_6_T0",
  "PA_Deep_0_6_T0",
  "PA_ICP_0_6_T0",
  "PA_InnerRetina_0_6_T0",
  "PA_NerveFiber_0_6_T0",
  # "PA_OuterRetina_0_6_T0",
  "PA_PED_0_6_T0",
  "PA_Retina_0_6_T0",
  "PA_Superficial_0_6_T0",
  "PA_SVP_0_6_T0",
  "PA_Vitreous_0_6_T0",

  # VD related
  "VD_DCP_0_6_T0",
  "VD_Deep_0_6_T0",
  "VD_ICP_0_6_T0",
  "VD_InnerRetina_0_6_T0",
  "VD_NerveFiber_0_6_T0",
  "VD_Superficial_0_6_T0",
  "VD_SVP_0_6_T0"
)

# 使用 model_data 中的数据进行PCA
retina_data <- model_data[, retina_vars]

# # 标准化数据（保留NA值）
# retina_scaled <- scale(retina_data)

# 检查相关性
cor_matrix <- cor(retina_data)
print("Correlation matrix of retinal features:")
print(round(cor_matrix[1:5, 1:5], 2))

# 最小-最大标准化函数
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# 对每列应用最小-最大标准化
retina_scaled <- as.data.frame(lapply(retina_data, min_max_scale))

# 然后进行插补
imp <- mice(retina_scaled, m=5)
complete_data <- complete(imp, 1)
pca_result <- prcomp(complete_data, scale. = FALSE)  # 因为已经标准化过了
# 使用第一个插补数据集
complete_data <- complete(imp, 1)
# 然后进行PCA
pca_result <- prcomp(complete_data, scale. = TRUE)

# 计算每个主成分解释的方差
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cum_var_explained <- cumsum(var_explained)

# 创建碎石图数据
scree_data <- data.frame(
  PC = 1:length(var_explained),
  Variance = var_explained,
  Cumulative = cum_var_explained
)

# 绘制碎石图
library(ggplot2)
p1 <- ggplot(scree_data, aes(x = PC)) +
  geom_line(aes(y = Variance), color = "blue") +
  geom_point(aes(y = Variance), color = "blue") +
  geom_line(aes(y = Cumulative), color = "red") +
  geom_point(aes(y = Cumulative), color = "red") +
  labs(
    title = "Scree Plot of PCA",
    x = "Principal Component",
    y = "Proportion of Variance Explained"
  ) +
  theme_bw() +
  scale_y_continuous(
    name = "Individual Variance Explained",
    sec.axis = sec_axis(~., name = "Cumulative Variance Explained")
  )

print(p1)

# 确定需要保留的主成分数（解释80%的方差）
n_components <- which(cum_var_explained >= 0.8)[1]
print(paste("\nNumber of components needed to explain 80% of variance:", n_components))

# 获取特征对主成分的贡献
loadings <- pca_result$rotation[, 1:min(5, n_components)]

# 提取主成分得分
pc_scores <- pca_result$x[, 1:n_components]
colnames(pc_scores) <- paste0("PC", 1:n_components)

# 创建载荷热图数据
loadings_df <- as.data.frame(loadings)
loadings_df$Feature <- rownames(loadings_df)
loadings_long <- tidyr::pivot_longer(loadings_df,
                                     cols = -Feature,
                                     names_to = "PC",
                                     values_to = "Loading")

# 创建热图
p2 <- ggplot(loadings_long, aes(x = PC, y = Feature, fill = Loading)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1)) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(title = "Feature Loadings on Principal Components",
       x = "Principal Component",
       y = "Feature")

print(p2)
