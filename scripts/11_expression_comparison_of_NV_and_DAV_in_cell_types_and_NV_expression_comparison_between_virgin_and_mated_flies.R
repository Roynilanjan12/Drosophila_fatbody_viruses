## Supplementary figure S8, S9, S10, S11

###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first

#–– 1. Load libraries
# ---- Load required R packages ----
# Packages are loaded explicitly to make dependencies clear for reproduction.
library(Seurat)
library(dplyr)
library(ggplot2)

#–– 2. Pull metadata + viral counts into one data.frame
meta <- fatbody_v.harmony@meta.data
# ---- Feature extraction ----
# Fetch expression values for specific genes/TEs/viral features from the Seurat object.
meta$DAV <- FetchData(fatbody_v.harmony,
                      vars = "KP969946.1-Drosophila-A-virus-isolate-LJ35")[,1]
meta$NV  <- FetchData(fatbody_v.harmony,
                      vars = "JX220408.1-Nora-virus-isolate-FR1")[,1]

meta <- meta %>%
  dplyr::select(Identity, type, cluster.ident, DAV, NV)
meta$Identity <- as.factor(meta$Identity)
meta$type <- as.factor(meta$type)

#–– 3. (Optional) focus only on non-zero NV values
meta_filtered <- meta %>%
  filter(NV > 0)

# to focus on NV values in all the cells uncomment the following
##meta_filtered <- meta

#–– 4. Overall ANOVA: test `type` while controlling for Identity & cluster
# ---- Statistical modeling (ANOVA/ANCOVA) ----
# Fit ANOVA models to compare expression across conditions while controlling for covariates.
aov_res <- aov(NV ~ type + Identity + cluster.ident+DAV,
               data = meta_filtered)
print(summary(aov_res))

#–– 5. Compute one p-value per cluster (controlling for Identity & DAV)
pvals_by_cluster <- meta_filtered %>%
  group_by(cluster.ident) %>%
  summarise(
    p_val = {
      grp <- cur_data()
      fit <- aov(NV ~ type + Identity + DAV, data = grp)
      summary(fit)[[1]]["type", "Pr(>F)"]
    },
    .groups = "drop"
  )

#–– 6. Merge p-values back into the main data
meta_merged <- meta_filtered %>%
  left_join(pvals_by_cluster, by = "cluster.ident")

#–– 7. Prep annotation positions for each cluster
ann_text <- meta_merged %>%
  group_by(cluster.ident) %>%
  summarise(
    x     = 1.5,
    y     = max(NV, na.rm = TRUE) * 1.05,
    label = paste0("p = ", signif(p_val[1], 3))
  )

#–– 8. Plot: faceted violin + boxplot, mean points/text, ANOVA p-values

# 2. Extract the overall F-test p-value for `type`
global_p     <- summary(aov_res)[[1]]["type", "Pr(>F)"]
global_label <- paste0("p = ", signif(global_p, 3))


# ---- Visualization ----
# Figures are generated with ggplot2; styling is kept consistent with manuscript themes.
ggplot(meta_filtered, aes(x = type, y = NV)) +
  geom_violin(trim = FALSE) +
  #geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 1) +
  stat_summary(fun = mean, geom = "text",
               aes(label = round(..y.., 2)), vjust    = -0.5, hjust = -0.2) +
  annotate(
    "text",
    x    = 1.5,
    y    = max(meta_filtered$NV, na.rm = TRUE) * 1.1,
    label = global_label
  ) +
  labs(
    x = NULL,
    y = "Normalized viral reads (mRNA)\nNV"
  ) +
  theme(text = element_text(size = 20))


library(dplyr)
library(ggplot2)
library(emmeans)
library(broom)

# 1. fit the full ANCOVA with an interaction
#focus only on non-zero NV values
fit_full <- aov(
  NV ~ type * cluster.ident + Identity + DAV,
  data = meta_filtered
)

# to focus on NV values in all the cells uncomment the following

#fit_full <- aov(
#NV ~ type * cluster.ident + Identity + DAV,
#data = meta
#)

# 2. extract the simple effect of `type` within each `cluster.ident`
# ---- Estimated marginal means & contrasts ----
# emmeans is used to compute condition contrasts within clusters, with multiple-testing adjustment.
emm <- emmeans(fit_full, ~ type | cluster.ident)

# 3. get pairwise contrasts (Mated vs Virgin) per cluster
pw <- contrast(emm, method = "pairwise", adjust = "tukey")

# 4. turn into a tidy data.frame
contrast_df <- broom::tidy(pw)
# contrast_df will have columns: cluster.ident, contrast, estimate, SE, df, t.ratio, p.value

# 5. (Optional) annotate your faceted plot with those p‐values:
ann_text <- contrast_df %>%
  mutate(
    x     = 1.5,
    y     = meta_filtered %>%
      group_by(cluster.ident) %>%
      summarise(y = max(NV, na.rm=TRUE)) %>%
      pull(y) * 1.05,
    label = paste0("p = ", signif(p.value, 3))
  )

# (1) Define the desired facet order
desired_levels <- c(
  "Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
  "Epithelial cells", "Muscle cells", "Oenocyte",
  "Hemocytes", "Female reproductive system", "Unknown 1",
  "Sensory neuron", "Unknown 2"
)

# (2) Re‐factor both your plotting data and ann_text (if it also has cluster.ident)
meta_filtered$cluster.ident <- factor(meta_filtered$cluster.ident, levels = desired_levels)
ann_text$cluster.ident     <- factor(ann_text$cluster.ident, levels = desired_levels)

# (3) Then run the same ggplot call:
ggplot(meta_filtered, aes(type, NV)) +
  geom_violin(trim = T) +
  # geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 1) +
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = round(..y.., 2)),
    color    = "blue",
    fontface = "bold",
    vjust    = -0.5, hjust=-0.5
  ) +
  geom_text(
    data = ann_text,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~ cluster.ident, scales = "free_y") +
  labs(
    y = "Normalized viral reads (mRNA)\nNV",
    x = NULL
  ) +
  theme(text = element_text(size = 20))




############

fatbody_v.filtered<- subset(x = fatbody_v.harmony, subset = type == "Mated")

coldata <- fatbody_v.filtered@meta.data
coldata$barcodes <- rownames(coldata)

NV<- FetchData(fatbody_v.filtered, vars = c('JX220408.1-Nora-virus-isolate-FR1'), slot = "data")
DAV<- FetchData(fatbody_v.filtered, vars = c('KP969946.1-Drosophila-A-virus-isolate-LJ35'), slot = "data")
colnames(NV) <- c('NV_Expression')
colnames(DAV) <- c('DAV_Expression')
NV$barcodes <- rownames(NV)
DAV$barcodes <- rownames(DAV)
v <- cbind(NV,DAV)
v$barcodes<-NULL
v$barcodes<-NULL
v$barcodes <- rownames(v)
df_metadata <- merge(x=coldata,y=v,by="barcodes",all=TRUE)
df_metadata$orig.ident <- NULL
df_metadata$nCount_RNA <- NULL
df_metadata$nFeature_RNA <- NULL
df_metadata$mito.percent <- NULL
df_metadata$RNA_snn_res.0.2 <- NULL
#df_metadata$type<-NULL
df_metadata$seurat_clusters <- NULL
df_metadata$RNA_snn_res.0.1 <- NULL
df_metadata$NV <- NULL
df_metadata$DAV <- NULL

##Reshape the data from wide to long format, this focuses on all the cells
df_long <- df_metadata %>%
  pivot_longer(cols = c(NV_Expression, DAV_Expression), 
               names_to = "Expression_Type", 
               values_to = "Expression_Value")

### uncomment the following if you want to focus on cells that have DAV and NV expression
#df_long <- df_long %>% filter(Expression_Value > 0)

df_long$Identity <- as.factor(df_long$Identity)
df_long$cluster.ident <- as.factor(df_long$cluster.ident)
df_long$Expression_Type <- as.factor(df_long$Expression_Type)

#–– 4. Overall ANOVA: test `type` while controlling for Identity & cluster
aov_res <- aov(Expression_Value ~ Expression_Type+ Identity + cluster.ident,
               data = df_long)
print(summary(aov_res))

# 2. Extract the overall F-test p-value for `type`
global_p     <- summary(aov_res)[[1]]["Expression_Type", "Pr(>F)"]
global_label <- paste0("p = ", signif(global_p, 3))
#global_label <- paste0("p = < 2e-16")

# First convert to character (so you’re not bound by old factor levels)
df_long$Expression_Type <- as.character(df_long$Expression_Type)

# Replace each full name with the short version
df_long$Expression_Type[df_long$Expression_Type == "NV_Expression"]  <- "NV"
df_long$Expression_Type[df_long$Expression_Type == "DAV_Expression"] <- "DAV"

# (Re-)turn it back into a factor if you need it as a factor
df_long$Expression_Type <- factor(df_long$Expression_Type)


ggplot(df_long, aes(x = Expression_Type, y = Expression_Value)) +
  geom_violin(trim = FALSE) +
  #geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 1) +
  stat_summary(fun = mean, geom = "text",
               aes(label = round(..y.., 2)), vjust    = -0.5, hjust = -0.2) +
  annotate(
    "text",
    x    = 1.5,
    y    = max(meta_filtered$NV, na.rm = TRUE) * 1.1,
    label = global_label
  ) +
  labs(
    x = NULL,
    y = "Normalized viral reads (mRNA)"
  ) +
  theme(text = element_text(size = 20))

#–– libraries
library(dplyr)
library(ggplot2)
library(emmeans)
library(broom)

#–– assume meta_long already exists, with columns:
#    cluster.ident, Identity, virus (factor with levels c("DAV","NV")), expression

#–– 1. Fit the full ANOVA with a virus × cluster.ident interaction
fit_full <- aov(
  Expression_Value ~ Expression_Type * cluster.ident + Identity,
  data = df_long
)

#–– 2. Compute the marginal means and pairwise (Tukey‐adjusted) contrasts of virus within each cluster
emm <- emmeans(fit_full, ~ Expression_Type | cluster.ident)
pw  <- contrast(emm, method = "pairwise", adjust = "tukey")

#–– 3. Tidy into a data.frame
pw_df <- broom::tidy(pw)
# pw_df has columns: cluster.ident, contrast, estimate, SE, df, t.ratio, p.value

#–– 4. Prepare annotation positions
#    get each cluster’s max expression to place the label
pos <- df_long %>%
  group_by(cluster.ident) %>%
  summarise(max_expr = max(Expression_Value, na.rm = TRUE))

ann_text <- pw_df %>%
  left_join(pos, by = "cluster.ident") %>%
  mutate(
    x     = 1.5,
    y     = max_expr * 1.05,
    label = paste0("p = ", signif(p.value, 3))
  )

# (1) Define the desired level order
desired_levels <- c(
  "Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
  "Epithelial cells", "Muscle cells", "Oenocyte",
  "Hemocytes", "Female reproductive system", "Unknown 1",
  "Sensory neuron", "Unknown 2"
)

# (2) Re‐factor your data (and also ann_text if it contains cluster.ident):
df_long$cluster.ident <- factor(df_long$cluster.ident, levels = desired_levels)
ann_text$cluster.ident <- factor(ann_text$cluster.ident, levels = desired_levels)

# (3) Then run the same ggplot call:
ggplot(df_long, aes(x = Expression_Type, y = Expression_Value)) +
  geom_violin(trim = F) +
  # geom_boxplot(width = 0.2, outlier.shape = NA) +
  stat_summary(fun = mean, geom = "point", size = 1) +
  stat_summary(
    fun = mean,
    geom = "text",
    aes(label = round(..y.., 2)),
    color    = "blue",
    fontface = "bold",
    vjust    = -0.5
  ) +
  geom_text(
    data = ann_text,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 4
  ) +
  facet_wrap(~ cluster.ident, scales = "free_y") +
  labs(
    x = NULL,
    y = "Normalized viral reads (mRNA)\nExpression"
  ) +
  theme(text = element_text(size = 20))





