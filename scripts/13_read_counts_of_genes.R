##Supplementary figure S20, S21
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first
###read_counts_in_virgin for NV
# ---- Subsetting / filtering ----
# Subset cells by QC thresholds, infection status, or condition.
fatbody_v<- subset(x = fatbody_v.harmony, subset = type == "Virgin")

uninfected <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0)
#DAV <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` > 0)
NV <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` > 0)

##change every "imd" to other genes to get their counts
##IMD genes: PGRP-LC, imd, Tak1, Rel, DptA, AttC
##Toll genes: SPZ, Tl, Myd88, pll, cact, Dif, dl, Drs, Mtk

# ---- Feature extraction ----
# Fetch expression values for specific genes/TEs/viral features from the Seurat object.
df_uninfected <- FetchData(object = uninfected, vars = c("imd"), layer = "data")
df_uninfected$condition <- "Uninfected cells"

#df_DAV <- FetchData(object = DAV, vars = c("imd"), layer = "count")
#df_DAV$condition <- "DAV infected cells"

df_NV <- FetchData(object = NV, vars = c("imd"), layer = "data")
df_NV$condition <- "NV infected cells"

df <- rbind(df_uninfected,df_NV)

#df <- df %>% filter(imd > 0)

# Load the required library
# ---- Load required R packages ----
# Packages are loaded explicitly to make dependencies clear for reproduction.
library(ggplot2)

# Assuming your dataframe is called df
# Convert the condition column to a factor for correct imdotting
df$condition <- factor(df$condition)

# Calculate mean and standard error
mean_data <- aggregate(`imd` ~ condition, data = df, FUN = mean)
stderr_data <- aggregate(`imd` ~ condition, data = df, FUN = function(x) sd(x)/sqrt(length(x)))

# Merge mean and standard error data
imd_data <- merge(mean_data, stderr_data, by = "condition")
imd_data<-imd_data %>% arrange(desc(condition))
imd_data$gene <- "imd"

levels(imd_data$condition)
imd_data$condition <- factor(imd_data$condition, levels = c('Uninfected cells', 'NV infected cells'))

# imdotting
# ---- Visualization ----
# Figures are generated with ggplot2; styling is kept consistent with manuscript themes.
ggplot(imd_data, aes(x = gene, y = `imd.x`, color = condition)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = `imd.x` - `imd.y`, ymax = `imd.x` + `imd.y`), width = 0.5, position = position_dodge(width = 0.5)) +
  labs(x = "", y = "Mean counts") +
  theme_classic() + 
  theme(text = element_text(size = 25)) +
  scale_color_manual(values = c("#999999", "blue"))

#######read counts in mated for both DAV and NV infetion
fatbody_v<- subset(x = fatbody_v.harmony, subset = type == "Mated")

uninfected <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` == 0)
DAV <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` > 0)
NV <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` > 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` == 0)

##change every "Mtk" to other genes to get their counts
##IMD genes: PGRP-LC, imd, Tak1, Rel, DptA, AttC
##Toll genes: SPZ, Tl, Myd88, pll, cact, Dif, dl, Drs, Mtk
df_uninfected <- FetchData(object = uninfected, vars = c("Mtk"), layer = "count")
df_uninfected$condition <- "Uninfected cells"

df_DAV <- FetchData(object = DAV, vars = c("Mtk"), layer = "data")
df_DAV$condition <- "DAV infected cells"

df_NV <- FetchData(object = NV, vars = c("Mtk"), layer = "data")
df_NV$condition <- "NV infected cells"

df <- rbind(df_uninfected,df_DAV,df_NV)

#df <- df %>% filter(dl > 0)

# Load the required library
library(ggplot2)

# Assuming your dataframe is called df
# Convert the condition column to a factor for correct plotting
df$condition <- factor(df$condition)

# Calculate mean and standard error
mean_data <- aggregate(`Mtk` ~ condition, data = df, FUN = mean)
stderr_data <- aggregate(`Mtk` ~ condition, data = df, FUN = function(x) sd(x)/sqrt(length(x)))

# Merge mean and standard error data
plot_data <- merge(mean_data, stderr_data, by = "condition")
plot_data<-plot_data %>% arrange(desc(condition))
plot_data$gene <- "Mtk"

levels(plot_data$condition)
plot_data$condition <- factor(plot_data$condition, levels = c('Uninfected cells', 'DAV infected cells', 'NV infected cells'))

# Plotting
ggplot(plot_data, aes(x = gene, y = `Mtk.x`, color = condition)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = `Mtk.x` - `Mtk.y`, ymax = `Mtk.x` + `Mtk.y`), width = 0.5, position = position_dodge(width = 0.5)) +
  labs(x = "", y = "Mean counts") +
  theme_classic() + 
  theme(text = element_text(size = 25)) +
  scale_color_manual(values = c("#999999", "red", "blue"))






