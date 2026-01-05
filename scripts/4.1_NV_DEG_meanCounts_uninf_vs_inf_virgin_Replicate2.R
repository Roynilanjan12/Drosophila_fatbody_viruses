##Supplementary figure S16
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first

# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v<- subset(x = fatbody_v.harmony, subset = Identity == "V_2")
uninfected <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0)
#DAV <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` > 0)
NV <- subset(x = fatbody_v, subset = `JX220408.1-Nora-virus-isolate-FR1` > 0)


### change every "Mlc2" to other genes like in supplementary figure S16 (alphaTry, betaTry, CG16826, DptA, DptB, AttB, AttC, Dro, Drs)


# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
df_uninfected <- FetchData(object = uninfected, vars = c("Mlc2"), layer = "data")
df_uninfected$condition <- "Uninfected cells"

#df_DAV <- FetchData(object = DAV, vars = c("Mlc2"), layer = "count")
#df_DAV$condition <- "DAV infected cells"

df_NV <- FetchData(object = NV, vars = c("Mlc2"), layer = "data")
df_NV$condition <- "NV infected cells"

df <- rbind(df_uninfected,df_NV)

#df <- df %>% filter(Mlc2 > 0)

# Load the required library

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(ggplot2)

# Assuming your dataframe is called df
# Convert the condition column to a factor for correct Mlc2otting
df$condition <- factor(df$condition)

# Calculate mean and standard error
mean_data <- aggregate(`Mlc2` ~ condition, data = df, FUN = mean)
stderr_data <- aggregate(`Mlc2` ~ condition, data = df, FUN = function(x) sd(x)/sqrt(length(x)))

# Merge mean and standard error data
Mlc2_data <- merge(mean_data, stderr_data, by = "condition")
Mlc2_data<-Mlc2_data %>% arrange(desc(condition))
Mlc2_data$gene <- "Mlc2"

levels(Mlc2_data$condition)
Mlc2_data$condition <- factor(Mlc2_data$condition, levels = c('Uninfected cells', 'NV infected cells'))

# Mlc2otting

# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
ggplot(Mlc2_data, aes(x = gene, y = `Mlc2.x`, color = condition)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = `Mlc2.x` - `Mlc2.y`, ymax = `Mlc2.x` + `Mlc2.y`), width = 0.5, position = position_dodge(width = 0.5)) +
  labs(x = "", y = "Mean counts") +
  theme_classic() + 
  theme(text = element_text(size = 25)) +
  scale_color_manual(values = c("#999999", "blue"))


