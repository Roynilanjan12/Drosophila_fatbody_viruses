###figure 6, 7
### Supplementary figure 27, 28

# ---- Load required R packages ----
# Packages are loaded explicitly to make dependencies clear for reproduction.
library(tidyverse)
## Change data from count_data_of _TE_piRNA.xlxs
## TE: SoloTE-DNA, SoloTE-RC, SoloTE-LTR, SoloTE-LINE
## piRNA genes: piwi, armi, vret, tj, Hen1, AGO3, aub, squ, vas, lncRNA:flam (flamenco)
## RNAi genes : AGO2, Dcr-2


## NV virgin and mated flies (figure 6)
data <- data.frame(
  Condition    = c("NV infected cells", "Uninfected cells", "Uninfected cells", "NV infected cells"),
  vas.x = c(1.771402053,        2.067750122,       1.489138002,       1.699536975),
  vas.y = c(0.028961597,        0.029528311,       0.0250479,         0.05848405),
  gene         = rep("SoloTE-LTR", 4),
  Treatment    = c("Virgin", "Virgin", "Mated", "Mated"),
  stringsAsFactors = FALSE
)

# Re-encode the Condition factor with new display labels
data$Treatment <- factor(data$Treatment, levels = c('Virgin', 'Mated'))
data$Condition <- factor(
  data$Condition,
  levels = c("Uninfected cells", "NV infected cells"),
  labels = c("NV(-)", "NV(+)")
)

# ---- Visualization ----
# Figures are generated with ggplot2; styling is kept consistent with manuscript themes.
ggplot(data, aes(x = Condition, y = vas.x,
                 color = Treatment, group = Treatment)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = vas.x - vas.y,
                    ymax = vas.x + vas.y),
                position = position_dodge(width = 0.5), width = 0.3) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.4) +
  labs(x = "", y = "Mean counts", title = "vas") +
  theme_classic() +theme(text = element_text(size = 20))+
  scale_color_manual(values = c("#999999", "black"))


###### DAV in mated flies (figure 7)
data <- data.frame(
  Condition = c("Uninfected cells", "DAV infected cells"),
  SoloTE_DNA_x = c(0.401011239, 0.492575172),
  SoloTE_DNA_y = c(0.014228689, 0.02373367),
  Gene = c("SoloTE-DNA", "SoloTE-DNA")
)

print(data)

levels(data$Condition)
data$Condition <- factor(data$Condition, levels = c('Uninfected cells', 'DAV infected cells'))

# Plotting
ggplot(data, aes(x = Gene, y = `SoloTE_DNA_x`, color = Condition)) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_errorbar(aes(ymin = `SoloTE_DNA_x` - `SoloTE_DNA_y`, ymax = `SoloTE_DNA_x` + `SoloTE_DNA_y`), width = 0.5, position = position_dodge(width = 0.5)) +
  labs(x = "", y = "Mean counts") +
  theme_classic() + 
  theme(text = element_text(size = 25)) +
  scale_color_manual(values = c("#999999", "red"))


