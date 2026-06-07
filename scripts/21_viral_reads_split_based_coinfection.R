##Supplementary figure S15
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first
library(multcompView)
fatbody_v.filtered<- subset(x = fatbody_v.harmony, subset = Identity == "M_2")

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

# After creating df_metadata and before pivot_longer
df_metadata <- df_metadata %>%
  mutate(
    Infection_Status = case_when(
      NV_Expression > 0 & DAV_Expression > 0 ~ "Co-infected",
      NV_Expression > 0 & DAV_Expression == 0 ~ "NV-only",
      NV_Expression == 0 & DAV_Expression > 0 ~ "DAV-only",
      TRUE ~ "Uninfected"
    )
  )

# Filter for infected cells only and reshape
df_plot <- df_metadata %>%
  filter(Infection_Status != "Uninfected") %>%
  pivot_longer(cols = c(NV_Expression, DAV_Expression), 
               names_to = "Virus", 
               values_to = "Load") %>%
  filter(Load > 0) %>% # Keep only the expression of the virus being measured
  mutate(Virus = recode(Virus, NV_Expression = "NV", DAV_Expression = "DAV"))

# Relevel for consistent plotting
df_plot$Infection_Status <- factor(df_plot$Infection_Status, 
                                   levels = c("NV-only", "DAV-only", "Co-infected"))

# Prepare data: Log-transform the load and filter for the virus of interest
# We only compare NV-only vs Co-infected for this specific test
nv_test_data <- df_plot %>%
  filter(Virus == "NV") %>%
  mutate(log_load = log1p(Load)) %>%
  filter(Infection_Status %in% c("NV-only", "Co-infected"))

# Fit the model: Load ~ Infection_Status + Cluster
# This tests the effect of co-infection while controlling for cell type
fit_nv <- lm(Load ~ Infection_Status + factor(cluster.ident), data = nv_test_data)

# Summary of results
summary(fit_nv)

# Prepare data for DAV
dav_test_data <- df_plot %>%
  filter(Virus == "DAV") %>%
  mutate(log_load = log1p(Load)) %>%
  # Filter to compare DAV alone vs. co-infection
  filter(Infection_Status %in% c("DAV-only", "Co-infected"))

# Fit the model: controlling for cluster (cell type)
fit_dav <- lm(Load ~ Infection_Status + factor(cluster.ident), data = dav_test_data)

# Summary of results
summary(fit_dav)

# Quick extraction of the co-infection effect for both
summary(fit_nv)$coefficients["Infection_StatusCo-infected", ]
summary(fit_dav)$coefficients["Infection_StatusCo-infected", ]

# Extract NV p-value
p_nv <- summary(fit_nv)$coefficients["Infection_StatusCo-infected", "Pr(>|t|)"]

# Extract DAV p-value
p_dav <- summary(fit_dav)$coefficients["Infection_StatusCo-infected", "Pr(>|t|)"]

# Print them to console
print(paste("NV Co-infection P-value:", round(p_nv, 4)))
print(paste("DAV Co-infection P-value:", round(p_dav, 4)))

ggplot(df_plot, aes(x = Infection_Status, y = Load, fill = Virus)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(size = 0.1, alpha = 0.2, position = position_jitter(0.2)) +
  facet_grid(Virus ~ type, scales = "free_y") + # Facet by Virus and Mated/Virgin status
  scale_fill_manual(values = c("NV" = "blue", "DAV" = "red")) +
  labs(
       y = "Viral normalized reads",
       x = "Cell Infection State") +
  theme_bw() +
  theme(text = element_text(size = 16), axis.text.x = element_text(angle = 45, hjust = 1))

df_pvals <- data.frame(
  Virus = c("NV", "DAV"),
  type  = c("Mated", "Mated"),  # adjust if needed
  Infection_Status = "Co-infected",  # THIS is your x
  y = Inf,                      # THIS is your y
  label = c(
    paste0("p = ", signif(p_nv, 3)),
    paste0("p = ", signif(p_dav, 3))
  )
)

ggplot(df_plot, aes(x = Infection_Status, y = Load, fill = Virus)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(size = 0.1, alpha = 0.2, position = position_jitter(0.2)) +
  facet_grid(Virus ~ type, scales = "free_y") +
  scale_fill_manual(values = c("NV" = "blue", "DAV" = "red")) +
  labs(
    y = "Viral normalized reads",
    x = "Cell Infection State"
  ) +
  geom_text(
    data = df_pvals,
    aes(x = Infection_Status, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 1.2,
    size = 5
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


####plot mean in the plot
# 1. Calculate the means for each group before plotting
df_means <- df_plot %>%
  group_by(Virus, type, Infection_Status) %>%
  summarise(mean_load = mean(Load), .groups = "drop")

# 2. Update the Plot
ggplot(df_plot, aes(x = Infection_Status, y = Load, fill = Virus)) +
  # --- Original Boxplot and Jitter ---
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(size = 0.1, alpha = 0.1, position = position_jitter(0.2)) +
  
  # --- NEW: Add Mean Points (White Diamond) ---
  geom_point(data = df_means, 
             aes(y = mean_load), 
             shape = 23,        # Diamond shape
             size = 3,          # Size of the diamond
             fill = "white",    # White fill to stand out against blue/red
             color = "black") + # Black outline
  
  # --- NEW: Add Mean Text Labels ---
  geom_text(data = df_means, 
            aes(y = mean_load, label = round(mean_load, 3)), # Round to 1 decimal
            vjust = -2, hjust = 1.5,
            color = "red", 
            fontface = "bold",
            size = 4) +
  
  # --- Rest of your original styling ---
  facet_grid(Virus ~ type, scales = "free_y") + 
  scale_fill_manual(values = c("NV" = "blue", "DAV" = "red")) +
  labs(
    y = "Viral normalized reads",
    x = "Cell Infection State"
  ) +
  
  # --- Your P-values ---
  geom_text(
    data = df_pvals,
    aes(x = Infection_Status, y = y, label = label),
    inherit.aes = FALSE,
    vjust = 1.2,
    size = 5
  ) +
  
  theme_bw() +
  theme(
    text = element_text(size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )




