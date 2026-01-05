####Figure 2b
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(multcompView)

# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered<- subset(x = fatbody_v.harmony, subset = type == "Mated")

coldata <- fatbody_v.filtered@meta.data
coldata$barcodes <- rownames(coldata)


# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
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


##Reshape the data from wide to long format
df_long <- df_metadata %>%
  pivot_longer(cols = c(NV_Expression, DAV_Expression), 
               names_to = "Expression_Type", 
               values_to = "Expression_Value")


##Reshape the data from wide to long format
df_long <- df_metadata %>%
  pivot_longer(cols = c(NV_Expression, DAV_Expression), 
               names_to = "Expression_Type", 
               values_to = "Expression_Value")

df_long %>%
  filter(Expression_Value > 0) %>%
  mutate(Expression_Type = recode_factor(Expression_Type,
                                         NV_Expression  = "NV",
                                         DAV_Expression = "DAV"),
         # if you still need that particular order:
         Expression_Type = fct_relevel(Expression_Type, "DAV", "NV")) %>%

# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
  ggplot(aes(x = factor(cluster.ident),
             y = Expression_Value,
             col = Expression_Type)) +
  geom_point(position = position_jitterdodge(), size = 0.1, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0,size=0.8) +
  labs(y = "Viral normalized reads", x="Cluster",
       color = "Virus") +
  scale_color_manual(values = c("NV" = "blue", "DAV" = "red")) +
  facet_wrap(~ type, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        text = element_text(size = 20))



# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered<- subset(x = fatbody_v.harmony, subset = type == "Virgin")

coldata <- fatbody_v.filtered@meta.data
coldata$barcodes <- rownames(coldata)


# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
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


##Reshape the data from wide to long format
df_long <- df_metadata %>%
  pivot_longer(cols = c(NV_Expression, DAV_Expression), 
               names_to = "Expression_Type", 
               values_to = "Expression_Value")


##Reshape the data from wide to long format
df_long <- df_metadata %>%
  pivot_longer(cols = c(NV_Expression, DAV_Expression), 
               names_to = "Expression_Type", 
               values_to = "Expression_Value")

df_long %>%
  filter(Expression_Value > 0) %>%
  mutate(Expression_Type = recode_factor(Expression_Type,
                                         NV_Expression  = "NV",
                                         DAV_Expression = "DAV"),
         # if you still need that particular order:
         Expression_Type = fct_relevel(Expression_Type, "DAV", "NV")) %>%

# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
  ggplot(aes(x = factor(cluster.ident),
             y = Expression_Value,
             col = Expression_Type)) +
  geom_point(position = position_jitterdodge(), size = 0.1, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0,size=0.8) +
  labs(y = "Viral normalized reads", x="Cluster",
       color = "Virus") +
  scale_color_manual(values = c("NV" = "blue", "DAV" = "red")) +
  facet_wrap(~ type, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        text = element_text(size = 20))







