####Supplementary figure S16

library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)
library(tidyverse)
library(dplyr)
library(magrittr)
library(Matrix)
library(purrr)
library(reshape2)
library(S4Vectors)
library(tibble)
library(SingleCellExperiment)
library(pheatmap)
library(png)
library(DESeq2)
library(RColorBrewer)
library(Seurat)
library(EnhancedVolcano)
library(Cairo)
library(ggplot2)
#library(Matrix.utils)
library(Nebulosa)


# Load the fatbody_v dataset
fatbody_v.data_v <- Read10X(data.dir = "./data/combined_counts_all_treatments_VU_VI_MU_MI")

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes
fatbody_v <- CreateSeuratObject(counts = fatbody_v.data_v, project = "10X_fatbody_v", assay = "RNA", min.cells = 3, min.features = 200)
#fatbody_v <- CreateSeuratObject(counts = fatbody_v.data_v, project = "10X_fatbody_v", assay = "RNA")
fatbody_v

####adding metadata
metadata<- fatbody_v@meta.data
metadata$barcodes <- rownames(metadata)
barcodes <- data.frame(do.call('rbind', strsplit(as.character(metadata$barcodes),'-',fixed=TRUE)))
barcodes$group <- barcodes$X2
barcodes$group[barcodes$X2 == 1 ] <- "VU_1"
barcodes$group[barcodes$X2 == 2 ] <- "MU_1"
barcodes$group[barcodes$X2 == 3 ] <- "VI_1"
barcodes$group[barcodes$X2 == 4 ] <- "MI_1"
barcodes$group[barcodes$X2 == 5 ] <- "VU_2"
barcodes$group[barcodes$X2 == 6 ] <- "VI_2"
barcodes$group[barcodes$X2 == 7 ] <- "MU_2"
barcodes$group[barcodes$X2 == 8 ] <- "MI_2"
fatbody_v@meta.data$treatment <- barcodes$group

fatbody_v@meta.data <- fatbody_v@meta.data %>%
  mutate(
    # replicate number (1, 2)
    replicate = sub("^.*_", "", treatment),
    
    # VU / VI / MU / MI
    condition = sub("_[0-9]+$", "", treatment),
    
    # V = Virgin, M = Mated
    mating_status = ifelse(
      substr(condition, 1, 1) == "V",
      "Virgin",
      "Mated"
    ),
    
    # U = Uninfected, I = Infected
    bacterial_infection_status = ifelse(
      substr(condition, 2, 2) == "U",
      "Uninfected",
      "Infected"
    )
  )

fatbody_v@meta.data <- fatbody_v@meta.data %>%
  mutate(
    replicate = factor(replicate, levels = c("1", "2")),
    mating_status = factor(mating_status, levels = c("Virgin", "Mated")),
    bacterial_infection_status = factor(
      bacterial_infection_status,
      levels = c("Uninfected", "Infected")
    ),
    condition = factor(condition, levels = c("VU", "VI", "MU", "MI"))
  )


fatbody_v <- JoinLayers(fatbody_v)


# QC and filtering
fatbody_v$mito.percent <- PercentageFeatureSet(fatbody_v, pattern = '^mt:')
#View(fatbody_v@meta.data)
# explore QC


#filter
fatbody_v
fatbody_v.filtered <- subset(fatbody_v, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mito.percent < 25)
fatbody_v.filtered <- NormalizeData(fatbody_v.filtered)
fatbody_v.filtered <- FindVariableFeatures(fatbody_v.filtered, selection.method = "vst", nfeatures = 2000)
top10_v <- head(VariableFeatures(fatbody_v.filtered), 10)
all.genes.v <- rownames(fatbody_v.filtered)

fatbody_v.filtered <- ScaleData(fatbody_v.filtered, features = rownames(fatbody_v.filtered))


fatbody_v.filtered <- RunPCA(fatbody_v.filtered, npcs = 150, verbose = FALSE, features = VariableFeatures(fatbody_v.filtered), nfeatures.print = 10)
ElbowPlot(object = fatbody_v.filtered, ndims = 50)

fatbody_v.filtered <- RunUMAP(fatbody_v.filtered, reduction = "pca", dims = 1:50)

before <- DimPlot(fatbody_v.filtered, reduction = 'umap', group.by = 'mating_status')
before

# run Harmony -----------
fatbody_v.harmony <- fatbody_v.filtered %>%
  RunHarmony(group.by.vars = 'treatment', plot_convergence = FALSE)

fatbody_v.harmony@reductions

fatbody_v.harmony.embed <- Embeddings(fatbody_v.harmony, "harmony")
fatbody_v.harmony.embed[1:10,1:10]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
set.seed(123) 
fatbody_v.harmony <- fatbody_v.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.2)

# visualize 
after <- DimPlot(fatbody_v.harmony, reduction = 'umap', group.by = 'mating_status')
before
after

DimPlot(fatbody_v.harmony, reduction = 'umap', group.by = 'mating_status',  cols = c('Virgin' = 'red', 'Mated' = 'blue'))
DimPlot(fatbody_v.harmony, reduction = 'umap', group.by = 'mating_status')

fatbody.markers <- FindAllMarkers(fatbody_v.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10_59

DimPlot(fatbody_v.harmony, reduction = 'umap', split.by = 'mating_status')

fatbody_v.harmony <- FindClusters(fatbody_v.harmony, resolution = 0.1)

fatbody.markers <- FindAllMarkers(fatbody_v.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10_59

fatbody_v.filtered <- fatbody_v.harmony

#new.cluster.ids <- c("Fatbody cells 1 (dhd,wisp)", "Fatbody cells 2 (yp1,yp3)", "Fatbody cells 3 (Pdfr,Egfr)",
#"Epithelial cells (grh,mgl)","Muscle cells (slo,sals)", "Oenocyte (FASN2,FASN3)", 
#"Hemocytes (Hml,Had2)", "Female reproductive system (Vm26Aa,trol)",
#"Unknown 1 (Oatp58Dc,Smvt)","Sensory neuron (scrt,eag)","Unknown 2 (sns,Cubn)")

fatbody_v.harmony <- subset(
  fatbody_v.harmony,
  subset = seurat_clusters %in% c("0", "1", "2")
)

new.cluster.ids <- c("Fatbody cells 2", "Fatbody cells 1", "Fatbody cells 3")

names(new.cluster.ids) <- levels(fatbody_v.harmony)
fatbody_v.harmony <- RenameIdents(fatbody_v.harmony, new.cluster.ids)
fatbody_v.harmony[["cluster.ident"]] <- Idents(object = fatbody_v.harmony)
DimPlot(fatbody_v.harmony,reduction = "umap", pt.size = 0.5, label = F, label.size = 5)
###split by replicates.
DimPlot(fatbody_v.harmony,reduction = "umap", pt.size = 0.5, label = F, label.size = 5, split.by = "mating_status")


#### DAV and NV comparisons in Mated files:
fatbody_v.harmony1<- subset(x = fatbody_v.harmony, subset = mating_status == "Mated")

coldata <- fatbody_v.harmony1@meta.data
coldata$barcodes <- rownames(coldata)

NV<- FetchData(fatbody_v.harmony1, vars = c('JX220408.1-Nora-virus-isolate-FR1'), slot = "data")
DAV<- FetchData(fatbody_v.harmony1, vars = c('KP969946.1-Drosophila-A-virus-isolate-LJ35'), slot = "data")
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

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggpubr)

###viral load

# 1. Clean Data & Set Order
# Define desired order
cluster_order <- c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3")

plot_data <- df_long %>%
  filter(Expression_Value > 0) %>%
  mutate(
    Expression_Type = recode_factor(Expression_Type,
                                    NV_Expression  = "NV",
                                    DAV_Expression = "DAV"),
    Expression_Type = fct_relevel(Expression_Type, "DAV", "NV"),
    condition = recode(condition, 
                       "MI" = "MBI", 
                       "MU" = "MBU"),
    
    # --- ORDERING APPLIED HERE ---
    cluster.ident = factor(cluster.ident, levels = cluster_order)

  )

# Calculate max y to set bracket positions dynamically
max_y <- max(plot_data$Expression_Value, na.rm = TRUE)

# 2. Plot
ggplot(plot_data, aes(x = condition, y = Expression_Value, col = Expression_Type)) +
  # Jitter and Boxplots with dodging
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             size = 0.5, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0, size = 0.8, 
               position = position_dodge(width = 0.75)) +
  
  #STATISTICAL LAYERS (Red to Red, Blue to Blue)
  
  # 1. Red Bracket (DAV)
  stat_compare_means(data = plot_data %>% filter(Expression_Type == "DAV"),
                     comparisons = list(c("MBU", "MBI")),
                     method = "wilcox.test",
                     label = "p.format",    
                     color = "red",         
                     tip.length = 0.01,
                     label.y = max_y * 1.05, 
                     step.increase = 0) +
  
  # 2. Blue Bracket (NV)
  stat_compare_means(data = plot_data %>% filter(Expression_Type == "NV"),
                     comparisons = list(c("MBU", "MBI")),
                     method = "wilcox.test",
                     label = "p.format",
                     color = "blue",        
                     tip.length = 0.01,
                     label.y = max_y * 1.2,  
                     step.increase = 0) +
  
  # -----------------------------------------------------

labs(y = "Viral normalized reads", 
     x = "Condition", 
     color = "Virus") +
  scale_color_manual(values = c("NV" = "blue", "DAV" = "red")) +
  
  # Facet by Cluster (now follows the specific factor order)
  facet_wrap(~ cluster.ident) + 
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        text = element_text(size = 20),
        strip.background = element_rect(fill = "white"))


####percentage of infected cells
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# 1. Define the Desired Order
cluster_order <- c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3")

# 2. Prepare Data & Force Order
pct_data <- df_long %>%
  # Logic: TRUE if infected (read > 0), FALSE otherwise
  mutate(is_infected = Expression_Value > 0) %>%
  
  # Recode and Factorize
  mutate(
    Expression_Type = recode_factor(Expression_Type,
                                    NV_Expression  = "NV",
                                    DAV_Expression = "DAV"),
    Expression_Type = fct_relevel(Expression_Type, "DAV", "NV"),
    condition = recode(condition, "MI" = "MBI", "MU" = "MBU"),
    
    # --- ORDERING STEP HERE ---
    cluster.ident = factor(cluster.ident, levels = cluster_order)

  ) %>%
  
  # Calculate %
  group_by(cluster.ident, condition, Expression_Type) %>%
  summarise(
    total_cells = n(),
    infected_cells = sum(is_infected),
    pct_infected = (infected_cells / total_cells) * 100,
    .groups = "drop"
  )

# 3. Plot
ggplot(pct_data, aes(x = condition, y = pct_infected, fill = Expression_Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 0.3) +
  
  scale_fill_manual(values = c("NV" = "blue", "DAV" = "red")) +
  labs(y = "% Infected Cells", 
       x = "Condition", 
       fill = "Virus") +
  
  facet_wrap(~ cluster.ident) + 
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        text = element_text(size = 20),
        strip.background = element_rect(fill = "white"))


#### viral load of NV in virgin files:
fatbody_v.harmony1<- subset(x = fatbody_v.harmony, subset = mating_status == "Virgin")

coldata <- fatbody_v.harmony1@meta.data
coldata$barcodes <- rownames(coldata)

NV<- FetchData(fatbody_v.harmony1, vars = c('JX220408.1-Nora-virus-isolate-FR1'), slot = "data")
DAV<- FetchData(fatbody_v.harmony1, vars = c('KP969946.1-Drosophila-A-virus-isolate-LJ35'), slot = "data")
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

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(ggpubr)

# 1. Define Desired Order
cluster_order <- c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3")

# 2. Clean Data & Apply Order
plot_data <- df_long %>%
  filter(Expression_Value > 0) %>%
  mutate(
    Expression_Type = recode_factor(Expression_Type,
                                    NV_Expression  = "NV",
                                    DAV_Expression = "DAV"),
    Expression_Type = fct_relevel(Expression_Type, "DAV", "NV"),
    condition = recode(condition, 
                       "VI" = "VBI", 
                       "VU" = "VBU"),
    
    # --- ORDERING APPLIED HERE ---
    cluster.ident = factor(cluster.ident, levels = cluster_order)

  )

# Calculate max y to set bracket positions dynamically
max_y <- max(plot_data$Expression_Value, na.rm = TRUE)

# 3. Plot
ggplot(plot_data, aes(x = condition, y = Expression_Value, col = Expression_Type)) +
  # Jitter and Boxplots
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             size = 0.5, alpha = 0.4) +
  geom_boxplot(outlier.shape = NA, alpha = 0, size = 0.8, 
               position = position_dodge(width = 0.75)) +
  
  # 1. Red Bracket (DAV)
  stat_compare_means(data = plot_data %>% filter(Expression_Type == "DAV"),
                     comparisons = list(c("VBU", "VBI")),
                     method = "wilcox.test",
                     label = "p.format",    
                     color = "red",         
                     tip.length = 0.01,
                     label.y = max_y * 1.05, 
                     step.increase = 0) +
  
  # 2. Blue Bracket (NV)
  stat_compare_means(data = plot_data %>% filter(Expression_Type == "NV"),
                     comparisons = list(c("VBU", "VBI")),
                     method = "wilcox.test",
                     label = "p.format",
                     color = "blue",        
                     tip.length = 0.01,
                     label.y = max_y * 1.2,  
                     step.increase = 0) +
  

labs(y = "Viral normalized reads", 
     x = "Condition", 
     color = "Virus") +
  scale_color_manual(values = c("NV" = "blue", "DAV" = "red")) +
  
  # Facet by Cluster (will now follow the specific factor order)
  facet_wrap(~ cluster.ident) + 
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        text = element_text(size = 20),
        strip.background = element_rect(fill = "white"))

###percentage of NV infected cells in virgin

library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# 1. Define Desired Order
cluster_order <- c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3")

# 2. Prepare Data (Start with full df_metadata)
pct_data_virgin <- df_metadata %>%
  pivot_longer(cols = c(NV_Expression, DAV_Expression), 
               names_to = "Expression_Type", 
               values_to = "Expression_Value") %>%
  
  # Determine Infection Status (TRUE/FALSE)
  mutate(is_infected = Expression_Value > 0) %>%
  
  # Recode and Apply Factor Order
  mutate(
    Expression_Type = recode_factor(Expression_Type,
                                    NV_Expression  = "NV",
                                    DAV_Expression = "DAV"),
    Expression_Type = fct_relevel(Expression_Type, "DAV", "NV"),
    
    # Recode Virgin Conditions
    condition = recode(condition, 
                       "VI" = "VBI", 
                       "VU" = "VBU"),
    
    # Apply specific cluster order
    cluster.ident = factor(cluster.ident, levels = cluster_order)
  ) %>%
  
  # Filter only for the Virgin conditions we want to plot
  filter(condition %in% c("VBI", "VBU")) %>%
  
  # Calculate Percentage
  group_by(cluster.ident, condition, Expression_Type) %>%
  summarise(
    total_cells = n(),
    infected_cells = sum(is_infected),
    pct_infected = (infected_cells / total_cells) * 100,
    .groups = "drop"
  )

# 3. Plot
ggplot(pct_data_virgin, aes(x = condition, y = pct_infected, fill = Expression_Type)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black", size = 0.3) +
  
  scale_fill_manual(values = c("NV" = "blue", "DAV" = "red")) +
  labs(y = "% Infected Cells", 
       x = "Condition", 
       fill = "Virus") +
  
  facet_wrap(~ cluster.ident) + 
  
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.4),
        text = element_text(size = 20),
        strip.background = element_rect(fill = "white"))

