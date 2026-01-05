####Figure 1c,d,e and Figure 2a
####Supplementary figure S3, S4, S5, S15, S17


# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)
library(readxl)
library(SeuratWrappers)
library(Nebulosa)


# ---- Load 10X Genomics count matrix ------------------------------------------------------
# Reads a CellRanger-formatted feature-barcode matrix into R.
# ----------------------------------------------------------------------------
fatbody_v.data_v_1 <- Read10X(data.dir = "/data/combined_virgin_bacterial_uninfected_VU")


# ---- Create Seurat object (raw counts) ------------------------------------------------------
# Initializes a Seurat object and applies minimal QC filters (min.cells/min.features).
# ----------------------------------------------------------------------------
sdata.1 <- CreateSeuratObject(fatbody_v.data_v_1, project = "Virgin",min.cells = 3, min.features = 200)

metadata<- sdata.1@meta.data
metadata$barcodes <- rownames(metadata)
barcodes <- data.frame(do.call('rbind', strsplit(as.character(metadata$barcodes),'-',fixed=TRUE)))
barcodes$group <- barcodes$X2
barcodes$group[barcodes$X2 == 1 ] <- "V_1"
barcodes$group[barcodes$X2 == 2 ] <- "V_2"
sdata.1@meta.data$Identity <- barcodes$group
sdata.1$type = "Virgin"

# ---- Join assay layers (Seurat v5) ------------------------------------------------------
# ----------------------------------------------------------------------------
sdata.1 <- JoinLayers(sdata.1)


# ---- Load 10X Genomics count matrix ------------------------------------------------------
# Reads a CellRanger-formatted feature-barcode matrix into R.
# ----------------------------------------------------------------------------
fatbody_v.data_v_2 <- Read10X(data.dir = "/data/combined_mated_bacterial_uninfected_MU")


# ---- Create Seurat object (raw counts) ------------------------------------------------------
# Initializes a Seurat object and applies minimal QC filters (min.cells/min.features).
# ----------------------------------------------------------------------------
sdata.2 <- CreateSeuratObject(fatbody_v.data_v_2, project = "Mated",min.cells = 3, min.features = 200)

metadata<- sdata.2@meta.data
metadata$barcodes <- rownames(metadata)
barcodes <- data.frame(do.call('rbind', strsplit(as.character(metadata$barcodes),'-',fixed=TRUE)))
barcodes$group <- barcodes$X2
barcodes$group[barcodes$X2 == 1 ] <- "M_1"
barcodes$group[barcodes$X2 == 2 ] <- "M_2"
sdata.2@meta.data$Identity <- barcodes$group

sdata.2$type = "Mated"

# ---- Join assay layers (Seurat v5) ------------------------------------------------------
# ----------------------------------------------------------------------------
sdata.2 <- JoinLayers(sdata.2)



# Merge datasets into one single seurat object
fatbody_v <- merge(sdata.2, c(sdata.1), add.cell.ids = c("M","V"))
rm(fatbody_v.data_v_1,fatbody_v.data_v_2,sdata.1,sdata.2)
gc()
head(fatbody_v@meta.data, 10)
fatbody_v <- JoinLayers(fatbody_v)



# QC and filtering

# ---- Compute QC metric: mitochondrial fraction ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v$mito.percent <- PercentageFeatureSet(fatbody_v, pattern = '^mt:')
#View(fatbody_v@meta.data)
# explore QC


#filter
fatbody_v

# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- subset(fatbody_v, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mito.percent < 25)

# ---- Normalize counts (e.g., CPM/relative counts/log-normalize) ------------------------------------------------------
# Normalization choices affect interpretability (e.g., CPM vs log-normalized counts).
# ----------------------------------------------------------------------------
fatbody_v.filtered <- NormalizeData(fatbody_v.filtered)

# ---- Identify highly variable genes ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- FindVariableFeatures(fatbody_v.filtered, selection.method = "vst", nfeatures = 2000)
top10_v <- head(VariableFeatures(fatbody_v.filtered), 10)
all.genes.v <- rownames(fatbody_v.filtered)


# ---- Scale/center expression values (for PCA/UMAP) ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- ScaleData(fatbody_v.filtered, features = rownames(fatbody_v.filtered))



# ---- Dimensionality reduction: PCA ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- RunPCA(fatbody_v.filtered, npcs = 150, verbose = FALSE, features = VariableFeatures(fatbody_v.filtered), nfeatures.print = 10)
ElbowPlot(object = fatbody_v.filtered, ndims = 150)


# ---- Dimensionality reduction: UMAP ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- RunUMAP(fatbody_v.filtered, reduction = "pca", dims = 1:50)

before <- DimPlot(fatbody_v.filtered, reduction = 'umap', group.by = 'type')
before

# run Harmony -----------
fatbody_v.harmony <- fatbody_v.filtered %>%

# ---- Batch correction/integration: Harmony ------------------------------------------------------
# Harmony is used to mitigate batch/replicate effects while preserving biology.
# ----------------------------------------------------------------------------
  RunHarmony(group.by.vars = 'Identity', plot_convergence = FALSE)

fatbody_v.harmony@reductions

fatbody_v.harmony.embed <- Embeddings(fatbody_v.harmony, "harmony")
fatbody_v.harmony.embed[1:10,1:10]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
set.seed(123) 
fatbody_v.harmony <- fatbody_v.harmony %>%

# ---- Dimensionality reduction: UMAP ------------------------------------------------------
# ----------------------------------------------------------------------------
  RunUMAP(reduction = 'harmony', dims = 1:50) %>%

# ---- Construct cell–cell graph (neighbors) ------------------------------------------------------
# ----------------------------------------------------------------------------
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%

# ---- Graph-based clustering ------------------------------------------------------
# ----------------------------------------------------------------------------
  FindClusters(resolution = 0.2)

# visualize 
after <- DimPlot(fatbody_v.harmony, reduction = 'umap', group.by = 'type')
before
after

DimPlot(fatbody_v.harmony, reduction = 'umap', group.by = 'type',  cols = c('Virgin' = 'red', 'Mated' = 'blue'))
DimPlot(fatbody_v.harmony, reduction = 'umap', group.by = 'type')


# ---- Marker gene discovery per cluster ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody.markers <- FindAllMarkers(fatbody_v.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10_59

DimPlot(fatbody_v.harmony, reduction = 'umap', split.by = 'type')


# ---- Graph-based clustering ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.harmony <- FindClusters(fatbody_v.harmony, resolution = 0.1)


# ---- Marker gene discovery per cluster ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody.markers <- FindAllMarkers(fatbody_v.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10_59

#new.cluster.ids <- c("Fatbody cells 1 (dhd,wisp)", "Fatbody cells 2 (yp1,yp3)", "Fatbody cells 3 (Pdfr,Egfr)",
                     #"Epithelial cells (grh,mgl)","Muscle cells (slo,sals)", "Oenocyte (FASN2,FASN3)", 
                     #"Hemocytes (Hml,Had2)", "Female reproductive system (Vm26Aa,trol)",
                     #"Unknown 1 (Oatp58Dc,Smvt)","Sensory neuron (scrt,eag)","Unknown 2 (sns,Cubn)")

new.cluster.ids <- c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
                     "Epithelial cells","Muscle cells", "Oenocyte", 
                     "Hemocytes", "Female reproductive system",
                     "Unknown 1","Sensory neuron","Unknown 2")
                      
names(new.cluster.ids) <- levels(fatbody_v.harmony)

# ---- Rename clusters to cell-type labels ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.harmony <- RenameIdents(fatbody_v.harmony, new.cluster.ids)
fatbody_v.harmony[["cluster.ident"]] <- Idents(object = fatbody_v.harmony)
DimPlot(fatbody_v.harmony,reduction = "umap", pt.size = 0.5, label = F, label.size = 5)
###split by replicates.
DimPlot(fatbody_v.harmony,reduction = "umap", pt.size = 0.5, label = F, label.size = 5, split.by = "Identity")



# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(scCustomize)

FeaturePlot_scCustom(fatbody_v.harmony, features = "JX220408.1-Nora-virus-isolate-FR1", split.by = "type", colors_use = c("lightgrey", "blue"))
FeaturePlot_scCustom(fatbody_v.harmony, features = "KP969946.1-Drosophila-A-virus-isolate-LJ35", split.by = "type", colors_use = c("lightgrey", "blue"))

Idents(fatbody_v.harmony) <- "cluster.ident"


# ---- Compute percent of cells expressing features above threshold ------------------------------------------------------
# ----------------------------------------------------------------------------
percent_stats <- Percent_Expressing(seurat_object = fatbody_v.harmony, features = c("JX220408.1-Nora-virus-isolate-FR1","KP969946.1-Drosophila-A-virus-isolate-LJ35"), threshold = 0, split_by = c("type"))
percent_stats <- as.data.frame(t(percent_stats))
percent_stats$con_cluster<-rownames(percent_stats)
percent_stats$con_cluster <- gsub("^X", "C", percent_stats$con_cluster)

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(reshape2)
percent_stats<-melt(percent_stats)

# Split con_cluster column into two separate columns
percent_stats <- cbind(percent_stats, do.call(rbind, strsplit(as.character(percent_stats$con_cluster), "_")))
colnames(percent_stats)[4:5] <- c("col1", "col2")

# Remove the original con_cluster column
percent_stats <- percent_stats[, -1]
colnames(percent_stats) <- c("virus","percent","cluster","condition")
df <- percent_stats

# Convert cluster to factor for correct ordering on the x-axis
df$cluster <- factor(df$cluster, levels = unique(df$cluster))

df$cluster <- gsub("\\.", " ", df$cluster)

# Plotting

# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
ggplot(df, aes(x = cluster, y = percent, color = virus)) +
  geom_point(size = 3) +
  geom_line(aes(group = virus)) +
  scale_color_manual(values = c("blue", "red")) + # Setting colors to blue and red
  facet_wrap(~ condition, scales = "free") +
  labs(x = "Cluster", y = "% of cells infected with virus", color = "Virus")+
  theme(text = element_text(size = 20),axis.text.x = element_text(angle = 90, hjust = 1))



## ──────────────────────────────────────────────────────────────
## 1.  Tag cells that express *both* viral genes
## ──────────────────────────────────────────────────────────────
mat <- GetAssayData(fatbody_v.harmony, assay = "RNA", slot = "data")

both_expr <- (mat["JX220408.1-Nora-virus-isolate-FR1", ] > 0) &
  (mat["KP969946.1-Drosophila-A-virus-isolate-LJ35", ] > 0)

fatbody_v.harmony$Both_Viruses <- both_expr              # add to metadata


## ──────────────────────────────────────────────────────────────
## 2.  % cells expressing each SINGLE virus  (your original code)
## ──────────────────────────────────────────────────────────────
Idents(fatbody_v.harmony) <- "cluster.ident"


# ---- Compute percent of cells expressing features above threshold ------------------------------------------------------
# ----------------------------------------------------------------------------
percent_stats <- Percent_Expressing(
  seurat_object = fatbody_v.harmony,
  features      = c("JX220408.1-Nora-virus-isolate-FR1",
                    "KP969946.1-Drosophila-A-virus-isolate-LJ35"),
  threshold     = 0,
  split_by      = "type"
)

percent_stats <- as.data.frame(t(percent_stats))
percent_stats$con_cluster <- rownames(percent_stats)
percent_stats$con_cluster <- gsub("^X", "C", percent_stats$con_cluster)


# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(reshape2)
percent_stats <- melt(percent_stats)

# split “condition_cluster” into two cols
percent_stats <- cbind(
  percent_stats,
  do.call(rbind, strsplit(as.character(percent_stats$con_cluster), "_"))
)
colnames(percent_stats)[4:5] <- c("cluster", "condition")
percent_stats <- percent_stats[, -1]          # drop con_cluster
colnames(percent_stats) <- c("virus", "percent", "cluster", "condition")


## ──────────────────────────────────────────────────────────────
## 3.  % cells co-expressing *both* viruses  (new bit)
## ──────────────────────────────────────────────────────────────
library(dplyr)

both_df <- fatbody_v.harmony@meta.data |>
  mutate(cluster   = fatbody_v.harmony$cluster.ident,
         condition = type) |>
  group_by(cluster, condition) |>
  summarise(percent = mean(Both_Viruses) * 100,
            .groups = "drop") |>
  mutate(virus = "Both")

## ──────────────────────────────────────────────────────────────
## 4.  Combine & tidy
## ──────────────────────────────────────────────────────────────
df <- bind_rows(percent_stats, both_df)

df$cluster <- factor(gsub("\\.", " ", df$cluster), levels = unique(gsub("\\.", " ", df$cluster)))

df <- df |>
  dplyr::filter(!(condition == "Virgin" &
                    virus %in% c("KP969946.1-Drosophila-A-virus-isolate-LJ35",
                                 "Both")))


## ──────────────────────────────────────────────────────────────
## 5.  Plot – blue (Nora), red (DAV), black (Both)
## ──────────────────────────────────────────────────────────────

# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
ggplot(df, aes(x = cluster, y = percent, color = virus, group = virus)) +
  geom_point(size = 3) +
  geom_line() +
  scale_color_manual(
    name   = "Virus",                               # legend title
    values = c(
      "JX220408.1-Nora-virus-isolate-FR1"         = "blue",
      "KP969946.1-Drosophila-A-virus-isolate-LJ35" = "red",
      "Both"                                      = "black"
    ),
    breaks = c("KP969946.1-Drosophila-A-virus-isolate-LJ35",
               "JX220408.1-Nora-virus-isolate-FR1",
               "Both"),
    labels = c("DAV", "NV", "Both")                # legend labels
  ) +
  facet_wrap(~ condition) +
  labs(x = "Cluster", y = "% of cells infected with virus") +
  theme_bw()+
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

###barplot of percentage of infected cells
df$virus <- factor(
  df$virus,
  levels = c(
    "KP969946.1-Drosophila-A-virus-isolate-LJ35",  # DAV → red
    "JX220408.1-Nora-virus-isolate-FR1",           # NV  → blue
    "Both"                                        # black
  )
)

ggplot(df, aes(x = cluster, y = percent, fill = virus)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(
    name = "Virus",
    values = c(
      "KP969946.1-Drosophila-A-virus-isolate-LJ35" = "red",
      "JX220408.1-Nora-virus-isolate-FR1"          = "blue",
      "Both"                                      = "black"
    ),
    labels = c("DAV", "NV", "Both")
  ) +
  scale_y_continuous(limits = c(0, 60)) +
  facet_wrap(~ condition) +
  labs(x = "Cluster", y = "% of cells infected \n with virus") +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
  )



####Nebulosa
fatbody_v.harmony@reductions$pca <- fatbody_v.harmony@reductions$harmony
fatbody_v.harmony@reductions$harmony <- NULL
plot_density(fatbody_v.harmony, features = "KP969946.1-Drosophila-A-virus-isolate-LJ35")+facet_grid(.~fatbody_v.harmony$type=="Mated")


# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
p <- SCpubr::do_BoxPlot(sample = subset(x = fatbody_v.harmony, subset = `JX220408.1-Nora-virus-isolate-FR1` > 0),
                        feature = c("JX220408.1-Nora-virus-isolate-FR1"),use_test = TRUE,comparisons = list(c("Epithelial cells", "Muscle cells")))

p


categorize_infection <- function(expression) {
  if (expression >= 0 & expression <= 0.5) {
    return("NE")
  } else if (expression > 0.5 & expression <= 3) {
    return("0.5-3")
  } else if (expression > 3) {
    return(">3")
  }
}

# Universal function for plotting gene expression categories remains unchanged,
# but now it uses the updated categorization logic.
plot_gene_expression <- function(seurat_obj, gene_name, reduction = 'umap') {
  # Fetch the expression data for the specified gene

# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
  gene_expr <- FetchData(seurat_obj, vars = gene_name)
  
  # Apply the categorization function to the gene expression data
  seurat_obj$virus_infection <- apply(gene_expr, 1, categorize_infection)
  
  # Define a custom color palette
  custom_colors <- c(
    "NE" = "#D8D8D8",
    "0.5-3" = "black",
    ">3" = "maroon1"
  )
  
  # Reorder the 'virus_infection' factor levels
  seurat_obj$virus_infection <- factor(
    seurat_obj$virus_infection,
    levels = c(">3", "0.5-3", "NE") # Explicit order
  )
  
  # Generate the DimPlot
  DimPlot(
    seurat_obj, 
    reduction = reduction, 
    group.by = 'virus_infection',
    split.by = "type",  ####change it to "Identity" if want to split by all replicates
    cols = custom_colors
  ) + ggtitle(gene_name)
}

plot_gene_expression(fatbody_v.harmony, "JX220408.1-Nora-virus-isolate-FR1")
plot_gene_expression(fatbody_v.harmony, "KP969946.1-Drosophila-A-virus-isolate-LJ35")

##looking at tropism at different expression threshold cutoff.
FeaturePlot(fatbody_v.harmony, "JX220408.1-Nora-virus-isolate-FR1" , split.by = "type", min.cutoff= 0.5)
FeaturePlot(fatbody_v.harmony, "JX220408.1-Nora-virus-isolate-FR1" , split.by = "type", min.cutoff= 1)
FeaturePlot(fatbody_v.harmony, "JX220408.1-Nora-virus-isolate-FR1" , split.by = "type", min.cutoff= 1.5)

FeaturePlot(fatbody_v.harmony, "KP969946.1-Drosophila-A-virus-isolate-LJ35" , split.by = "type", min.cutoff= 0.5)
FeaturePlot(fatbody_v.harmony, "KP969946.1-Drosophila-A-virus-isolate-LJ35" , split.by = "type", min.cutoff= 1)
FeaturePlot(fatbody_v.harmony, "KP969946.1-Drosophila-A-virus-isolate-LJ35" , split.by = "type", min.cutoff= 1.5)


#### UMAP-based comparison of differentially expressed genes between uninfected and Nora virus infected flies (virgin) in the fat body.
####UMAP plot for VU

# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered1 = subset(fatbody_v.harmony, subset = type == "Virgin")

FeaturePlot_scCustom(fatbody_v.filtered1, features = "JX220408.1-Nora-virus-isolate-FR1", split.by = "Identity", colors_use = c("lightgrey", "blue"))
FeaturePlot_scCustom(fatbody_v.filtered1, features = "DptA", split.by = "Identity", colors_use = c("lightgrey", "blue"))
FeaturePlot_scCustom(fatbody_v.filtered1, features = "AttC", split.by = "Identity", colors_use = c("lightgrey", "blue"))
FeaturePlot_scCustom(fatbody_v.filtered1, features = "Drs", split.by = "Identity", colors_use = c("lightgrey", "blue"))
FeaturePlot_scCustom(fatbody_v.filtered1, features = "Mtk", split.by = "Identity", colors_use = c("lightgrey", "blue"))
FeaturePlot_scCustom(fatbody_v.filtered1, features = "Dro", split.by = "Identity", colors_use = c("lightgrey", "blue"))

####UMAP plot for MU

fatbody_v.filtered1 = subset(fatbody_v.harmony, subset = type == "Mated")
uninfected <- subset(x = fatbody_v.filtered1, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` == 0)
DAV <- subset(x = fatbody_v.filtered1, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` > 0)
NV <- subset(x = fatbody_v.filtered1, subset = `JX220408.1-Nora-virus-isolate-FR1` > 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` == 0)

FeaturePlot(uninfected, features = c("CG33926")) | FeaturePlot(DAV, features = c("CG33926"))
FeaturePlot(uninfected, features = c("armi")) | FeaturePlot(DAV, features = c("armi"))
FeaturePlot(uninfected, features = c("aub")) | FeaturePlot(DAV, features = c("aub"))
FeaturePlot(uninfected, features = c("Vago")) | FeaturePlot(DAV, features = c("Vago"))


