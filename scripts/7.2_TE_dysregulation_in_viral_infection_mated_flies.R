# ---- Load required R packages ----
# Packages are loaded explicitly to make dependencies clear for reproduction.
library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)
library(readxl)
library(SeuratWrappers)
library(Nebulosa)

fatbody_v.data_v_1 <- Read10X(data.dir = "/data/SoloTE_MU_replicate_1_only")
fatbody_v.data_v_2 <- Read10X(data.dir = "/data/SoloTE_MU_replicate_2_only")

# ---- Construct Seurat objects ----
# Create per-replicate Seurat objects from count matrices; thresholds follow the manuscript QC.
sdata.1 <- CreateSeuratObject(fatbody_v.data_v_1, project = "mu_rep1",min.cells = 3, min.features = 200)
sdata.2 <- CreateSeuratObject(fatbody_v.data_v_2, project = "mu_rep2",min.cells = 3, min.features = 200)

# add metadata
sdata.1$type = "mu_rep1"
sdata.2$type = "mu_rep2"


# Merge datasets into one single seurat object
# ---- Merge replicates into a single object ----
# Merging enables integrated QC, normalization, and batch correction downstream.
fatbody_v <- merge(sdata.1, c(sdata.2), add.cell.ids = c("rep1","rep2"))
fatbody_v <- JoinLayers(fatbody_v)
rm(fatbody_v.data_v_1,fatbody_v.data_v_2,sdata.1,sdata.2)
gc()
as.data.frame(fatbody_v@assays$RNA$counts[1:10, 1:2])
head(fatbody_v@meta.data, 10)


# QC and filtering
# ---- Quality control metrics ----
# Compute mitochondrial fraction and other QC indicators used for filtering.
fatbody_v$mito.percent <- PercentageFeatureSet(fatbody_v, pattern = '^mt:')
#View(fatbody_v@meta.data)
# explore QC

fatbody_v$log10GenesPerUMI <- log10(fatbody_v$nFeature_RNA) / log10(fatbody_v$nCount_RNA)

#filter
fatbody_v
fatbody_v.filtered <- subset(fatbody_v, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mito.percent < 25)
fatbody_v.filtered <- NormalizeData(fatbody_v.filtered)
fatbody_v.filtered <- FindVariableFeatures(fatbody_v.filtered, selection.method = "vst", nfeatures = 2000)
top10_v <- head(VariableFeatures(fatbody_v.filtered), 10)
all.genes.v <- rownames(fatbody_v.filtered)

fatbody_v.filtered <- ScaleData(fatbody_v.filtered, features = rownames(fatbody_v.filtered))

# ---- Dimensionality reduction (PCA) ----
# PCA is used to capture major sources of variation prior to UMAP/clustering.
fatbody_v.filtered <- RunPCA(fatbody_v.filtered, npcs = 150, verbose = FALSE, features = VariableFeatures(fatbody_v.filtered), nfeatures.print = 10)
ElbowPlot(object = fatbody_v.filtered, ndims = 150)

# ---- UMAP embedding and neighborhood graph ----
# UMAP provides a low-dimensional visualization; neighbors/clusters are computed for downstream grouping.
fatbody_v.filtered <- RunUMAP(fatbody_v.filtered, reduction = "pca", dims = 1:50)

before <- DimPlot(fatbody_v.filtered, reduction = 'umap', group.by = 'type')


# run Harmony -----------
fatbody_v.harmony <- fatbody_v.filtered %>%
# ---- Batch correction / integration (Harmony) ----
# Harmony is applied to mitigate replicate/batch effects while preserving biological signal.
  RunHarmony(group.by.vars = 'type', plot_convergence = FALSE)

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
after <- DimPlot(fatbody_v.harmony, reduction = 'umap', group.by = 'type')

before|after

p1<-FeaturePlot(fatbody_v.harmony, features = "JX220408.1-Nora-virus-isolate-FR1")
VlnPlot(fatbody_v.harmony, features = "JX220408.1-Nora-virus-isolate-FR1")
fatbody_v.harmony@reductions$pca <- fatbody_v.harmony@reductions$harmony
fatbody_v.harmony@reductions$harmony <- NULL
p2<-plot_density(fatbody_v.harmony, "JX220408.1-Nora-virus-isolate-FR1")

#before|after|p1|p2
DimPlot(fatbody_v.harmony,reduction = "umap", pt.size = 0.5, label = T, label.size = 4)

# ---- Marker detection (optional) ----
# Marker genes support cell-type annotation and cluster interpretation.
fatbody.markers <- FindAllMarkers(fatbody_v.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10_59


stim_cells = Cells(fatbody_v.harmony)[which(fatbody_v.harmony$type == "mu_rep1")]
ctrl_cells = Cells(fatbody_v.harmony)[which(fatbody_v.harmony$type == "mu_rep2")]
slct_stim_cells = sample(stim_cells, size = 4000)
slct_ctrl_cells = sample(ctrl_cells, size = 4000)
fatbody_v.filtered1 = subset(fatbody_v.harmony, cells = c(slct_stim_cells, slct_ctrl_cells))


# ---- Subsetting / filtering ----
# Subset cells by QC thresholds, infection status, or condition.
uninfected <- subset(x = fatbody_v.filtered1, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` == 0)
DAV <- subset(x = fatbody_v.filtered1, subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` > 0)
NV <- subset(x = fatbody_v.filtered1, subset = `JX220408.1-Nora-virus-isolate-FR1` > 0 & `KP969946.1-Drosophila-A-virus-isolate-LJ35` == 0)

### change "vas" to any transposable element (TE), piRNA and RNAi pathway genes to get the counts
## TE: SoloTE-DNA, SoloTE-RC, SoloTE-LTR, SoloTE-LINE
## piRNA genes: piwi, armi, vret, tj, Hen1, AGO3, aub, squ, vas, lncRNA:flam (flamenco)
## RNAi genes: AGO2, Dcr-2

# ---- Feature extraction ----
# Fetch expression values for specific genes/TEs/viral features from the Seurat object.

df_uninfected <- FetchData(object = uninfected, vars = c("vas"), layer = "data")
df_uninfected$condition <- "Uninfected cells"

df_DAV <- FetchData(object = DAV, vars = c("vas"), layer = "data")
df_DAV$condition <- "DAV infected cells"

df_NV <- FetchData(object = NV, vars = c("vas"), layer = "data")
df_NV$condition <- "NV infected cells"

df <- rbind(df_uninfected,df_DAV,df_NV)

df <- df %>% filter(`vas` > 0)

# Load the required library
library(ggplot2)

# Assuming your dataframe is called df
# Convert the condition column to a factor for correct plotting
df$condition <- factor(df$condition)

# Calculate mean and standard error
mean_data <- aggregate(`vas` ~ condition, data = df, FUN = mean)
stderr_data <- aggregate(`vas` ~ condition, data = df, FUN = function(x) sd(x)/sqrt(length(x)))

# Merge mean and standard error data
plot_data <- merge(mean_data, stderr_data, by = "condition")
plot_data<-plot_data %>% arrange(desc(condition))
plot_data$gene <- "vas"

levels(plot_data$condition)
plot_data$condition <- factor(plot_data$condition, levels = c('Uninfected cells', 'DAV infected cells', 'NV infected cells'))
# ---- Output export ----
# Optional export of intermediate tables/figures for downstream use.
clipr::write_clip(plot_data)
### plot_data will give counts data of TE, piRNA, RNAi pathway genes separately that was recored in xlxs file (count_data_of _TE_piRNA.xlxs) for plotting
### this code will give counts of DAV and NV in mated conditions. 


