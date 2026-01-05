# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(ggplot2)
library(BiocParallel)
library(readxl)
library(SeuratWrappers)
library(Nebulosa)
library(car)



# ---- Load 10X Genomics count matrix ------------------------------------------------------
# Reads a CellRanger-formatted feature-barcode matrix into R.
# ----------------------------------------------------------------------------
fatbody_v.data_v_1 <- Read10X(data.dir = "/data/SoloTE_VU_replicate_1_only")
fatbody_v.data_v_2 <- Read10X(data.dir = "/data/SoloTE_VU_replicate_2_only")



# ---- Create Seurat object (raw counts) ------------------------------------------------------
# Initializes a Seurat object and applies minimal QC filters (min.cells/min.features).
# ----------------------------------------------------------------------------
sdata.1 <- CreateSeuratObject(fatbody_v.data_v_1, project = "vu_1",min.cells = 3, min.features = 200)
sdata.2 <- CreateSeuratObject(fatbody_v.data_v_2, project = "vu_2",min.cells = 3, min.features = 200)

#sdata.1 <- CreateSeuratObject(fatbody_v.data_v_1, project = "vu_1", min.cells = 3, min.features = 50)
#sdata.2 <- CreateSeuratObject(fatbody_v.data_v_2, project = "vu_2", min.cells = 3, min.features = 50)

# add metadata
sdata.1$type = "NV(-)"
sdata.2$type = "NV(+)"



# Merge datasets into one single seurat object
fatbody_v <- merge(sdata.1, c(sdata.2), add.cell.ids = c("vu_1","vu_2"))
rm(fatbody_v.data_v_1,fatbody_v.data_v_2,sdata.1,sdata.2)
gc()
#as.data.frame(fatbody_v@assays$RNA@counts[1:10, 1:2])
head(fatbody_v@meta.data, 10)


# filter
fatbody_v

# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- subset(fatbody_v, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

# ---- Join assay layers (Seurat v5) ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered<- JoinLayers(object = fatbody_v.filtered)
rm(fatbody_v)

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
fatbody_v.filtered <- RunPCA(fatbody_v.filtered, npcs = 50, verbose = FALSE, features = VariableFeatures(fatbody_v.filtered), nfeatures.print = 10)
ElbowPlot(object = fatbody_v.filtered, ndims = 50)
#set.seed(123) 

# ---- Dimensionality reduction: UMAP ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- RunUMAP(fatbody_v.filtered, reduction = "pca", dims = 1:9)

# ---- Construct cell–cell graph (neighbors) ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- FindNeighbors(fatbody_v.filtered, reduction = "pca", dims = 1:9)

# ---- Graph-based clustering ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- FindClusters(fatbody_v.filtered, resolution = 0.2)
DimPlot(fatbody_v.filtered,reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 5)

DimPlot(fatbody_v.filtered,reduction = "umap", pt.size = 0.5, label = TRUE, label.size = 5, group.by = "type")




# ---- Marker gene discovery per cluster ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered.markers <- FindAllMarkers(fatbody_v.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody_v.filtered.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody_v.filtered.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10_59

new.cluster.ids <- c("Fatbody cells 1 (yp1,yp3)", "Fatbody cells 2 (Hsp27,dhd)", "Fatbody cells 3 (Pdfr,Egfr)",
                     "Epithelial cells (mgh,grh)", "Muscle cells (bt,sals)", "Hemocytes (Hml,Had2)", "Oenocyte (FASN2,CG7910)", "Fatbody cells 4 (Mur18b,CG14292)", "Fatbody cells 5 (Mal-A8,Jon65Aiv)")

#new.cluster.ids <- c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3","Epithelial cells", "Muscle cells", "Hemocytes", "Oenocyte", "Fatbody cells 4", "Fatbody cells 5")

names(new.cluster.ids) <- levels(fatbody_v.filtered)

# ---- Rename clusters to cell-type labels ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- RenameIdents(fatbody_v.filtered, new.cluster.ids)
fatbody_v.filtered[["cluster.ident"]] <- Idents(object = fatbody_v.filtered)
DimPlot(fatbody_v.filtered, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()

stim_cells = Cells(fatbody_v.filtered)[which(fatbody_v.filtered$type == "NV(+)")]
ctrl_cells = Cells(fatbody_v.filtered)[which(fatbody_v.filtered$type == "NV(-)")]
slct_stim_cells = sample(stim_cells, size = 4000)
slct_ctrl_cells = sample(ctrl_cells, size = 4000)

# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered1 = subset(fatbody_v.filtered, cells = c(slct_stim_cells, slct_ctrl_cells))

uninfected <- subset(x = fatbody_v.filtered1, subset = type == "NV(-)")

NV <- subset(x = fatbody_v.filtered1, subset = type == "NV(+)")

### change every "vas" to any transposable element (TE), piRNA and RNAi pathway genes to get the counts
## TE: SoloTE-DNA, SoloTE-RC, SoloTE-LTR, SoloTE-LINE
## piRNA genes: piwi, armi, vret, tj, Hen1, AGO3, aub, squ, vas, lncRNA:flam (flamenco)
## RNAi genes: AGO2, Dcr-2


# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
df_uninfected <- FetchData(object = uninfected, vars = c("vas"), layer = "data")
df_uninfected$condition <- "NV(-)"

#df_DAV <- FetchData(object = DAV, vars = c("vas"), layer = "count")
#df_DAV$condition <- "DAV infected cells"

df_NV <- FetchData(object = NV, vars = c("vas"), layer = "data")
df_NV$condition <- "NV(+)"

df <- rbind(df_uninfected,df_NV)

df <- df %>% filter(`vas` > 0)

# Load the required library

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
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
plot_data$condition <- factor(plot_data$condition, levels = c('NV(-)', 'NV(+)'))
clipr::write_clip(plot_data)
### plot_data will give counts data of TE, piRNA, RNAi pathway genes separately that was recored in xlxs file (count_data_of _TE_piRNA.xlxs) for plotting
### this code will give counts of NV in virgin conditions.


