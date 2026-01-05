####Figure 1a
####Supplementary figure S1


# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
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

# ---- Load 10X Genomics count matrix ------------------------------------------------------
# Reads a CellRanger-formatted feature-barcode matrix into R.
# ----------------------------------------------------------------------------
fatbody_v.data_v <- Read10X(data.dir = "/data/combined_counts_all_treatments_VU_VI_MU_MI")

# Initialize the Seurat object with the raw (non-normalized data).  Keep all
# genes expressed in >= 3 cells (~0.1% of the data). Keep all cells with at
# least 200 detected genes

# ---- Create Seurat object (raw counts) ------------------------------------------------------
# Initializes a Seurat object and applies minimal QC filters (min.cells/min.features).
# ----------------------------------------------------------------------------
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
fatbody_v@meta.data$Treatment <- barcodes$group

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats

# ---- Compute QC metric: mitochondrial fraction ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v[["percent.mt"]] <- PercentageFeatureSet(fatbody_v, pattern = "^mt:")
# Show QC metrics for the first 5 cells
head(fatbody_v@meta.data, 5)

# Visualize QC metrics as a violin plot
#VlnPlot(fatbody_v, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

#plot1 <- FeatureScatter(fatbody_v, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(fatbody_v, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#plot1 + plot2


# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v <- subset(fatbody_v, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
#fatbody_v <- subset(fatbody_v, subset = percent.mt < 25)
head(fatbody_v@meta.data, 5)



# ---- Normalize counts (e.g., CPM/relative counts/log-normalize) ------------------------------------------------------
# Normalization choices affect interpretability (e.g., CPM vs log-normalized counts).
# ----------------------------------------------------------------------------
seurat_obj <- NormalizeData(
  object             = fatbody_v,
  assay              = "RNA",
  normalization.method = "RC",      # "RC" = Relative Counts
  scale.factor       = 1e6          # scale each cell so that total counts = 1,000,000
)

# 1. Get normalized data (CPM, since method = "RC" with scale.factor = 1e6)
cpm_mat <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

# Define genes and desired column names
genes_of_interest <- c("JX220408.1-Nora-virus-isolate-FR1", "KP969946.1-Drosophila-A-virus-isolate-LJ35", "DCV-gp1", "Newfield-gp1")
gene_labels <- c(
  "JX220408.1-Nora-virus-isolate-FR1" = "Nora virus",
  "KP969946.1-Drosophila-A-virus-isolate-LJ35" = "Drosophila A virus",
  "DCV-gp1" = "Drosophila C virus",
  "Newfield-gp1" = "Newfield virus"
)

# Get CPM matrix
cpm_mat <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")

# Check which genes are present
genes_present <- intersect(genes_of_interest, rownames(cpm_mat))
if (length(genes_present) < length(genes_of_interest)) {
  warning("Missing genes: ", paste(setdiff(genes_of_interest, genes_present), collapse = ", "))
}

# Fetch data: CPM values + metadata

# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
df_cpm <- FetchData(
  object = seurat_obj,
  vars   = c(genes_present, "Treatment"),
  slot   = "data"
)


# Calculate mean only for those present genes
mean_cpm <- df_cpm %>%
  group_by(Treatment) %>%
  summarise(across(all_of(genes_present), mean), .groups = "drop")

mean_cpm$Tissue <- "Fatbody"
mean_cpm <- mean_cpm[, c("Tissue", setdiff(names(mean_cpm), "Tissue"))]

colnames(mean_cpm) <- c("Tissue", "Treatment", "NV", "DAV")


# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(ggbreak)

# mean_cpm already in memory:
# Tissue | Treatment | NV | DAV

plot_df <- mean_cpm %>%
  dplyr::select(Treatment, NV, DAV) %>%
  pivot_longer(cols = c(NV, DAV), names_to = "Virus", values_to = "CPM") %>%
  separate(Treatment, into = c("Condition", "Replicate"), sep = "_", remove = FALSE) %>%
  mutate(
    Condition = factor(Condition, levels = c("VU","VI","MU","MI")),
    Replicate = factor(Replicate, levels = c("1","2"),
                       labels = c("Replicate 1","Replicate 2")),
    Virus = factor(Virus, levels = c("DAV","NV"))  # DAV first, then NV
  ) %>%
  complete(Condition, Replicate, Virus, fill = list(CPM = 0)) %>%
  arrange(Condition)

# Mean CPM (white tick) per panel
mean_df <- plot_df %>%
  group_by(Replicate, Condition) %>%
  summarize(mean_cpm = mean(CPM, na.rm = TRUE), .groups = "drop")



# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
p <- ggplot(plot_df, aes(x = Condition, y = CPM, fill = Virus)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = ifelse(CPM < 5, "", round(CPM, 1))),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3.5
  ) +
  facet_wrap(~ Replicate, nrow = 1) +
  scale_fill_manual(values = c(DAV = "red", NV = "blue")) +
  labs(
    x = "Conditions",
    y = "Mean Viral CPM",
    fill = "Virus"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    strip.background = element_rect(color = "black", fill = NA),
    strip.text = element_text(face = "bold")
  )

print(p)


# keep only VU and MU
plot_df_vm <- plot_df %>%
  filter(Condition %in% c("VU","MU")) %>%
  mutate(Condition = factor(Condition, levels = c("VU","MU"))) %>%
  droplevels()


p <- ggplot(plot_df_vm, aes(x = Condition, y = CPM, fill = Virus)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = ifelse(CPM == 0, "", round(CPM, 1))),
    position = position_dodge(width = 0.8),
    vjust = -0.6,   # push labels a bit higher above bars
    size  = 5       # bigger text
  ) +
  facet_wrap(~ Replicate, nrow = 1) +
  scale_fill_manual(values = c(DAV = "red", NV = "blue")) +
  labs(
    x = "Conditions",
    y = "Mean Viral CPM",
    fill = "Virus"
  ) +
  # give some extra space on top so labels are not cut
  expand_limits(y = max(plot_df_vm$CPM, na.rm = TRUE) * 1.15) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "right",
    strip.background = element_rect(color = "black", fill = NA),
    strip.text = element_text(face = "bold")
  )

print(p)
