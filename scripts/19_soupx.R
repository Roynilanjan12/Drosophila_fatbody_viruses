##Supplementary figure S6

suppressPackageStartupMessages({
  library(harmony)
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(readxl)
  library(SeuratWrappers)
  library(Nebulosa)
  library(scCustomize)
  library(SingleCellExperiment)
  library(SoupX)     
})

set.seed(123)

#----------------------------
# Helper: add replicate Identity from barcode "-1/-2" suffix
#----------------------------
add_identity_from_barcode_suffix <- function(seu, mapping_named_vector) {
  md <- seu@meta.data
  md$barcodes <- rownames(md)
  parts <- data.frame(do.call("rbind", strsplit(as.character(md$barcodes), "-", fixed = TRUE)))
  grp <- as.character(parts$X2)
  grp2 <- ifelse(grp %in% names(mapping_named_vector), mapping_named_vector[grp], grp)
  seu@meta.data$Identity <- grp2
  return(seu)
}

#----------------------------
# Helper: Run SoupX 
#----------------------------
run_soupx_clean <- function(seu) {
  # 1. Run quick clustering (required for autoEstCont)
  message("Running temporary clustering for SoupX...")
  DefaultAssay(seu) <- "RNA"
  seu <- suppressWarnings({
    NormalizeData(seu, verbose = FALSE) |>
      FindVariableFeatures(verbose = FALSE) |>
      ScaleData(verbose = FALSE) |>
      RunPCA(verbose = FALSE) |>
      FindNeighbors(dims = 1:20, verbose = FALSE) |>
      FindClusters(resolution = 0.5, verbose = FALSE)
  })
  
  # 2. Extract Data
  counts_matrix <- GetAssayData(seu, assay = "RNA", layer = "counts")
  
  # Initialize SoupChannel
  # calcSoupProfile = FALSE prevents it from looking for empty droplets automatically
  sc <- SoupChannel(counts_matrix, counts_matrix, calcSoupProfile = FALSE)
  
  # 3.Estimate soup profile from global average
  avg_counts <- rowMeans(counts_matrix)
  
  soup_profile <- data.frame(
    est = avg_counts / sum(avg_counts),  # Normalized to sum to 1
    counts = avg_counts,                 # Raw average counts
    row.names = rownames(counts_matrix)
  )
  
  # Apply the custom profile
  sc <- setSoupProfile(sc, soup_profile)
  
  # 4. Add clusters
  sc <- setClusters(sc, seu$seurat_clusters)
  
  # 5. Estimate contamination (rho)
  # forceAccept=TRUE ensures it runs even if rho is very low/high
  message("Estimating contamination fraction...")
  sc <- autoEstCont(sc, doPlot = FALSE, forceAccept = TRUE)
  
  # 6. Clean the data
  # roundToInt=TRUE ensures integer counts compatible with Seurat
  message("Adjusting counts...")
  out_counts <- adjustCounts(sc, roundToInt = TRUE, verbose = FALSE)
  
  # 7. Add corrected counts as a NEW Assay named "soupX"
  seu[["soupX"]] <- CreateAssayObject(counts = out_counts)
  
  # 8. ROBUST METADATA EXTRACTION
  # Try 'fit' slot first, then fall back to 'metaData' table
  rho_val <- sc$fit$rhoGlobal
  
  if (is.null(rho_val)) {
    # Fallback: Check if rho is stored in the metadata table (common in custom profiles)
    if ("rho" %in% colnames(sc$metaData)) {
      rho_val <- unique(sc$metaData$rho)[1] # Take the first value (usually uniform)
    } else {
      rho_val <- NA
      message("Warning: Could not extract rho value for metadata (but correction ran successfully).")
    }
  }
  
  # Assign to Seurat metadata
  seu$soupX_rho <- rho_val
  
  # Print the result
  if (!is.na(rho_val)) {
    message(paste("Estimated global contamination (rho):", round(rho_val, 4)))
  }
  
  return(seu)
}


#-----------------------------------------------------------
# 1) Read filtered matrices + Create Seurat objects
#-----------------------------------------------------------
fatbody_v.data_v_1 <- Read10X(data.dir = "./data/combined_virgin_bacterial_uninfected_VU")
sdata.1 <- CreateSeuratObject(fatbody_v.data_v_1, project = "Virgin", min.cells = 3, min.features = 200)
sdata.1$type <- "Virgin"
sdata.1 <- JoinLayers(sdata.1)

fatbody_v.data_v_2 <- Read10X(data.dir = "./data/combined_mated_bacterial_uninfected_MU")
sdata.2 <- CreateSeuratObject(fatbody_v.data_v_2, project = "Mated", min.cells = 3, min.features = 200)
sdata.2$type <- "Mated"
sdata.2 <- JoinLayers(sdata.2)

sdata.1 <- add_identity_from_barcode_suffix(sdata.1, c("1" = "V_1", "2" = "V_2"))
sdata.2 <- add_identity_from_barcode_suffix(sdata.2, c("1" = "M_1", "2" = "M_2"))

#-----------------------------------------------------------
# 2) Run SoupX (Adds 'soupX' assay)
#-----------------------------------------------------------
print("Running SoupX on sdata.1...")
sdata.1 <- run_soupx_clean(sdata.1)

print("Running SoupX on sdata.2...")
sdata.2 <- run_soupx_clean(sdata.2)

#-----------------------------------------------------------
# 3) Set Default Assay & Merge
#-----------------------------------------------------------
# IMPORTANT: Switch default to the corrected counts
DefaultAssay(sdata.1) <- "soupX"
DefaultAssay(sdata.2) <- "soupX"

fatbody_v <- merge(sdata.2, y = c(sdata.1), add.cell.ids = c("M", "V"))
rm(fatbody_v.data_v_1, fatbody_v.data_v_2)
gc()

#----------------------------
# 4) Your Seurat QC filtering (as before)
#----------------------------
# Note: QC is usually done on the raw RNA to get accurate Mito %, 
# but we filter the whole object (which removes cells from both RNA and soupX assays)
DefaultAssay(fatbody_v) <- "RNA" 
fatbody_v$mito.percent <- PercentageFeatureSet(fatbody_v, pattern = "^mt:")

# Apply filter
fatbody_v.filtered <- subset(
  fatbody_v,
  subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & mito.percent < 25
)

# SWITCH BACK TO soupX FOR ANALYSIS
DefaultAssay(fatbody_v.filtered) <- "soupX"

#----------------------------
# 5) Normalize / HVG / Scale / PCA / UMAP
#----------------------------
# All these functions will now use the "soupX" assay by default
fatbody_v.filtered <- NormalizeData(fatbody_v.filtered)
fatbody_v.filtered <- FindVariableFeatures(fatbody_v.filtered, selection.method = "vst", nfeatures = 2000)

fatbody_v.filtered <- ScaleData(fatbody_v.filtered) 

fatbody_v.filtered <- RunPCA(
  fatbody_v.filtered,
  npcs = 150,
  verbose = FALSE
)
ElbowPlot(object = fatbody_v.filtered, ndims = 150)

fatbody_v.filtered <- RunUMAP(fatbody_v.filtered, reduction = "pca", dims = 1:50)

before <- DimPlot(fatbody_v.filtered, reduction = "umap", group.by = "type")
print(before)

#----------------------------
# 6) Harmony integration
#----------------------------
fatbody_v.harmony <- fatbody_v.filtered %>%
  RunHarmony(group.by.vars = "Identity", plot_convergence = FALSE)

set.seed(123)
fatbody_v.harmony <- fatbody_v.harmony %>%
  RunUMAP(reduction = "harmony", dims = 1:50) %>%
  FindNeighbors(reduction = "harmony", dims = 1:50) %>%
  FindClusters(resolution = 0.1)

after <- DimPlot(fatbody_v.harmony, reduction = "umap", group.by = "type")
print(after)

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

fatbody.markers <- FindAllMarkers(fatbody_v.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10_59

DimPlot(fatbody_v.harmony, reduction = 'umap', split.by = 'type')

fatbody_v.harmony <- FindClusters(fatbody_v.harmony, resolution = 0.1)

fatbody.markers <- FindAllMarkers(fatbody_v.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10_59

fatbody_subset <- subset(
  fatbody_v.harmony,
  subset = seurat_clusters %in% c("0", "1", "2")
)

#new.cluster.ids <- c("Fatbody cells 1 (dhd,wisp)", "Fatbody cells 2 (yp1,yp3)", "Fatbody cells 3 (Pdfr,Egfr)",
#"Epithelial cells (grh,mgl)","Muscle cells (slo,sals)", "Oenocyte (FASN2,FASN3)", 
#"Hemocytes (Hml,Had2)", "Female reproductive system (Vm26Aa,trol)",
#"Unknown 1 (Oatp58Dc,Smvt)","Sensory neuron (scrt,eag)","Unknown 2 (sns,Cubn)")

new.cluster.ids <- c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3")

names(new.cluster.ids) <- levels(fatbody_subset)
fatbody_subset <- RenameIdents(fatbody_subset, new.cluster.ids)
fatbody_subset[["cluster.ident"]] <- Idents(object = fatbody_subset)
DimPlot(fatbody_subset,reduction = "umap", pt.size = 0.5, label = F, label.size = 5)
###split by replicates.
DimPlot(fatbody_subset,reduction = "umap", pt.size = 0.5, label = F, label.size = 5, split.by = "Identity")


fatbody_v.filtered<- subset(x = fatbody_subset, subset = type == "Mated")

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



####for virgins
fatbody_v.filtered<- subset(x = fatbody_subset, subset = type == "Virgin")

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

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
})

# -----------------------------
# 0) Settings / feature names
# -----------------------------
nv  <- "JX220408.1-Nora-virus-isolate-FR1"
dav <- "KP969946.1-Drosophila-A-virus-isolate-LJ35"

DefaultAssay(fatbody_subset) <- "RNA"

# IMPORTANT (Seurat v5): join layers so GetAssayData works
fatbody_subset <- JoinLayers(fatbody_subset, assay = "RNA")

# -----------------------------
# 1) Tag cells that express BOTH viral genes
# -----------------------------
mat <- GetAssayData(fatbody_subset, assay = "RNA", layer = "data")

both_expr <- (mat[nv, ] > 0) & (mat[dav, ] > 0)
fatbody_subset$Both_Viruses <- both_expr

# -----------------------------
# 2) % cells expressing each SINGLE virus (Percent_Expressing)
# -----------------------------
Idents(fatbody_subset) <- "cluster.ident"

percent_stats <- Percent_Expressing(
  seurat_object = fatbody_subset,
  features      = c(nv, dav),
  threshold     = 0,
  split_by      = "type"
)

percent_stats <- as.data.frame(t(percent_stats))
percent_stats$con_cluster <- rownames(percent_stats)
percent_stats$con_cluster <- gsub("^X", "C", percent_stats$con_cluster)

percent_stats <- melt(percent_stats)

# split “condition_cluster” into two cols
percent_stats <- cbind(
  percent_stats,
  do.call(rbind, strsplit(as.character(percent_stats$con_cluster), "_"))
)
colnames(percent_stats)[4:5] <- c("cluster", "condition")

percent_stats <- percent_stats[, -1]  # drop con_cluster
colnames(percent_stats) <- c("virus", "percent", "cluster", "condition")

# -----------------------------
# 3) % cells co-expressing BOTH viruses
# -----------------------------
both_df <- fatbody_subset@meta.data %>%
  mutate(
    cluster   = fatbody_subset$cluster.ident,
    condition = fatbody_subset$type
  ) %>%
  group_by(cluster, condition) %>%
  summarise(percent = mean(Both_Viruses) * 100, .groups = "drop") %>%
  mutate(virus = "Both")

# -----------------------------
# 4) Combine & tidy
# -----------------------------
df <- bind_rows(percent_stats, both_df)

df$cluster <- gsub("\\.", " ", df$cluster)
df$cluster <- factor(df$cluster, levels = unique(df$cluster))

# Optional: drop DAV and Both from Virgin (your original filter)
df <- df %>%
  filter(!(condition == "Virgin" & virus %in% c(dav, "Both")))

# Set virus order: DAV (red) then NV (blue) then Both (black)
df$virus <- factor(df$virus, levels = c(dav, nv, "Both"))

# Correct color mapping (IMPORTANT!)
virus_cols <- setNames(
  c("red", "blue", "black"),
  c(dav, nv, "Both")
)

# -----------------------------
# 5) BARPLOT
# -----------------------------
ggplot(df, aes(x = cluster, y = percent, fill = virus)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black"
  ) +
  scale_fill_manual(
    name   = "Virus",
    values = virus_cols,
    breaks = c(dav, nv, "Both"),
    labels = c("DAV", "NV", "Both")
  ) +
  facet_wrap(~ condition) +
  labs(
    x = "Cluster",
    y = "% of cells infected with virus"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 20),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )










