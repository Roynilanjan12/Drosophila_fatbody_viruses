##### Figure 5b
##### Supplementary figure S26
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first


# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v<- subset(x = fatbody_v.harmony, subset = type == "Mated")
fatbody_v <- subset(x = fatbody_v, idents = c("Fatbody cells 1", "Fatbody cells 2","Fatbody cells 3"))
##### making metadata

coldata <- fatbody_v@meta.data
coldata$barcodes <- rownames(coldata)

# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
NV<- FetchData(fatbody_v, vars = c('JX220408.1-Nora-virus-isolate-FR1'), slot = "data")
DAV<- FetchData(fatbody_v, vars = c('KP969946.1-Drosophila-A-virus-isolate-LJ35'), slot = "data")
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
df_metadata$type<-NULL
df_metadata$seurat_clusters <- NULL
df_metadata$RNA_snn_res.0.1 <- NULL
df_metadata$Both_Viruses <- NULL
df_metadata$cluster_sctype<- NULL
colnames(df_metadata) <- c("barcodes","identity", "clusters", "virus_NV", "virus_DAV")
df_metadata$virus_NV <- as.numeric(df_metadata$virus_NV > 0)
df_metadata$virus_NV <- as.factor(df_metadata$virus_NV)
df_metadata$virus_DAV <- as.numeric(df_metadata$virus_DAV > 0)
df_metadata$virus_DAV <- as.factor(df_metadata$virus_DAV)
rownames(df_metadata) <- df_metadata$barcodes

raw_counts <- fatbody_v@assays$RNA$data
ordering <- colnames(raw_counts)
df_metadata <- df_metadata[match(ordering, df_metadata$barcodes),]

all(colnames(raw_counts) %in% rownames(df_metadata))
all(colnames(raw_counts) == rownames(df_metadata))


count_norm<-as.data.frame(raw_counts+0.01)

assoc_results <- data.frame(Pvalue = numeric())

for (i in 1:nrow(count_norm)) {
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),df_metadata$virus_NV, df_metadata$virus_DAV,df_metadata$identity, df_metadata$clusters)
  colnames(data) <- c("gene","NV","DAV","identity","clusters")
  data$condition <- as.factor(data$NV)
  data$condition <- as.factor(data$DAV)
  data$identity <- as.factor(data$identity)
  data$clusters <- as.factor(data$clusters)
  #try((ancova_model <- aov(gene ~ DAV*clusters + clusters+NV+identity, data = data)),silent = T)

# ---- Model fitting (ANOVA/ANCOVA) ------------------------------------------------------
# The model below tests main effects and interaction terms; Type-III sums of squares
# are used to obtain p-values when multiple factors are present.
# ----------------------------------------------------------------------------
  try((ancova_model <- aov(gene ~ NV*clusters+DAV*clusters+clusters+identity+DAV+NV, data = data)),silent = T)

# ---- Type-III ANOVA table (car::Anova) ------------------------------------------------------
# The model below tests main effects and interaction terms; Type-III sums of squares
# are used to obtain p-values when multiple factors are present.
# ----------------------------------------------------------------------------
  try(d<- Anova(ancova_model, type="III"),silent = T) 
  assoc_results[length(assoc_results$Pvalue)+1, ] = c(d$`Pr(>F)`[7]) ##Pvalue: 6 to get interaction effect genes of NV (cluster effect on infection), 7 for DAV
  print(i)
}


assoc_results$Pvalue <- as.numeric(assoc_results$Pvalue)

padj=p.adjust(assoc_results$Pvalue, method = "fdr")

conditions <- as.data.frame(df_metadata$virus_DAV) ##change it to virus_DAV for DAV or virus_NV for NV
conditions <- t(conditions)
conditions <- as.factor(conditions)

conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))



avona_results <- cbind(foldChanges,assoc_results, padj)
rownames(avona_results)=rownames(count_norm)
avona_results$gene <- row.names(avona_results)
avona_results$foldChanges<- round(avona_results$foldChanges,2)

avona_results$gene[avona_results$gene == 'JX220408.1-Nora-virus-isolate-FR1'] <- 'NV'
#avona_results <- avona_results[!(row.names(avona_results) %in% c("JX220408.1-Nora-virus-isolate-FR1")),]
avona_results <- avona_results[!(row.names(avona_results) %in% c("Mf")),]


# Subset the significant results
sig_res_anova <- dplyr::filter(avona_results, padj < 0.05) %>%
  dplyr::arrange(padj)

#write.csv(sig_res_anova, "DEGsUV_controlling_cluster.csv", quote=F, row.names = T,col.names = T)
avona_results[avona_results == -Inf] <- NA
avona_results <- na.omit(avona_results)

# Load required libraries

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ComplexHeatmap)
library(circlize)
library(viridisLite)

immune_set <- read.csv("/data/List_of_immune_genes_updated.csv", header=T)
immune_set$gene <- immune_set$Symbol
immune_go <- merge(sig_res_anova, immune_set, by = "gene", all.x = TRUE)
immune_go<-NaRV.omit(immune_go)
#immune_go <- immune_go[ order( immune_go$Immune.Process ), ]

ord <- order(immune_go$padj, na.last = NA)
immune_top20 <- immune_go[ ord[seq_len(20)], ]
immune_top20 <- immune_top20[ order( immune_top20$Immune.Process ), ]
# keep only the first term before ";" (and trim spaces)
immune_top20$Immune.Process <- trimws(sub(";.*$", "", immune_top20$Immune.Process))
# 2) recode selected processes to "Cellular"
immune_top20$Immune.Process <- ifelse(
  immune_top20$Immune.Process %in% c("Phagocytosis","Encapsulation","Coagulation","Epithelial"),
  "Cellular",
  immune_top20$Immune.Process
)
genes  <- immune_top20$gene ####important: sig_res_anova$gene to get cell type specific log2FC for all the significant interaction effect genes (Data created: ST5_interaction_effect_genes_in_the_fatbody_cells)

# 1) Define your gene set
#genes <- c("DptA", "DptB", "AttB", "AttC", "Dro", "Drs","Mlc2","up")
#genes <- sig_res_anova$gene
#genes <- genes[-1]
#genes <- genes[ !grepl(":", genes) ]
#genes <- genes[ !grepl("NV", genes) ]
#genes <- immune_go$gene

# 2) Add infection metadata
fatbody_v$DAV <- df_metadata$virus_DAV ##change it to virus_DAV for DAV or virus_NV for NV

# 3) Fetch counts + metadata

# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
data <- FetchData(
  object = fatbody_v,
  vars   = c(genes, "DAV", "cluster.ident"), ## change it to DAV for DAV analysis
  layer  = "data"
)
colnames(data) <- c(genes, "conditions", "clusters")

# 4) Add a pseudocount to avoid zeros in log2
data[genes] <- data[genes] + 0.01

# 5) Pivot to long format
long <- data %>%
  pivot_longer(
    cols      = all_of(genes),
    names_to  = "gene",
    values_to = "count"
  )

# 6) Summarize per (cluster, gene)
summary_df <- long %>%
  group_by(clusters, gene) %>%
  summarise(
    prop_infected = mean(conditions == 1),
    log2FC        = log2(
      mean(count[conditions == 1]) /
        mean(count[conditions == 0])
    ),
    .groups = "drop"
  )

# 7) Build a matrix of log2FC (rows = clusters, cols = genes)
wide_fc <- summary_df %>%
  dplyr::select(clusters, gene, log2FC) %>%
  pivot_wider(names_from = gene, values_from = log2FC)
log2fc_mat <- tibble::column_to_rownames(wide_fc, var = "clusters")
log2fc_mat <- log2fc_mat[, colSums(log2fc_mat != 0) > 0]

# 8) Extract proportion-infected as a named vector
prop_df <- summary_df %>%
  dplyr::select(clusters, prop_infected) %>%
  distinct() %>%
  arrange(clusters)
prop_vec <- prop_df$prop_infected
names(prop_vec) <- prop_df$clusters

# 9) Create a row annotation for proportion infected
ha <- ComplexHeatmap::rowAnnotation(
  "Prop infected" = prop_vec,
  col = list(
    "Prop infected" = circlize::colorRamp2(
      c(min(prop_vec), max(prop_vec)),
      c("grey", "black")
    )
  ),
  width = grid::unit(4, "mm"),
  annotation_legend_param = list(
    title = "Proportion infected"
  )
)

# ——————————————————————————————————————————
# 0) Load libraries
# ——————————————————————————————————————————

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(ComplexHeatmap)
library(circlize)
library(viridisLite)
library(grid)          # for unit()




# ——————————————————————————————————————————
# 1) Build named vector of Immune.Process per gene
# ——————————————————————————————————————————

immune_top20 <- immune_top20 %>%
  filter(gene %in% colnames(log2fc_mat))


gene_proc_vec <- setNames(
  immune_top20$Immune.Process,
  immune_top20$gene
)

# ——————————————————————————————————————————
# 2) Define the order of your processes
# ——————————————————————————————————————————
proc_levels <- unique(gene_proc_vec)

# ——————————————————————————————————————————
# 3) Create a vector of genes ordered by process
# ——————————————————————————————————————————
genes_by_proc <- unlist(lapply(proc_levels, function(p) {
  names(gene_proc_vec)[ gene_proc_vec == p ]
}), use.names = FALSE)

# reorder annotation vector
gene_proc_vec <- gene_proc_vec[ genes_by_proc ]
stopifnot( identical(names(gene_proc_vec), genes_by_proc) )

# ——————————————————————————————————————————
# 4) Define a colour palette for the processes
# ——————————————————————————————————————————
process_cols <- setNames(
  viridisLite::turbo(length(proc_levels)),
  proc_levels
)

# ——————————————————————————————————————————
# 5) Create the top‐of‐heatmap annotation
# ——————————————————————————————————————————
ha_top <- columnAnnotation(
  "Immunological Process" = gene_proc_vec,
  col            = list("Immunological Process" = process_cols),
  annotation_legend_param = list(
    title          = "Immunological Process",
    title_position = "topcenter"
  )
)

# ——————————————————————————————————————————
# 6) Compute square dimensions
# ——————————————————————————————————————————
cell_size     <- unit(5, "mm")                     # <-- each cell 5 mm × 5 mm
total_width   <- length(genes_by_proc) * cell_size # columns × cell_size
total_height  <- nrow(log2fc_mat)     * cell_size  # rows    × cell_size

# ——————————————————————————————————————————
# 7) Draw the heatmap with square cells
# ——————————————————————————————————————————



# ---- Heatmap visualization ------------------------------------------------------
# Heatmaps summarize per-cluster log2 fold-changes for selected gene sets.
# ----------------------------------------------------------------------------
ComplexHeatmap::Heatmap(
  log2fc_mat[, genes_by_proc],
  name               = "Log2 FC",
  col                = circlize::colorRamp2(
    c(min(log2fc_mat), 0, max(log2fc_mat)),
    c("blue", "white","red")
  ),
  cluster_rows       = FALSE,
  cluster_columns    = FALSE,
  column_order       = genes_by_proc,
  top_annotation     = ha_top,
  right_annotation   = ha,    # if you have a `ha` defined elsewhere
  width              = total_width,
  height             = total_height,
  heatmap_legend_param = list(
    title_position = "topcenter"
  )
)


## sig_res_anova (cell type by NV interaction) dataframe was stored in ST5_interaction_effect_genes_in_the_fatbody_cells.xlxs
## from that list of genes GO enrichment was done using flyenricher webtool.
## below is the result:

# Recreate the data from the shown plot
data_plot4 <- data.frame(
  Process = c(
    "nonassociative learning (GO:0046958)",
    "negative regulation of cellular macromolecule biosynthetic process (GO:2000113)",
    "oogenesis (GO:0048477)",
    "nuclear import (GO:0051170)",
    "female gamete generation (GO:0007292)",
    "regulation of mitotic cell cycle (GO:0007346)",
    "positive regulation of cellular process (GO:0048522)",
    "habituation (GO:0046959)",
    "protein import (GO:0017038)",
    "regulation of oskar mRNA translation (GO:0046011)"
  ),
  Score = c(
    52,
    49,
    46,
    42,
    41,
    40,
    39,
    38,
    38,
    37
  )
)

# Reorder factor levels by score
data_plot4$Process <- factor(
  data_plot4$Process,
  levels = data_plot4$Process[order(data_plot4$Score)]
)

# Plot (same style as previous figures)

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(ggplot2)


# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
ggplot(data_plot4, aes(x = Score, y = Process, fill = Score)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    axis.text.x = element_text(angle = 0, hjust = 1)
  ) +
  labs(
    x = "Enrichment Score (-log(P)*Z-score)",
    y = "",
    title = ""
  ) +
  guides(fill = FALSE)



