##Figure 4a
##Supplementary figure S13, S14


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
library(monocle3)


# ---- Load 10X Genomics count matrix ------------------------------------------------------
# Reads a CellRanger-formatted feature-barcode matrix into R.
# ----------------------------------------------------------------------------
fatbody_v.data_v_1 <- Read10X(data.dir = "/data/VU_replicate_1_only/filtered_feature_bc_matrix")
fatbody_v.data_v_2 <- Read10X(data.dir = "/data/VU_replicate_2_only/filtered_feature_bc_matrix")


# ---- Create Seurat object (raw counts) ------------------------------------------------------
# Initializes a Seurat object and applies minimal QC filters (min.cells/min.features).
# ----------------------------------------------------------------------------
sdata.1 <- CreateSeuratObject(fatbody_v.data_v_1, project = "vu_rep1",min.cells = 3, min.features = 200)
sdata.2 <- CreateSeuratObject(fatbody_v.data_v_2, project = "vu_rep2",min.cells = 3, min.features = 200)

# add metadata
sdata.1$type = "NV(-)"
sdata.2$type = "NV(+)"


# Merge datasets into one single seurat object
fatbody_v <- merge(sdata.1, c(sdata.2), add.cell.ids = c("rep1","rep2"))

# ---- Join assay layers (Seurat v5) ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v <- JoinLayers(fatbody_v)
rm(fatbody_v.data_v_1,fatbody_v.data_v_2,sdata.1,sdata.2)
gc()
as.data.frame(fatbody_v@assays$RNA$counts[1:10, 1:2])
head(fatbody_v@meta.data, 10)


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

fatbody_v.filtered$log10GenesPerUMI <- log10(fatbody_v$nFeature_RNA) / log10(fatbody_v$nCount_RNA)



# ---- Dimensionality reduction: PCA ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- RunPCA(fatbody_v.filtered, npcs = 150, verbose = FALSE, features = VariableFeatures(fatbody_v.filtered), nfeatures.print = 10)
ElbowPlot(object = fatbody_v.filtered, ndims = 150)


# ---- Dimensionality reduction: UMAP ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.filtered <- RunUMAP(fatbody_v.filtered, reduction = "pca", dims = 1:50)

before <- DimPlot(fatbody_v.filtered, reduction = 'umap', group.by = 'type')


# run Harmony -----------
fatbody_v.harmony <- fatbody_v.filtered %>%

# ---- Batch correction/integration: Harmony ------------------------------------------------------
# Harmony is used to mitigate batch/replicate effects while preserving biology.
# ----------------------------------------------------------------------------
  RunHarmony(group.by.vars = 'type', plot_convergence = FALSE)

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

before|after

p1<-FeaturePlot(fatbody_v.harmony, features = "JX220408.1-Nora-virus-isolate-FR1")
VlnPlot(fatbody_v.harmony, features = "JX220408.1-Nora-virus-isolate-FR1")
fatbody_v.harmony@reductions$pca <- fatbody_v.harmony@reductions$harmony
fatbody_v.harmony@reductions$harmony <- NULL
p2<-plot_density(fatbody_v.harmony, "JX220408.1-Nora-virus-isolate-FR1")

#before|after|p1|p2
DimPlot(fatbody_v.harmony,reduction = "umap", pt.size = 0.5, label = T, label.size = 4)


# ---- Marker gene discovery per cluster ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody.markers <- FindAllMarkers(fatbody_v.harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10_59


new.cluster.ids <- c("Fatbody cells 1 (dhd,wisp)", "Fatbody cells 2 (yp1,yp3)", "Fatbody cells 3 (Pdfr,Egfr)",
                     "Muscle cells (bt,sals)", "Epithelial cells (grh,mgl)", "Hemocytes (Hml,Had2)", 
                     "Oenocyte (FASN2,FASN3)", "Female reproductive system (Fas2,Six4)")
names(new.cluster.ids) <- levels(fatbody_v.harmony)

# ---- Rename clusters to cell-type labels ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v.harmony <- RenameIdents(fatbody_v.harmony, new.cluster.ids)
fatbody_v.harmony[["cluster.ident"]] <- Idents(object = fatbody_v.harmony)
DimPlot(fatbody_v.harmony,reduction = "umap", pt.size = 0.5, label = F, label.size = 5)

fatbody_v.harmony[["cluster.ident"]] <- Idents(object = fatbody_v.harmony)

df <- as.data.frame(melt(table(fatbody_v.harmony@meta.data$cluster.ident, fatbody_v.harmony@meta.data$type)))
colnames(df) <- c("Cell Types", "Sample", "Ncells")


# ---- Compute percent of cells expressing features above threshold ------------------------------------------------------
# ----------------------------------------------------------------------------
percent_stats <- Percent_Expressing(seurat_object = fatbody_v.harmony,features = c("JX220408.1-Nora-virus-isolate-FR1"), threshold = 0, group_by = "cluster.ident")
percent_stats <- melt(percent_stats)
colnames(percent_stats) <- c("Cell Types", "% of NV Expression")

count_norm <- as.data.frame(fatbody_v.harmony@assays$RNA$data)

conditions <- as.data.frame(fatbody_v.harmony$type)
colnames(conditions) <- c("type")
conditions <- as.factor(conditions$type)

clusters <- as.data.frame(fatbody_v.harmony$seurat_clusters)
colnames(clusters) <- c("clusters")
clusters <- as.factor(clusters$clusters)


assoc_results <- data.frame(Pvalue = numeric())

for (i in 1:nrow(count_norm)) {
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions, clusters)
  #data$gene <- data$gene+0.0001

# ---- Model fitting (ANOVA/ANCOVA) ------------------------------------------------------
# The model below tests main effects and interaction terms; Type-III sums of squares
# are used to obtain p-values when multiple factors are present.
# ----------------------------------------------------------------------------
  try((ancova_model <- aov(gene ~ conditions+conditions*clusters + clusters, data = data)),silent = T)

# ---- Type-III ANOVA table (car::Anova) ------------------------------------------------------
# The model below tests main effects and interaction terms; Type-III sums of squares
# are used to obtain p-values when multiple factors are present.
# ----------------------------------------------------------------------------
  try(d<- Anova(ancova_model, type="III"),silent = T) 
  assoc_results[length(assoc_results$Pvalue)+1, ] = c(d$`Pr(>F)`[2]) ## get p-value: 2 for condition (NV-0/1), 3 for cluster effect, 4 for cluster by NV intereaction effect
  print(i)
}

assoc_results$Pvalue <- as.numeric(assoc_results$Pvalue)
padj=p.adjust(assoc_results$Pvalue, method = "fdr")

conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
## foldchange here is only for NV uninfected vs infected condition. 
## conditions*clusters term here gives us which genes have interaction effect (cluster by NV intereaction effect).
## later foldchange of interaction effect genes by different cell types (cluster) are calculated using script ######

ancova_results <- cbind(foldChanges,assoc_results, padj)
rownames(ancova_results)=rownames(count_norm)
ancova_results$gene <- row.names(ancova_results)

# Subset the significant results
significant_results <- dplyr::filter(ancova_results, padj < 0.05, abs(foldChanges)>0.5)


# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(IDPmisc)
ancova_results <- NaRV.omit(ancova_results)
##ancova_results is in DEGs_VU1_vs_VU2.xlsx

library(readxl)
library(EnhancedVolcano)
library(ggpubr)
library(IDPmisc)

volcano <- read_excel("/data/DEGs_VU1_vs_VU2.xlsx")
volcano$gene <- gsub("JX220408\\.1-Nora-virus-isolate-FR1", "NV", volcano$gene)
volcano <- as.data.frame(volcano)
volcano$foldChanges <- as.numeric(volcano$foldChanges)
go <- filter(volcano,padj<1e-05, abs(foldChanges)>0.5)

immune_set <- read.csv("/data/List_of_immune_genes_updated.csv", header=T)
immune_set$Gene <- immune_set$Symbol
immune_go <- merge(go, immune_set, by = "Gene", all.x = TRUE)
immune_go<-NaRV.omit(immune_go)
ggbarplot(immune_go, x = "Gene", y = "foldChanges",
          fill = "Immune.Process",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE, 
          rotate = TRUE,
          ggtheme = theme_minimal()
          
)+theme(text = element_text(size = 15))+ylab("Log2FC")


# Create a vector of colors based on conditions
keyvals <- ifelse(volcano$foldChanges < -0.5 & volcano$Pvalue < 5.298294e-06, 'royalblue',
                  ifelse(volcano$foldChanges > 0.5 & volcano$Pvalue < 5.298294e-06, 'red2',
                         'grey30'))

# Replace NA values with 'grey30'
keyvals[is.na(keyvals)] <- 'grey30'

# Assign names to the vector elements based on their values
names(keyvals) <- ifelse(keyvals == 'red2', 'Upregulated',
                         ifelse(keyvals == 'grey30', 'Not Significant',
                                ifelse(keyvals == 'royalblue', 'Downregulated', NA)))





# ---- Volcano plot for differential expression results ------------------------------------------------------
# ----------------------------------------------------------------------------
EnhancedVolcano(volcano,
                x = 'foldChanges',
                y = 'padj',
                lab = volcano$gene,
                selectLab = c('DptA','Drs','betaTry',"alphaTry","Dro","DptB", "Mtk",
                              "AttB","AttC",
                              "Mlc2", "CG16826",
                              "Mhc","Tm2","up",
                              "Cyp6w1","Cyp12d1-p",
                              "dnr1","NV"),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = '',
                subtitle = "",
                pCutoff = 5.298294e-06,
                FCcutoff = 0.5,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                drawConnectors = TRUE,
                pointSize = 2.0,
                labSize = 4.0,
                widthConnectors = 0.5,
                legendPosition = 'right',
                colCustom = keyvals,
                labCol = 'black',
                boxedLabels = TRUE)




# Load necessary library

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(ggplot2)

## base on pvalue <1e-05 (NV in virgin flies) from DEGs_VU1_vs_VU2.xlsx/ST_6_DEGs_VU1_vs_VU2.xlsx data GO enrichment analysis was
## done using flyenricher webtool. The result from the flyenricher webtool is below for plotting

# Create a dataframe with the provided data
data <- data.frame(
  Process = c("Antibacterial Humoral Response", "Defense Response to Gram-positive Bacterium", 
              "Muscle Contraction", "Actomyosin Structure Organization", 
              "Striated Muscle Cell Development", "Oligosaccharide Catabolic Process", 
              "Defense Response to Bacterium", "Proteolysis Involved in Cellular Protein Catabolic Process", 
              "Proteasomal Ubiquitin-independent Protein Catabolic Process", "Polyol Biosynthetic Process"),
  Score = c(69.09, 57.96, 34.39, 33.65, 31.51, 29.65, 28.57, 27.27, 26.84, 24.96)
)

# Reorder the levels of the Process variable based on the Score
data$Process <- factor(data$Process, levels = data$Process[order(data$Score)])

# Create the plot

# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
ggplot(data, aes(x = Score, y = Process, fill = Score)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Enrichment Score (-log(P)*Z-score)", y = "", title = "")+theme(text = element_text(size = 18))+
  guides(fill = FALSE)



