####Supplementary figure S31a, b

library(harmony)
library(Seurat)
library(tidyverse)
library(ggplot2)
library(scDblFinder)
library(BiocParallel)
library(readxl)
library(SeuratWrappers)
library(Nebulosa)

fatbody_v.data_v_1 <- Read10X(data.dir = "./data/combined_virgin_bacterial_uninfected_VU")

sdata.1 <- CreateSeuratObject(fatbody_v.data_v_1, project = "Virgin",min.cells = 3, min.features = 200)

metadata<- sdata.1@meta.data
metadata$barcodes <- rownames(metadata)
barcodes <- data.frame(do.call('rbind', strsplit(as.character(metadata$barcodes),'-',fixed=TRUE)))
barcodes$group <- barcodes$X2
barcodes$group[barcodes$X2 == 1 ] <- "V_1"
barcodes$group[barcodes$X2 == 2 ] <- "V_2"
sdata.1@meta.data$Identity <- barcodes$group
sdata.1$type = "Virgin"
sdata.1 <- JoinLayers(sdata.1)

fatbody_v.data_v_2 <- Read10X(data.dir = "./data/combined_mated_bacterial_uninfected_MU")

sdata.2 <- CreateSeuratObject(fatbody_v.data_v_2, project = "Mated",min.cells = 3, min.features = 200)

metadata<- sdata.2@meta.data
metadata$barcodes <- rownames(metadata)
barcodes <- data.frame(do.call('rbind', strsplit(as.character(metadata$barcodes),'-',fixed=TRUE)))
barcodes$group <- barcodes$X2
barcodes$group[barcodes$X2 == 1 ] <- "M_1"
barcodes$group[barcodes$X2 == 2 ] <- "M_2"
sdata.2@meta.data$Identity <- barcodes$group

sdata.2$type = "Mated"
sdata.2 <- JoinLayers(sdata.2)



# Merge datasets into one single seurat object
fatbody_v <- merge(sdata.2, c(sdata.1), add.cell.ids = c("M","V"))
rm(fatbody_v.data_v_1,fatbody_v.data_v_2,sdata.1,sdata.2)
gc()
head(fatbody_v@meta.data, 10)
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
ElbowPlot(object = fatbody_v.filtered, ndims = 150)

fatbody_v.filtered <- RunUMAP(fatbody_v.filtered, reduction = "pca", dims = 1:50)

before <- DimPlot(fatbody_v.filtered, reduction = 'umap', group.by = 'type')
before


# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
set.seed(123) 
fatbody_v.filtered <- fatbody_v.filtered %>%
  RunUMAP(reduction = 'pca', dims = 1:50) %>%
  FindNeighbors(reduction = "pca", dims = 1:50) %>%
  FindClusters(resolution = 0.2)

# visualize 
after <- DimPlot(fatbody_v.filtered, reduction = 'umap', group.by = 'type')
before
after

DimPlot(fatbody_v.filtered, reduction = 'umap', group.by = 'type',  cols = c('Virgin' = 'red', 'Mated' = 'blue'))
DimPlot(fatbody_v.filtered, reduction = 'umap', group.by = 'type')

fatbody.markers <- FindAllMarkers(fatbody_v.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
fatbody.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

fatbody.markers %>%
  group_by(cluster) %>%
  top_n(n = 100, wt = avg_log2FC) -> top10_59

DimPlot(fatbody_v.filtered, reduction = 'umap', split.by = 'type')

fatbody_v.filtered <- FindClusters(fatbody_v.filtered, resolution = 0.1)

fatbody.markers <- FindAllMarkers(fatbody_v.filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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

names(new.cluster.ids) <- levels(fatbody_v.filtered)
fatbody_v.filtered <- RenameIdents(fatbody_v.filtered, new.cluster.ids)
fatbody_v.filtered[["cluster.ident"]] <- Idents(object = fatbody_v.filtered)
DimPlot(fatbody_v.filtered,reduction = "umap", pt.size = 0.5, label = F, label.size = 5)
###split by replicates.
DimPlot(fatbody_v.filtered,reduction = "umap", pt.size = 0.5, label = F, label.size = 5, split.by = "Identity")



fatbody_v<- subset(x = fatbody_v.filtered, subset = type == "Mated") ### change between mated and virgin conditions
##### making metadata
coldata <- fatbody_v@meta.data
coldata$barcodes <- rownames(coldata)
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


count_norm<-as.data.frame(raw_counts)

assoc_results <- data.frame(Pvalue = numeric())

for (i in 1:nrow(count_norm)) {
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),df_metadata$virus_NV, df_metadata$virus_DAV,df_metadata$identity, df_metadata$clusters)
  colnames(data) <- c("gene","NV","DAV","identity","clusters")
  data$condition <- as.factor(data$NV)
  data$condition <- as.factor(data$DAV)
  data$identity <- as.factor(data$identity)
  data$clusters <- as.factor(data$clusters)
  #try((ancova_model <- aov(gene ~ DAV*clusters + clusters+NV+identity, data = data)),silent = T)
  try((ancova_model <- aov(gene ~ NV+DAV+DAV*clusters+NV*clusters+clusters+identity, data = data)),silent = T) ## get p-value: 2 for NV (0/1), 3 for DAV (0/1), 4 for clusters (cell type), 5 for identity (replicates), 6 for DAV:clusters (cluster by DAV intereaction effect), 7 for NV:clusters (cluster by NV intereaction effect).
  try(d<- Anova(ancova_model, type="III"),silent = T) 
  assoc_results[length(assoc_results$Pvalue)+1, ] = c(d$`Pr(>F)`[6])
  print(i)
}


assoc_results$Pvalue <- as.numeric(assoc_results$Pvalue)

padj=p.adjust(assoc_results$Pvalue, method = "fdr")

conditions <- as.data.frame(df_metadata$virus_DAV) ###change it to "virus_NV" when calculating log2fold change for NV
conditions <- t(conditions)
conditions <- as.factor(conditions)

conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))

## foldchange here is only for NV uninfected vs infected cells and DAV uninfected vs infected cells.
## DAV*clusters and NV*clusters term here give  which genes have interaction effect (cluster by NV interaction effect andcluster by DAV interaction).
## later foldchange of interaction effect genes in DAV and NV infection by different cell types (cluster) are calculated using script ######



avona_results <- cbind(foldChanges,assoc_results, padj)
rownames(avona_results)=rownames(count_norm)
avona_results$gene <- row.names(avona_results)
avona_results$foldChanges<- round(avona_results$foldChanges,2)

avona_results$gene[avona_results$gene == 'JX220408.1-Nora-virus-isolate-FR1'] <- 'NV'
#avona_results <- avona_results[!(row.names(avona_results) %in% c("JX220408.1-Nora-virus-isolate-FR1")),]
avona_results <- avona_results[!(row.names(avona_results) %in% c("Mf")),]


# Subset the significant results
sig_res_anova <- dplyr::filter(avona_results, padj < 0.05, abs(foldChanges)>0.5) %>%
  dplyr::arrange(padj)

#write.csv(sig_res_anova, "DEGsUV_controlling_cluster.csv", quote=F, row.names = T,col.names = T)
avona_results[avona_results == -Inf] <- NA
avona_results <- na.omit(avona_results)

## we compare this results with the one with harmony integration and make the venn diagram