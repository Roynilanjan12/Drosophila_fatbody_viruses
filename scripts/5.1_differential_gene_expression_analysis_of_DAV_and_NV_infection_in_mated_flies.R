##Figure 4b
##Supplementary figure S18, S19, S22, S23, S24

###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first

# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v<- subset(x = fatbody_v.harmony, subset = type == "Mated")
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

# ---- Model fitting (ANOVA/ANCOVA) ------------------------------------------------------
# The model below tests main effects and interaction terms; Type-III sums of squares
# are used to obtain p-values when multiple factors are present.
# ----------------------------------------------------------------------------
  try((ancova_model <- aov(gene ~ NV+DAV+DAV*clusters+NV*clusters+clusters+identity, data = data)),silent = T) ## get p-value: 2 for NV (0/1), 3 for DAV (0/1), 4 for clusters (cell type), 5 for identity (replicates), 6 for DAV:clusters (cluster by DAV intereaction effect), 7 for NV:clusters (cluster by NV intereaction effect).

# ---- Type-III ANOVA table (car::Anova) ------------------------------------------------------
# The model below tests main effects and interaction terms; Type-III sums of squares
# are used to obtain p-values when multiple factors are present.
# ----------------------------------------------------------------------------
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



# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(readxl)
library(EnhancedVolcano)
library(ggpubr)
library(IDPmisc)

### DAV in mated

volcano <- read_excel("/data/DEGs_DAV_NV_mated_condition.xlsx", sheet = "DAV_MU")
volcano <- as.data.frame(volcano)
go <- filter(volcano,padj<0.05, abs(foldChanges)>0.5)

immune_set <- read.csv("/data/List_of_immune_genes_updated.csv", header=T)
immune_set$gene <- immune_set$Symbol
immune_go <- merge(go, immune_set, by = "gene", all.x = TRUE)
immune_go<-NaRV.omit(immune_go)
ggbarplot(immune_go, x = "gene", y = "foldChanges",
          fill = "Immune.Process",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE, 
          rotate = TRUE,
          ggtheme = theme_minimal()
          
)+theme(text = element_text(size = 15))+ylab("Log2FC")


# Create a vector of colors based on conditions
keyvals <- ifelse(volcano$foldChanges < -0.5 & volcano$padj < 0.05, 'royalblue',
                  ifelse(volcano$foldChanges > 0.5 & volcano$padj < 0.05, 'red2',
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
                selectLab = c('DAV','ref(2)P',
                              "NimB1",
                              "DptA",
                              "Vago","DNApol-alpha73","spd-2","Nup154","CG17224","CG13641","CG33926"),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = '',
                subtitle = "",
                pCutoff = 0.05,
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
                boxedLabels = TRUE,
                maxoverlapsConnectors = Inf)




# Load necessary library

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(ggplot2)
## base on pvalue <1e-05 (for DAV in mated flies) from ST7_DEGs_DAV_NV_mated_condition.xlsx data GO enrichment analysis was
## done using flyenricher webtool. The result from the flyenricher webtool is below for plotting


# Recreate the data from the shown plot
data_plot2 <- data.frame(
  Process = c(
    "cytoplasmic translation (GO:0002181)",
    "translation (GO:0006412)",
    "ribosome assembly (GO:0042255)",
    "ATP synthesis coupled proton transport (GO:0015986)",
    "ribosomal large subunit assembly (GO:0000027)",
    "proton transport (GO:0015992)",
    "energy coupled proton transport, down electrochemical gradient (GO:0015985)",
    "hydrogen ion transmembrane transport (GO:1902600)",
    "ribosomal small subunit assembly (GO:0000028)",
    "ATP biosynthetic process (GO:0006754)"
  ),
  Score = c(
    270,  # cytoplasmic translation
    180,  # translation
    70,   # ribosome assembly
    65,   # ATP synthesis coupled proton transport
    55,   # ribosomal large subunit assembly
    52,   # proton transport
    50,   # energy coupled proton transport
    48,   # hydrogen ion transmembrane transport
    47,   # ribosomal small subunit assembly
    45    # ATP biosynthetic process
  )
)

# Reorder factor levels by score
data_plot2$Process <- factor(
  data_plot2$Process,
  levels = data_plot2$Process[order(data_plot2$Score)]
)

# Plot (same style as before)
library(ggplot2)


# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
ggplot(data_plot2, aes(x = Score, y = Process, fill = Score)) +
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


## For NV in mated
volcano <- read_excel("/data/DEGs_DAV_NV_mated_condition.xlsx", sheet = "NV_MU")
volcano <- as.data.frame(volcano)
go <- filter(volcano,padj<0.05, abs(foldChanges)>0.5)

immune_set <- read.csv("/data/List_of_immune_genes_updated.csv", header=T)
immune_set$gene <- immune_set$Symbol
immune_go <- merge(go, immune_set, by = "gene", all.x = TRUE)
immune_go<-NaRV.omit(immune_go)
ggbarplot(immune_go, x = "gene", y = "foldChanges",
          fill = "Immune.Process",               # change fill color by cyl
          color = "white",            # Set bar border colors to white
          # jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE, 
          rotate = TRUE,
          ggtheme = theme_minimal()
          
)+theme(text = element_text(size = 15))+ylab("Log2FC")

# Create a vector of colors based on conditions
keyvals <- ifelse(volcano$foldChanges < -0.5 & volcano$padj < 0.05, 'royalblue',
                  ifelse(volcano$foldChanges > 0.5 & volcano$padj < 0.05, 'red2',
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
                selectLab = c('CYLD','CG32039','Atg13',"eater","CG4461","Karl", "Pxn",
                              "MKP-4","Dredd",
                              "ft", "Myo95E",
                              "TotA","CG9008","Mvb12",
                              "Sbf","Traf4","sing","NV","NimB3"),
                xlab = bquote(~Log[2]~ 'fold change'),
                title = '',
                subtitle = "",
                pCutoff = 0.05,
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
                boxedLabels = TRUE,
                maxoverlapsConnectors = Inf)



## base on pvalue <1e-05 (for NV in mated flies) from ST7_DEGs_DAV_NV_mated_condition.xlsx data GO enrichment analysis was
## done using flyenricher webtool. The result from the flyenricher webtool is below for plotting

# Recreate the data from the plot
data_plot <- data.frame(
  Process = c(
    "cytoplasmic translation (GO:0002181)",
    "translation (GO:0006412)",
    "ATP synthesis coupled proton transport (GO:0015986)",
    "proton transport (GO:0015992)",
    "energy coupled proton transport, down electrochemical gradient (GO:0015985)",
    "hydrogen ion transmembrane transport (GO:1902600)",
    "ATP biosynthetic process (GO:0006754)",
    "cellular protein complex disassembly (GO:0043624)",
    "muscle contraction (GO:0006936)",
    "negative regulation of MAP kinase activity (GO:0043407)"
  ),
  Score = c(
    148,
    105,
    55,
    52,
    50,
    48,
    46,
    36,
    35,
    34
  )
)

# Reorder factor levels by score
data_plot$Process <- factor(
  data_plot$Process,
  levels = data_plot$Process[order(data_plot$Score)]
)

# Plot

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(ggplot2)


# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
ggplot(data_plot, aes(x = Score, y = Process, fill = Score)) +
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

