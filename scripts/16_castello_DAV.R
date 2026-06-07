##supplementary figure S33 (b, b')
library(readxl)
library(dplyr)
library(org.Dm.eg.db)
library(AnnotationDbi)
library(ggpubr)

# --- Read + map Castello FBgn -> gene symbol ---
castello <- read_excel(".data/pbio.3003437.s003.xlsx")

castello <- castello %>%
  dplyr::mutate(
    gene_symbol = mapIds(
      org.Dm.eg.db,
      keys = Id,
      column = "SYMBOL",
      keytype = "FLYBASE",
      multiVals = "first"
    )
  )

castello$gene <- castello$gene_symbol
castello$gene_symbol <- NULL

# (recommended) ensure numeric columns are numeric (readxl sometimes imports as character)
castello$log2FoldChange <- as.numeric(castello$log2FoldChange)
castello$padj <- as.numeric(castello$padj)
castello$DPE <- as.numeric(castello$DPE)

# --- Read your DEG file ---
my <- read_excel("./supplementary_files/ST7_DEGs_DAV_NV_mated_condition.xlsx", sheet = "DAV_MU")
my$...1 <- NULL
# your file has foldChanges; create log2FoldChange
my$log2FoldChange <- as.numeric(my$foldChanges)
my$padj <- as.numeric(my$padj)

# filter DEGs
#my <- my %>%
    #dplyr::filter(padj < 0.05 & abs(log2FoldChange) > 0.5)
my <- my %>%
  dplyr::filter(padj < 0.05)
castello <- castello %>%
  dplyr::filter(padj < 0.05)

# filter castello to overlapping genes and condition
castello <- castello %>%
  dplyr::filter(gene %in% my$gene)

castello <- castello %>%
  dplyr::filter(Treatment == "DAV" & DPE == 12)  ## 1, 12, 25

# define castello_sub (you used this name later)
castello_sub <- castello

# --- Merge on gene (inner join keeps overlap only) ---
merged_df <- dplyr::inner_join(
  castello_sub %>% dplyr::select(gene, log2FoldChange_castello = log2FoldChange),
  my %>% dplyr::select(gene, log2FoldChange_my = log2FoldChange),
  by = "gene"
)

# --- Plot: Castello log2FC vs My log2FC ---
p <- ggscatter(
  merged_df,
  x = "log2FoldChange_castello",
  y = "log2FoldChange_my",
  add = "reg.line",
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "spearman",   # change to "pearson" if you want
  xlab = "Castello et al., 2025, Log2 Fold Change in \n Drosophila A virus infection, 12 day post infection)",
  ylab = "Log2 Fold Change in \n Drosophila A virus infection \n (this study)",
  size = 2
) +
  theme_classic()+theme(text = element_text(size = 16))

p

library(dplyr)
library(ggplot2)

# -----------------------------
# 1️⃣ Define threshold
# -----------------------------
threshold <- 0.5

# -----------------------------
# 2️⃣ Create logical significance flags
# -----------------------------
merged_df <- merged_df %>%
  mutate(
    sig_castello = abs(log2FoldChange_castello) > threshold,
    sig_my = abs(log2FoldChange_my) > threshold
  )

# -----------------------------
# 3️⃣ Build 2×2 contingency matrix manually
# -----------------------------
a <- sum( merged_df$sig_castello &  merged_df$sig_my )   # both significant
b <- sum( merged_df$sig_castello & !merged_df$sig_my )   # castello only
c <- sum(!merged_df$sig_castello &  merged_df$sig_my )   # ours only
d <- sum(!merged_df$sig_castello & !merged_df$sig_my )   # neither

fisher_matrix <- matrix(
  c(d, c,
    b, a),
  nrow = 2,
  byrow = TRUE
)

rownames(fisher_matrix) <- c("Castello_NotSig", "Castello_Sig")
colnames(fisher_matrix) <- c("Ours_NotSig", "Ours_Sig")

print(fisher_matrix)

# -----------------------------
# 4️⃣ Fisher’s Exact Test
# -----------------------------
fisher_res <- fisher.test(fisher_matrix)

p_value <- fisher_res$p.value
odds_ratio <- fisher_res$estimate

print(fisher_res)

# -----------------------------
# 5️⃣ Convert matrix for plotting
# -----------------------------
plot_df <- as.data.frame(as.table(fisher_matrix))
colnames(plot_df) <- c("Castello", "Ours", "Count")

# Clean labels for axes
plot_df$Castello <- ifelse(plot_df$Castello == "Castello_Sig",
                           "Significant",
                           "Not Significant")

plot_df$Ours <- ifelse(plot_df$Ours == "Ours_Sig",
                       "Significant",
                       "Not Significant")

# Optional: control order (NotSig first, Sig second)
plot_df$Castello <- factor(plot_df$Castello,
                           levels = c("Not Significant", "Significant"))

plot_df$Ours <- factor(plot_df$Ours,
                       levels = c("Not Significant", "Significant"))

# -----------------------------
# 6️⃣ Plot contingency table (numbers only)
# -----------------------------
ggplot(plot_df, aes(x = Castello, y = Ours)) +
  geom_tile(fill = "white", color = "black") +
  geom_text(aes(label = Count), size = 7) +
  labs(
    subtitle = paste0(
      "Significant if |log2FC| > ", threshold, # Added comma, removed extra quote
      "\nFisher Exact Test: p = ", signif(p_value, 3)
    ),
    x = "Castello et al., 2025 \n Drosophila A virus infection, 12 day post infection",
    y = "This Study"
  ) +
  theme_classic() + theme(text = element_text(size = 16))



