#supplementary figure S22
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first

## Permutation Analysis for ONE gene (Null Distribution Check)


suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

set.seed(123)
# ---- Filter/subset cells or features
fatbody_v<- subset(x = fatbody_v.harmony, subset = Identity == "V_2")

# SETTINGS

GENE    <- "AttC"  #change the genes accordingly to DEGs
N_PERM  <- 10000   # Number of times to shuffle labels
LAYER   <- "data" # Normalized data

# -------------------------
# 0) Prepare Data
# -------------------------
# Using your existing subsetted replicate (fatbody_v)
# We need to know which cells are NV+ and which are NV-
total_cells <- FetchData(fatbody_v, vars = c(GENE, "JX220408.1-Nora-virus-isolate-FR1"), layer = LAYER)
colnames(total_cells) <- c("expr", "virus")

# Define the "True" status
total_cells$status <- ifelse(total_cells$virus > 0, "Infected", "Uninfected")

# Calculate OBSERVED Log2FC
obs_inf   <- mean(total_cells$expr[total_cells$status == "Infected"])
obs_uninf <- mean(total_cells$expr[total_cells$status == "Uninfected"])
# Adding a small pseudocount to avoid division by zero/log(0) if using raw counts
obs_l2fc  <- log2((obs_inf + 0.01) / (obs_uninf + 0.01))

#Permutation Loop

perm_l2fc <- numeric(N_PERM)
n_inf     <- sum(total_cells$status == "Infected")

for (i in seq_len(N_PERM)) {
  # SHUFFLE the labels: randomly pick n_inf cells to be "Infected"
  shuffled_indices <- sample(seq_len(nrow(total_cells)), n_inf)
  
  perm_inf   <- mean(total_cells$expr[shuffled_indices])
  perm_uninf <- mean(total_cells$expr[-shuffled_indices])
  
  perm_l2fc[i] <- log2((perm_inf + 0.01) / (perm_uninf + 0.01))
}


# Calculate P-value (Empirical)

# The fraction of times the random L2FC was more extreme than our real L2FC
p_val <- mean(abs(perm_l2fc) >= abs(obs_l2fc))

cat("Observed L2FC:", obs_l2fc, "\n")
cat("Empirical p-value:", p_val, "\n")


#Bell curve plot

perm_df <- data.frame(l2fc = perm_l2fc)

# Dynamically set limits based on whichever is further from zero
max_val <- max(abs(c(perm_l2fc, obs_l2fc))) * 1.1

ggplot(perm_df, aes(x = l2fc)) +
  geom_histogram(
    bins = 50, 
    fill = "grey70", 
    color = "white"
  ) +
  geom_vline(
    xintercept = obs_l2fc, 
    color = "red", 
    linewidth = 1.2
  ) +
  # This ensures the red line is never on the "edge" of the image
  coord_cartesian(xlim = c(-max_val, max_val)) + 
  labs(
    title = paste0(
      GENE),
    x = "Log2 Fold Change",
    y = "Frequency"
  ) +
  theme_classic(base_size = 18)
