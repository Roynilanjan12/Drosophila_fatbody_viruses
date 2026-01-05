##### Figure 4c
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first
# --- Packages -----------------------------------------------------------

# ---- Dependencies ------------------------------------------------------
# ----------------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)  # to stitch per-gene panels together
library(clipr)

# --- Data subsets -------------------------------------------------------

# ---- Filter/subset cells or features ------------------------------------------------------
# ----------------------------------------------------------------------------
fatbody_v <- subset(x = fatbody_v.harmony, subset = type == "Mated")

uninfected <- subset(
  x = fatbody_v,
  subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 &
    `KP969946.1-Drosophila-A-virus-isolate-LJ35` == 0
)

DAV <- subset(
  x = fatbody_v,
  subset = `JX220408.1-Nora-virus-isolate-FR1` == 0 &
    `KP969946.1-Drosophila-A-virus-isolate-LJ35` > 0
)

NV <- subset(
  x = fatbody_v,
  subset = `JX220408.1-Nora-virus-isolate-FR1` > 0 &
    `KP969946.1-Drosophila-A-virus-isolate-LJ35` == 0
)

Both <- subset(
  x = fatbody_v,
  subset = `JX220408.1-Nora-virus-isolate-FR1` > 0 &
    `KP969946.1-Drosophila-A-virus-isolate-LJ35` > 0
)

# --- Gene lists & order -------------------------------------------------
genes_reception    <- c("spz","Tl")
genes_transduction <- c("Myd88","pll","cact","Dif","dl")
genes_response     <- c("Drs","Mtk","IM14","IM3")
genes_all <- c(genes_reception, genes_transduction, genes_response)

# --- Helper: fetch counts to long format --------------------------------
fetch_long <- function(obj, condition_label) {

# ---- Extract expression/metadata for downstream analyses ------------------------------------------------------
# ----------------------------------------------------------------------------
  dd <- FetchData(object = obj, vars = genes_all, layer = "data")
  dd$cell <- rownames(dd)
  dd %>%
    pivot_longer(cols = all_of(genes_all),
                 names_to = "gene", values_to = "count") %>%
    mutate(condition = condition_label)
}

df <- bind_rows(
  fetch_long(uninfected, "Uninfected cells"),
  fetch_long(DAV,         "DAV infected cells"),
  fetch_long(NV,          "NV infected cells"),
  fetch_long(Both,        "Both")
)

# --- Summary stats (mean ± SE) ------------------------------------------
sum_df <- df %>%
  group_by(condition, gene) %>%
  summarise(mean = mean(count), se = sd(count)/sqrt(n()), .groups = "drop")

sum_df$condition <- factor(
  sum_df$condition,
  levels = c("Uninfected cells", "DAV infected cells", "NV infected cells", "Both")
)

# --- Small plotting function (one panel per gene, no legend) ------------
make_gene_plot <- function(gene_name, show_y_title = FALSE) {

# ---- Plotting with ggplot2 ------------------------------------------------------
# ----------------------------------------------------------------------------
  ggplot(sum_df %>% filter(gene == gene_name),
         aes(x = gene, y = mean, color = condition)) +
    geom_point(position = position_dodge(width = 0.6), size = 3) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  width = 0.35, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = c("#999999", "red", "blue", "black")) +
    labs(x = NULL, y = if (show_y_title) "Mean counts" else NULL) +
    theme_classic() +
    theme(
      text = element_text(size = 22),
      axis.text.x = element_text(size = 22, angle = 0, hjust = 0.5, vjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 20, face = "bold"),
      axis.title.y = element_text(size = 22, margin = margin(r = 10)),
      legend.position = "none"
    )
}

# --- Plot 1: Reception + Transduction -----------------------------------
genes_plot1 <- c(genes_reception, genes_transduction)
plot_list1 <- lapply(seq_along(genes_plot1), function(i) {
  make_gene_plot(genes_plot1[i], show_y_title = (i == 1))
})

# ---- Combine panels (patchwork) ------------------------------------------------------
# ----------------------------------------------------------------------------
plot1 <- wrap_plots(plot_list1, nrow = 1)

# --- Plot 2: Response ---------------------------------------------------
genes_plot2 <- genes_response
plot_list2 <- lapply(seq_along(genes_plot2), function(i) {
  make_gene_plot(genes_plot2[i], show_y_title = (i == 1))
})
plot2 <- wrap_plots(plot_list2, nrow = 1)

# --- Show the two plots -------------------------------------------------
plot1
plot2

# Optional: export
# ggsave("reception_transduction_no_legend.png", plot1, width = 14, height = 4, dpi = 300)
# ggsave("response_no_legend.png",              plot2, width = 8,  height = 4, dpi = 300)

# Optional: copy stats table
clipr::write_clip(sum_df)
