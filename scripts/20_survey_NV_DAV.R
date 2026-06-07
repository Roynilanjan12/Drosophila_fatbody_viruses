##Supplementary figure S3

library(tidyverse)

# 1. Load Data
RPM_fatbody_virus <- read.table(
  ".data/viral_cov_RPM100_fatbody.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

bioproject <- read_csv("/Users/nilanjan/Downloads/Dr._Unckless_Lab/Paper Manuscript/Codes/final/first_round_review/bioproject.csv")
bioproject$SRR <- bioproject$Run

# 2. Process Viral Data and Join BioProject info
RPM_fatbody_virus <- RPM_fatbody_virus %>%
  mutate(
    virus = case_when(
      grepl("Nora", rname, ignore.case = TRUE) ~ "Nora virus",
      grepl("Drosophila_A", rname, ignore.case = TRUE) ~ "Drosophila A virus",
      TRUE ~ NA_character_
    ),
    log10_RPM = log10(RPM)
  ) %>%
  filter(!is.na(virus), is.finite(log10_RPM)) %>%
  left_join(bioproject %>% dplyr::select(SRR, BioProject), by = "SRR")

# 3. Calculate Stats per Sample (SRR)
srr_counts <- RPM_fatbody_virus %>%
  group_by(SRR, BioProject) %>%
  summarise(
    has_DAV = any(virus == "Drosophila A virus"),
    has_NV  = any(virus == "Nora virus"),
    .groups = "drop"
  ) %>%
  mutate(group = case_when(
    has_DAV & has_NV ~ "Both",
    has_DAV          ~ "Drosophila A virus",
    has_NV           ~ "Nora virus"
  ))

# 4. Summarize and build the final dataframe with ordered labels
plot_stats <- srr_counts %>%
  group_by(group) %>%
  summarise(
    n = n(),
    n_proj = n_distinct(BioProject),
    .groups = "drop"
  )

# Add Uninfected
uninfected_n <- 1362

df_final <- plot_stats %>%
  add_row(group = "Uninfected", n = uninfected_n) %>%
  mutate(
    prop = n / sum(n),
    # Generate the dynamic label string
    legend_label_str = case_when(
      group == "Uninfected" ~ paste0(group, " — ", n, " (", round(100 * prop, 1), "%)"),
      TRUE ~ paste0(group, " — ", n, " (", round(100 * prop, 1), "%), ", n_proj, " BioProjects")
    )
  )

# --- THE KEY FIX: FACTOR ORDERING ---
# We extract the labels in the specific group order you requested
ordered_labels <- df_final$legend_label_str[match(
  c("Drosophila A virus", "Nora virus", "Both", "Uninfected"), 
  df_final$group
)]

df_final <- df_final %>%
  mutate(legend_label = factor(legend_label_str, levels = ordered_labels))

# 5. Create a named color vector matching the factor levels
plot_colors <- setNames(
  c("red", "blue", "black", "grey80"),
  ordered_labels
)

# 6. Plot
ggplot(df_final, aes(x = "", y = prop, fill = legend_label)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = plot_colors) +
  labs(fill = "Infection status") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10)
  )

# compute means
mean_df <- RPM_fatbody_virus %>%
  group_by(virus) %>%
  summarise(mean_log10_RPM = mean(log10_RPM, na.rm = TRUE))

ggplot(RPM_fatbody_virus, aes(x = log10_RPM, fill = virus)) +
  geom_density(alpha = 0.5) +
  geom_vline(
    data = mean_df,
    aes(xintercept = mean_log10_RPM, color = virus),
    linetype = "dashed",
    linewidth = 1,
    show.legend = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Drosophila A virus" = "red",
      "Nora virus" = "blue"
    )
  ) +
  scale_color_manual(
    values = c(
      "Drosophila A virus" = "red",
      "Nora virus" = "blue"
    )
  ) +
  labs(
    x = expression(log[10]~"(RPM)"),
    y = "Density",
    fill = "Virus"
  ) +
  theme_classic()


