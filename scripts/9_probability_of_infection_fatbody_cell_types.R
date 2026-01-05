### Supplementary figure S6

###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first
# ---- Load required R packages ----
# Packages are loaded explicitly to make dependencies clear for reproduction.
library(dplyr)
library(purrr)

#### Mated flies
meta<- fatbody_v.harmony@meta.data
# ---- Feature extraction ----
# Fetch expression values for specific genes/TEs/viral features from the Seurat object.
DAV<-FetchData(fatbody_v.harmony, vars = "KP969946.1-Drosophila-A-virus-isolate-LJ35") 
NV<-FetchData(fatbody_v.harmony, vars = "JX220408.1-Nora-virus-isolate-FR1") 
meta$DAV<- DAV$`KP969946.1-Drosophila-A-virus-isolate-LJ35`
meta$NV<- NV$`JX220408.1-Nora-virus-isolate-FR1`
meta <- meta[ , c("Identity", "type","cluster.ident","DAV","NV")]
meta$DAV_status <- ifelse(meta$DAV > 0, "Infected", "Uninfected")
meta$NV_status <- ifelse(meta$NV > 0, "Infected", "Uninfected")
meta <- filter(meta, type=="Mated")
meta$DAV <- NULL
meta$NV <- NULL

### for DAV (mated flies)
meta2 <- meta %>%
  # add a numeric NV flag
  mutate(NV_flag = ifelse(NV_status == "Uninfected", 0, 1))

# 1) Summarize by Identity, cluster, and NV_flag
dav <- meta2 %>%
  group_by(Identity, cluster.ident, NV_flag) %>%
  summarize(
    DAV_infected   = sum(DAV_status   == "Infected",   na.rm = TRUE),
    DAV_uninfected = sum(DAV_status   == "Uninfected", na.rm = TRUE),
    .groups = "drop"
  )

dav$Identity <- as.factor(dav$Identity)
dav$NV_flag <- as.factor(dav$NV_flag)


# 2) Fit the logistic‐regression (binomial) model
# ---- Statistical modeling (GLM) ----
# Fit generalized linear models (binomial/logit) to quantify infection probability differences.
mod <- glm(cbind(DAV_infected, DAV_uninfected) ~ cluster.ident+NV_flag+Identity,
           data = dav,
           family = binomial(link = "logit"))
summary(mod)

# 3) Overall test of cluster effect via likelihood‐ratio (Chi‐square)
anova(mod, test = "Chisq")


# 4) Post-hoc Tukey comparisons

library(multcomp)

tuk <- glht(mod,
            linfct = mcp(cluster.ident = "Tukey"))
summary(tuk)  


#######plotting

# 1) Run your Tukey post-hoc as before
#tuk <- glht(mod, linfct = mcp(cluster.ident = "Tukey"))
tuk_sum <- summary(tuk)

# 2) Extract contrast labels, estimates and p-values
comp  <- names(tuk_sum$test$coefficients)
est   <- tuk_sum$test$coefficients
pval  <- tuk_sum$test$pvalues

# 3) Clean up the contrast names (drop the " == 0")
comp_clean <- gsub(" == 0", "", comp)
parts      <- strsplit(comp_clean, " - ")

# 4) Build a tidy data.frame of pairwise contrasts
tuk_df <- data.frame(
  cluster1 = sapply(parts, `[`, 2),  # reference
  cluster2 = sapply(parts, `[`, 1),  # comparison
  estimate = est,
  p.value  = pval,
  stringsAsFactors = FALSE
) %>%
  # 5) Add significance stars
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

rownames(tuk_df) <- 1:dim(tuk_df)[1]

# assume you already have your desired cluster_order vector defined:
cluster_order <- c(
  "Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
  "Epithelial cells", "Muscle cells", "Oenocyte",
  "Hemocytes", "Female reproductive system",
  "Unknown 1", "Sensory neuron", "Unknown 2"
)

# build a diagonal of zeros / NAs / empty sig
diag_df <- tibble(
  cluster1 = cluster_order,
  cluster2 = cluster_order,
  estimate = 0,
  p.value  = NA_real_,
  sig      = ""
)

# bind it on to tuk_df
tuk_df <- bind_rows(tuk_df, diag_df) %>% 
  # (optional) if you want to keep only one row per pair in case any duplicates:
  distinct(cluster1, cluster2, .keep_all = TRUE)


# 6) Mirror to get a full square matrix
tuk_full <- bind_rows(
  tuk_df,
  tuk_df %>% transmute(
    cluster1 = cluster2,
    cluster2 = cluster1,
  )
)

# 7) Define your desired ordering
cluster_order <- c(
  "Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
  "Epithelial cells", "Muscle cells", "Oenocyte",
  "Hemocytes", "Female reproductive system",
  "Unknown 1", "Sensory neuron", "Unknown 2"
)

# 8) Re-factor both axes
tuk_full <- tuk_full %>%
  mutate(
    cluster1 = factor(cluster1, levels = cluster_order),
    cluster2 = factor(cluster2, levels = cluster_order)
  )

tuk_full <- tuk_full %>%
  filter(!is.na(estimate))

# 9) Plot with reference (cluster1) on y, comparison (cluster2) on x
##Pairwise comparison of DAV infection probability\nin the fat body cell types (Mated)
tuk_full <- tuk_full %>%
  filter(
    cluster1 %in% c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3"),
    cluster2 %in% c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3")
  )
# ---- Visualization ----
# Figures are generated with ggplot2; styling is kept consistent with manuscript themes.
ggplot(tuk_full, aes(x = cluster2, y = cluster1, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = sig), vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(
    midpoint = 0,
    low      = "blue",
    mid      = "grey95",
    high     = "red",
    name     = "Log-odds\ndifference"
  ) +
  labs(
    x     = "Compared Cluster",
    y     = "Reference Cluster",
    title = "DAV (Mated)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title  = element_text(face = "bold")
  )





############ NV (mated flies)
meta2 <- meta %>%
  # add a numeric NV flag
  mutate(DAV_flag = ifelse(DAV_status == "Uninfected", 0, 1))

# 1) Summarize by Identity, cluster, and NV_flag
nv <- meta2 %>%
  group_by(Identity, cluster.ident, DAV_flag) %>%
  summarize(
    NV_infected   = sum(NV_status   == "Infected",   na.rm = TRUE),
    NV_uninfected = sum(NV_status   == "Uninfected", na.rm = TRUE),
    .groups = "drop"
  )

nv$Identity <- as.factor(nv$Identity)
nv$DAV_flag <- as.factor(nv$DAV_flag)


# 2) Fit the logistic‐regression (binomial) model
mod <- glm(cbind(NV_infected, NV_uninfected) ~ cluster.ident+DAV_flag+Identity,
           data = nv,
           family = binomial(link = "logit"))
summary(mod)

# 3) Overall test of cluster effect via likelihood‐ratio (Chi‐square)
anova(mod, test = "Chisq")


# 4) Post-hoc Tukey comparisons

library(multcomp)

tuk <- glht(mod,
            linfct = mcp(cluster.ident = "Tukey"))
summary(tuk)  


#######plotting

# 1) Run your Tukey post-hoc as before
#tuk <- glht(mod, linfct = mcp(cluster.ident = "Tukey"))
tuk_sum <- summary(tuk)

# 2) Extract contrast labels, estimates and p-values
comp  <- names(tuk_sum$test$coefficients)
est   <- tuk_sum$test$coefficients
pval  <- tuk_sum$test$pvalues

# 3) Clean up the contrast names (drop the " == 0")
comp_clean <- gsub(" == 0", "", comp)
parts      <- strsplit(comp_clean, " - ")

# 4) Build a tidy data.frame of pairwise contrasts
tuk_df <- data.frame(
  cluster1 = sapply(parts, `[`, 2),  # reference
  cluster2 = sapply(parts, `[`, 1),  # comparison
  estimate = est,
  p.value  = pval,
  stringsAsFactors = FALSE
) %>%
  # 5) Add significance stars
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

rownames(tuk_df) <- 1:dim(tuk_df)[1]

# assume you already have your desired cluster_order vector defined:
cluster_order <- c(
  "Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
  "Epithelial cells", "Muscle cells", "Oenocyte",
  "Hemocytes", "Female reproductive system",
  "Unknown 1", "Sensory neuron", "Unknown 2"
)

# build a diagonal of zeros / NAs / empty sig
diag_df <- tibble(
  cluster1 = cluster_order,
  cluster2 = cluster_order,
  estimate = 0,
  p.value  = NA_real_,
  sig      = ""
)

# bind it on to tuk_df
tuk_df <- bind_rows(tuk_df, diag_df) %>% 
  # (optional) if you want to keep only one row per pair in case any duplicates:
  distinct(cluster1, cluster2, .keep_all = TRUE)


# 6) Mirror to get a full square matrix
tuk_full <- bind_rows(
  tuk_df,
  tuk_df %>% transmute(
    cluster1 = cluster2,
    cluster2 = cluster1,
  )
)

# 7) Define your desired ordering
cluster_order <- c(
  "Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
  "Epithelial cells", "Muscle cells", "Oenocyte",
  "Hemocytes", "Female reproductive system",
  "Unknown 1", "Sensory neuron", "Unknown 2"
)

# 8) Re-factor both axes
tuk_full <- tuk_full %>%
  mutate(
    cluster1 = factor(cluster1, levels = cluster_order),
    cluster2 = factor(cluster2, levels = cluster_order)
  )

tuk_full <- tuk_full %>%
  filter(!is.na(estimate))

# 9) Plot with reference (cluster1) on y, comparison (cluster2) on x
##Pairwise comparison of NV infection probability\nin the fat body cell types (Mated)
tuk_full <- tuk_full %>%
  filter(
    cluster1 %in% c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3"),
    cluster2 %in% c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3")
  )
ggplot(tuk_full, aes(x = cluster2, y = cluster1, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = sig), vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(
    midpoint = 0,
    low      = "blue",
    mid      = "grey95",
    high     = "red",
    name     = "Log-odds\ndifference"
  ) +
  labs(
    x     = "Compared Cluster",
    y     = "Reference Cluster",
    title = "NV (Mated)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title  = element_text(face = "bold")
  )





####virgin flies
meta<- fatbody_v.harmony@meta.data
DAV<-FetchData(fatbody_v.harmony, vars = "KP969946.1-Drosophila-A-virus-isolate-LJ35") 
NV<-FetchData(fatbody_v.harmony, vars = "JX220408.1-Nora-virus-isolate-FR1") 
meta$DAV<- DAV$`KP969946.1-Drosophila-A-virus-isolate-LJ35`
meta$NV<- NV$`JX220408.1-Nora-virus-isolate-FR1`
meta <- meta[ , c("Identity", "type","cluster.ident","DAV","NV")]
meta$DAV_status <- ifelse(meta$DAV > 0, "Infected", "Uninfected")
meta$NV_status <- ifelse(meta$NV > 0, "Infected", "Uninfected")
meta <- filter(meta, type=="Virgin")
meta$DAV <- NULL
meta$NV <- NULL


############ NV (virgin flies)
meta2 <- meta

# 1) Summarize by Identity, cluster, and NV_flag
nv <- meta2 %>%
  group_by(Identity, cluster.ident) %>%
  summarize(
    NV_infected   = sum(NV_status   == "Infected",   na.rm = TRUE),
    NV_uninfected = sum(NV_status   == "Uninfected", na.rm = TRUE),
    .groups = "drop"
  )

nv$Identity <- as.factor(nv$Identity)

# 2) Fit the logistic‐regression (binomial) model
mod <- glm(cbind(NV_infected, NV_uninfected) ~ cluster.ident+Identity,
           data = nv,
           family = binomial(link = "logit"))
summary(mod)

# 3) Overall test of cluster effect via likelihood‐ratio (Chi‐square)
anova(mod, test = "Chisq")


# 4) Post-hoc Tukey comparisons

library(multcomp)

tuk <- glht(mod,
            linfct = mcp(cluster.ident = "Tukey"))
summary(tuk)  

#######plotting

# 1) Run your Tukey post-hoc as before
#tuk <- glht(mod, linfct = mcp(cluster.ident = "Tukey"))
tuk_sum <- summary(tuk)

# 2) Extract contrast labels, estimates and p-values
comp  <- names(tuk_sum$test$coefficients)
est   <- tuk_sum$test$coefficients
pval  <- tuk_sum$test$pvalues

# 3) Clean up the contrast names (drop the " == 0")
comp_clean <- gsub(" == 0", "", comp)
parts      <- strsplit(comp_clean, " - ")

# 4) Build a tidy data.frame of pairwise contrasts
tuk_df <- data.frame(
  cluster1 = sapply(parts, `[`, 2),  # reference
  cluster2 = sapply(parts, `[`, 1),  # comparison
  estimate = est,
  p.value  = pval,
  stringsAsFactors = FALSE
) %>%
  # 5) Add significance stars
  mutate(sig = case_when(
    p.value < 0.001 ~ "***",
    p.value < 0.01  ~ "**",
    p.value < 0.05  ~ "*",
    TRUE            ~ ""
  ))

rownames(tuk_df) <- 1:dim(tuk_df)[1]

# assume you already have your desired cluster_order vector defined:
cluster_order <- c(
  "Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
  "Epithelial cells", "Muscle cells", "Oenocyte",
  "Hemocytes", "Female reproductive system",
  "Unknown 1", "Sensory neuron", "Unknown 2"
)

# build a diagonal of zeros / NAs / empty sig
diag_df <- tibble(
  cluster1 = cluster_order,
  cluster2 = cluster_order,
  estimate = 0,
  p.value  = NA_real_,
  sig      = ""
)

# bind it on to tuk_df
tuk_df <- bind_rows(tuk_df, diag_df) %>% 
  # (optional) if you want to keep only one row per pair in case any duplicates:
  distinct(cluster1, cluster2, .keep_all = TRUE)


# 6) Mirror to get a full square matrix
tuk_full <- bind_rows(
  tuk_df,
  tuk_df %>% transmute(
    cluster1 = cluster2,
    cluster2 = cluster1,
  )
)

# 7) Define your desired ordering
cluster_order <- c(
  "Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3",
  "Epithelial cells", "Muscle cells", "Oenocyte",
  "Hemocytes", "Female reproductive system",
  "Unknown 1", "Sensory neuron", "Unknown 2"
)

# 8) Re-factor both axes
tuk_full <- tuk_full %>%
  mutate(
    cluster1 = factor(cluster1, levels = cluster_order),
    cluster2 = factor(cluster2, levels = cluster_order)
  )

tuk_full <- tuk_full %>%
  filter(!is.na(estimate))

# 9) Plot with reference (cluster1) on y, comparison (cluster2) on x
## Pairwise comparison of NV infection probability\nin the fat body cell types (virgin)
tuk_full <- tuk_full %>%
  filter(
    cluster1 %in% c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3"),
    cluster2 %in% c("Fatbody cells 1", "Fatbody cells 2", "Fatbody cells 3")
  )

ggplot(tuk_full, aes(x = cluster2, y = cluster1, fill = estimate)) +
  geom_tile() +
  geom_text(aes(label = sig), vjust = 0.5, hjust = 0.5) +
  scale_fill_gradient2(
    midpoint = 0,
    low      = "blue",
    mid      = "grey95",
    high     = "red",
    name     = "Log-odds\ndifference"
  ) +
  labs(
    x     = "Compared Cluster",
    y     = "Reference Cluster",
    title = "NV (Virgin)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title  = element_text(face = "bold")
  )




