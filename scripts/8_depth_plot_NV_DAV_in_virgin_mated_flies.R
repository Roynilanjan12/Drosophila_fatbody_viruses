##Supplementary figure S2
# ---- Load required R packages ----
# Packages are loaded explicitly to make dependencies clear for reproduction.
library(tidyverse)
library(reshape2)
library(Seurat)

# =========================
# 1) Read all depth files
# =========================

# Mated Uninfected:   SRS8175458, SRS8175463
# Mated Infected:     SRS8175460, SRS8175464
# Virgin Uninfected:  SRS8175457, SRS8175461
# Virgin Infected:    SRS8175459, SRS8175462

## Virgin Uninfected
Virgin_Uninfected_1 <- read.table(
  "/data/VU1.Depth",
  header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE
)
colnames(Virgin_Uninfected_1) <- c("Chr", "BP", "Depth")
Virgin_Uninfected_1$Sample <- "VU1"

Virgin_Uninfected_2 <- read.table(
  "/data/VU2.Depth",
  header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE
)
colnames(Virgin_Uninfected_2) <- c("Chr", "BP", "Depth")
Virgin_Uninfected_2$Sample <- "VU2"


## Virgin Infected
Virgin_Infected_1 <- read.table(
  "/data/VI1.Depth",
  header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE
)
colnames(Virgin_Infected_1) <- c("Chr", "BP", "Depth")
Virgin_Infected_1$Sample <- "VI1"

Virgin_Infected_2 <- read.table(
  "/data/VI2.Depth",
  header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE
)
colnames(Virgin_Infected_2) <- c("Chr", "BP", "Depth")
Virgin_Infected_2$Sample <- "VI2"


## Mated Uninfected
Mated_Uninfected_1 <- read.table(
  "/data/MU1.Depth",
  header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE
)
colnames(Mated_Uninfected_1) <- c("Chr", "BP", "Depth")
Mated_Uninfected_1$Sample <- "MU1"

Mated_Uninfected_2 <- read.table(
  "/data/MU2.Depth",
  header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE
)
colnames(Mated_Uninfected_2) <- c("Chr", "BP", "Depth")
Mated_Uninfected_2$Sample <- "MU2"


## Mated Infected
Mated_Infected_1 <- read.table(
  "/data/MI1.Depth",
  header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE
)
colnames(Mated_Infected_1) <- c("Chr", "BP", "Depth")
Mated_Infected_1$Sample <- "MI1"

Mated_Infected_2 <- read.table(
  "/data/MI2.Depth",
  header = FALSE, sep = "\t", na.strings = "NA", dec = ".", strip.white = TRUE
)
colnames(Mated_Infected_2) <- c("Chr", "BP", "Depth")
Mated_Infected_2$Sample <- "MI2"


# =========================
# 2) Combine & clean
# =========================

all_sample_depth <- rbind(
  Virgin_Uninfected_1,
  Virgin_Uninfected_2,
  Virgin_Infected_1,
  Virgin_Infected_2,
  Mated_Uninfected_1,
  Mated_Uninfected_2,
  Mated_Infected_1,
  Mated_Infected_2
)

rm(
  Virgin_Uninfected_1,
  Virgin_Uninfected_2,
  Virgin_Infected_1,
  Virgin_Infected_2,
  Mated_Uninfected_1,
  Mated_Uninfected_2,
  Mated_Infected_1,
  Mated_Infected_2
)
gc()   # optional


# =========================
# 3) Recode virus names
# =========================

library(dplyr)

all_sample_depth <- all_sample_depth %>%
  mutate(
    Chr = case_when(
      Chr == "JX220408.1_Nora_virus_isolate_FR1" ~ "NV",
      Chr == "KP969946.1_Drosophila_A_virus_isolate_LJ35" ~ "DAV",
      TRUE ~ Chr
    )
  )

virus_depth <- all_sample_depth %>%
  filter(Chr %in% c("NV", "DAV")) %>%
  mutate(
    Virus  = factor(Chr, levels = c("NV", "DAV")),
    Sample = factor(Sample, levels = c("VU1", "VU2", "VI1", "VI2",
                                       "MU1", "MU2", "MI1", "MI2"))
  )


# =========================
# 4) Plot: coverage per sample × virus (no log transform)
# =========================


p <- ggplot(virus_depth, aes(x = BP, y = Depth, fill = Virus)) +
  geom_area() +
  scale_fill_manual(
    values = c(
      "NV"  = "blue",
      "DAV" = "red"
    )
  ) +
  facet_grid(Sample ~ Virus, scales = "free") +
  xlab("Virus genome position (BP)") +
  ylab("Depth") +
  theme(text = element_text(size = 16))

# Optional to remove legend:
# p <- p + NoLegend()

p
