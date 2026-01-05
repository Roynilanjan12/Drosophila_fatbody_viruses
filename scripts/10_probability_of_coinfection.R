## Supplementary figure S7
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first
########### try out the logistic regression to control for Identity + cluster.ident+type
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
#meta <- filter(meta, type=="Mated")
meta$DAV <- NULL
meta$NV <- NULL


meta <- meta %>%
  mutate(
    # Make DAV_status a factor with “Uninfected” as the reference level
    DAV_status    = factor(DAV_status, levels = c("Uninfected", "Infected")),
    # Make NV_status a factor (we’ll also create a 0/1 version below)
    NV_status     = factor(NV_status,  levels = c("Uninfected", "Infected")),
    # Identity as a factor (reference level will be its first alphabetic level unless you override)
    Identity      = factor(Identity),
    # cluster.ident as a factor
    cluster.ident = factor(cluster.ident),
    type = factor(type)
  )

meta <- meta %>%
  mutate(
    NV01 = ifelse(NV_status == "Infected", 1L, 0L)
  )

# ---- Statistical modeling (GLM) ----
# Fit generalized linear models (binomial/logit) to quantify infection probability differences.
fit_logit <- glm(
  NV01 ~ DAV_status + Identity + cluster.ident+type,
  data   = meta,
  family = binomial(link = "logit")
)

# 6A) Summarize the results
summary(fit_logit)


### predictive plotting based on fit_logit
# ---- Load required R packages ----
# Packages are loaded explicitly to make dependencies clear for reproduction.
library(ggplot2)

## Coefficients from your model
b0    <- -1.79856   # (Intercept)
b_DAV <-  0.39936   # DAV(0/1)

## Create a small "fake" dataset: DAV = 0 and DAV = 1
plot_dat <- data.frame(
  DAV = c(0, 1)
)

## Predicted probability of NV infection for each DAV value
plot_dat$logit <- b0 + b_DAV * plot_dat$DAV
plot_dat$prob  <- plogis(plot_dat$logit)   # logistic transform

plot_dat
#   DAV     logit      prob
# 1   0 -1.79856 0.1420264  (~14% when DAV=0)
# 2   1 -1.39920 0.1979431  (~20% when DAV=1)

## Plot
# ---- Visualization ----
# Figures are generated with ggplot2; styling is kept consistent with manuscript themes.
ggplot(plot_dat, aes(x = factor(DAV), y = prob)) +
  geom_point(size = 4) +
  geom_line(aes(group = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "DAV infection (0 = uninfected, 1 = infected)",
    y = "Predicted probability of NV infection",
    title = "Effect of DAV on NV infection (from logistic regression)"
  ) +
  theme_classic()



