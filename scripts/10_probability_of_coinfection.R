## Supplementary figure S9
###running 2_cell_type_cluster_NV_DAV_tropism_with_infection_percentage_in_cell_types.R script first
########### try out the logistic regression to control for Identity + cluster.ident+type
meta<- fatbody_v.harmony@meta.data
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
meta <- filter(meta, Identity=="M_2")

meta <- meta %>%
  mutate(
    NV01 = ifelse(NV_status == "Infected", 1L, 0L)
  )

fit_logit <- glm(
  NV01 ~ DAV_status + cluster.ident,
  data   = meta,
  family = binomial(link = "logit")
)

# 6A) Summarize the results
summary(fit_logit)


### predictive plotting based on fit_logit
library(ggplot2)

## Coefficients from your model
b0    <- -1.78348   # (Intercept)
b_DAV <-  0.41171   # DAV(0/1)

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
ggplot(plot_dat, aes(x = factor(DAV), y = prob)) +
  geom_point(size = 4) +
  geom_line(aes(group = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "DAV infection \n (0 = uninfected, 1 = infected)",
    y = "Predicted probability of NV infection"
  ) +
  theme_classic()



# ... [Assume previous data fetching and factor creation is done as above] ...

# 1. Create the binary 0/1 variable for DAV
meta <- meta %>%
  mutate(
    DAV01 = ifelse(DAV_status == "Infected", 1L, 0L)
  )

# 2. Run the Logistic Regression
# Response: DAV01 (0 or 1)
# Predictors: NV_status + Identity + cluster.ident + type
fit_logit_DAV <- glm(
  DAV01 ~ NV_status + cluster.ident ,
  data   = meta,
  family = binomial(link = "logit")
)

# 3. Summarize the results
summary(fit_logit_DAV)

library(ggplot2)

## 1. Define Coefficients from your new model
b0   <- 0.508307  # (Intercept)
b_NV <- 0.411706  # Coefficient for NV_statusInfected

## 2. Create the dataset with NV = 0 and NV = 1
plot_dat <- data.frame(
  NV = c(0, 1)
)

## 3. Calculate predicted probability of DAV infection
# Formula: logit = Intercept + (Effect of NV * NV status)
plot_dat$logit <- b0 + b_NV * plot_dat$NV
plot_dat$prob  <- plogis(plot_dat$logit)   # logistic transform

# View the calculated probabilities
plot_dat
#   NV    logit      prob
# 1  0  0.508307  0.6244105  (~62.4% chance of DAV when NV=0)
# 2  1  0.920013  0.7150436  (~71.5% chance of DAV when NV=1)

## 4. Plot
ggplot(plot_dat, aes(x = factor(NV), y = prob)) +
  geom_point(size = 4) +
  geom_line(aes(group = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x = "NV infection \n (0 = uninfected, 1 = infected)",
    y = "Predicted probability of DAV infection",
  ) +
  theme_classic()


# Create a dummy dataset for prediction
# We must include 'cluster.ident' because it's in the model. 
# We pick the most common cluster or the reference one.
ref_cluster <- levels(meta$cluster.ident)[1] 

newdata <- data.frame(
  DAV_status = factor(c("Uninfected", "Infected"), levels = c("Uninfected", "Infected")),
  cluster.ident = factor(ref_cluster, levels = levels(meta$cluster.ident))
)

# Predict returns the fit AND the Standard Error (se.fit)
preds <- predict(fit_logit_NV, newdata = newdata, type = "link", se.fit = TRUE)

# Calculate Confidence Intervals on the Logit scale, then convert to Probability
newdata$fit_logit <- preds$fit
newdata$lower_logit <- preds$fit - (1.96 * preds$se.fit)
newdata$upper_logit <- preds$fit + (1.96 * preds$se.fit)

# Convert all to Probability (Inverse Logit)
newdata$prob <- plogis(newdata$fit_logit)
newdata$lower_prob <- plogis(newdata$lower_logit)
newdata$upper_prob <- plogis(newdata$upper_logit)

# Plot
ggplot(newdata, aes(x = DAV_status, y = prob)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = lower_prob, ymax = upper_prob), width = 0.2) +
  geom_line(aes(group = 1)) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "DAV Status", y = "Probability of NV") +
  theme_classic()



library(emmeans)
library(ggplot2)

# --- MODEL 1: NV dependent on DAV ---
fit_logit_NV <- glm(
  NV01 ~ DAV_status + cluster.ident,
  data = meta,
  family = binomial(link = "logit")
)

# Calculate expected probabilities (averaged over clusters)
# type="response" converts logit -> probability automatically
emm_NV <- emmeans(fit_logit_NV, ~ DAV_status, type = "response")
emm_NV_df <- as.data.frame(emm_NV)

# Plot with Confidence Intervals
p1 <- ggplot(emm_NV_df, aes(x = DAV_status, y = prob)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) + # Add Error Bars
  geom_line(aes(group = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Effect of DAV on NV",
    x = "DAV Status",
    y = "Predicted Probability of NV Infection"
  ) +
  theme_classic()

print(p1)

# --- MODEL 2: DAV dependent on NV ---

# 1. Fit the model
fit_logit_DAV <- glm(
  DAV01 ~ NV_status + cluster.ident,
  data = meta,
  family = binomial(link = "logit")
)

# 2. Calculate expected probabilities (averaged over clusters)
# type="response" converts logit -> probability automatically
emm_DAV <- emmeans(fit_logit_DAV, ~ NV_status, type = "response")
emm_DAV_df <- as.data.frame(emm_DAV)

# 3. Plot with Confidence Intervals
p2 <- ggplot(emm_DAV_df, aes(x = NV_status, y = prob)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) + # Add Error Bars
  geom_line(aes(group = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Effect of NV on DAV",
    x = "NV Status",
    y = "Predicted Probability of DAV Infection"
  ) +
  theme_classic()

print(p2)







