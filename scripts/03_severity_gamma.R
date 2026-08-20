# ==============================================================================
# Script: 03_severity_gamma.R
# ==============================================================================

# 1. Setup & Load Data ------------------------------------------------------
library(tidyverse)
library(interactions)
rm(list = ls()) # Clear completely the global environment for reproducibility

# Load the subset strictly for the Severity model (only Claims > 0)
insurance_severity <- readRDS("data/insurance_severity.rds")

# 2. Exploratory Data Analysis (EDA) ----------------------------------------
# Checking the influence of covariates on average claim size
fig_km <- glm(average_claim_size ~ Kilometres, family = Gamma(link = "inverse"), data = insurance_severity)

# Plotting empirical log means
q1 <- cat_plot(fig_km, pred = Kilometres, geom = "line", outcome.scale = "link", y.label = "Empirical Log Means") + 
  theme_light() +
  labs(title = "Effect of Kilometres on Average Claim Size")
print(q1)

# 3. Fit Gamma Regression Model ---------------------------------------------
# Modeling severity with a log-link and using Claims as weights
gamma_model <- glm(formula = average_claim_size ~ Kilometres + Zone + Bonus + Make, 
                   data = insurance_severity, 
                   family = Gamma(link = "log"), 
                   weights = Claims)

summary(gamma_model)
# Note: Dispersion parameter is around 2.95 -> may indicate overdispersion.

# 4. Residual Deviance Test -------------------------------------------------
# Test at alpha = 0.05
alpha <- 0.05
deviance <- gamma_model$deviance
df <- gamma_model$df.residual

# Extract the exact dispersion parameter from the model summary
dispersion_param <- summary(gamma_model)$dispersion

# Reject H0 iff scaled deviance > qchisq(1 - alpha, df)
reject_H0 <- (deviance / dispersion_param) > qchisq(1 - alpha, df = df)
print(paste("Reject H0 (Model is misspecified):", reject_H0)) # FALSE

# Conclusion: We do not reject the Hypothesis H_0 that the model is correctly specified.

# 5. Save Final Severity Model ----------------------------------------------
saveRDS(gamma_model, "output/model_severity_gamma.rds")
