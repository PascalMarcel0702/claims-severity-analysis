# ==============================================================================
# Script: 02_frequency_poisson.R
# ==============================================================================

# 1. Setup & Load Data ------------------------------------------------------
library(tidyverse)
rm(list = ls()) # Clear completely the global environment for reproducibility

insurance <- readRDS("data/insurance_frequency.rds")

# 2. Stepwise Covariate Selection (Main Effects) ----------------------------
# Select the Poisson model for the data for which the covariates are ordered by significance (AIC).
# Response: Claims | Covariates: Kilometres, Zone, Bonus, Make | Offset: log(Insured)

model_null <- glm(formula = Claims ~ 1 + offset(log(Insured)), data = insurance, family = poisson(link = "log"))
model_full <- glm(formula = Claims ~ Kilometres + Zone + Bonus + Make + offset(log(Insured)), data = insurance, family = poisson(link = "log"))

# Stepwise forward selection
model_main <- step(model_null, scope = list(lower = model_null, upper = model_full), direction = "forward")
summary(model_main)
# Answer: Bonus + Zone + Kilometres + Make are all significant main effects.

# 3. Interaction Effects ----------------------------------------------------
# Check for significant interaction effects at alpha = 0.05
model_inter <- step(model_main, ~.^2, direction = "forward")
summary(model_inter)

# 4. Partial Deviance Test --------------------------------------------------
# Perform a partial deviance test to compare model_main and model_inter
df_main <- model_main$df.residual
dev_main <- model_main$deviance

df_inter <- model_inter$df.residual
dev_inter <- model_inter$deviance

delta_dev <- dev_main - dev_inter
delta_df <- df_main - df_inter

# Reject H0: \beta_2 = 0 at level alpha = 0.05 iff delta_dev > qchisq(0.95, delta_df)
delta_dev > qchisq(0.95, df = delta_df) # TRUE
# Interpretation: We reject H0. That indicates that the added covariates/interactions 
# are significant. Therefore we decide to use the interaction model (model_inter).

# 5. Model Diagnostics & Overdispersion -------------------------------------
# Residual deviance test for model_inter
dev_inter > qchisq(0.95, df = df_inter) # TRUE
# We reject H0 that the model assumptions of the specified GLM are satisfied. 
# A reason for this might be overdispersion.

# Rule of thumb for overdispersion: deviance > df
dev_inter > df_inter # TRUE -> indication of overdispersion.

# Compute and plot deviance residuals to check variance structure
res_dev <- residuals(model_inter, type = "deviance")
n <- length(res_dev)

p_dev <- ggplot(data.frame(y = res_dev, x = 1:n), aes(x = x, y = y)) +
  geom_point() +
  theme_light() +
  labs(title = "Deviance Residuals of Frequency Model",
       x = "Observation Number",
       y = "Deviance Residuals")
print(p_dev)

# Export plot for documentation
ggsave("output/poisson_deviance_residuals.png", plot = p_dev, width = 6, height = 4)

# Interpretation: The scatter plot of the deviance residuals shows a random pattern 
# around their expected zero mean which indicates that the model reflects the variability 
# of the data well.

# 6. Save Final Frequency Model ---------------------------------------------
saveRDS(model_inter, "output/model_frequency_poisson.rds")