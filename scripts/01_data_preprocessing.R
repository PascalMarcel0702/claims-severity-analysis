# ==============================================================================
# Script: 01_data_preprocessing_and_eda.R
# ==============================================================================
# Clear completely the global environment for reproducibility
rm(list = ls())
# 1. Load Packages ----------------------------------------------------------
library(tidyverse)
library(interactions)

# 2. Load Data --------------------------------------------------------------
insurance <- read.table(file = "data/insurance.txt", header = TRUE, sep = "")

# 3. Initial Formatting -----------------------------------------------------
# Define factor variables based on data scale (Kilometres: ordinal, Zone/Make: nominal)
insurance$Kilometres <- as.factor(insurance$Kilometres)
insurance$Zone <- as.factor(insurance$Zone)
insurance$Make <- as.factor(insurance$Make)

# 4. Exploratory Data Analysis (EDA) ----------------------------------------
ti <- insurance$Insured
f1 <- glm(data = insurance, family = poisson(link = "log"), formula = Claims ~ Kilometres + offset(log(Insured)))

# Plot empirical log-means to check for non-linear effects
p1 <- cat_plot(f1, pred = Kilometres, data = insurance, outcome.scale = "response", 
               geom = "line", y.label = "Log-Mean", set.offset = ti, 
               plot.points = FALSE) + 
  theme_light() +
  labs(title = "Effect of Kilometres on Claim Frequency",
       subtitle = "Empirical Log-Means before Category Merging")
print(p1)

# Export plot for documentation
dir.create("output", showWarnings = FALSE)
ggsave("output/eda_kilometres_before_merge.png", plot = p1, width = 6, height = 4)

# Conclusion from EDA: Kilometres shows a non-linear effect. 
# Categories 2, 3, and 4 exhibit similar log-means. 
# Merging them linearizes the relationship.

# 5. Feature Engineering (Adjusting based on EDA) ---------------------------

# Adjusting Kilometres
levels(insurance$Kilometres)[c(2, 3, 4)] = 2
levels(insurance$Kilometres)

# Adjusting Make (shifting 8 to 7, and 9 to 8)
Make.new = insurance$Make
Make.new[insurance$Make == 8] = 7
Make.new[insurance$Make == 9] = 8
insurance$Make = factor(Make.new)

# Calculate average claim size for Gamma Regression (Severity)
insurance <- insurance |> mutate(average_claim_size = Payment / Claims)

# Create subset strictly for the Severity (Gamma) model (Claims > 0)
insurance_severity <- insurance |> filter(Claims > 0)

# 6. Save Processed Data & Cleanup ------------------------------------------
# Save two distinct datasets for downstream modeling
saveRDS(insurance, "data/insurance_frequency.rds")
saveRDS(insurance_severity, "data/insurance_severity.rds")

# Clean up temporary variables
rm(ti, f1, p1, Make.new)

# Quick check to ensure everything worked
summary(insurance)