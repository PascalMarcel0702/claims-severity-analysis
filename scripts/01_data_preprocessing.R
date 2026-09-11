# ==============================================================================
# Script: 01_data_preprocessing_and_eda.R
# ==============================================================================
# Clear completely the global environment for reproducibility
# Clear environment for reproducibility
rm(list = ls())

# 1. Setup ---------------------------------------------------------------------
library(tidyverse)
library(patchwork)   # For professional side-by-side EDA plots
library(mgcv)        # For potential GAM checks (if continuous variables existed)

# Create output directories
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# 2. Data Import & Inspection --------------------------------------------------
insurance <- read.table(file = "data/insurance.txt", header = TRUE, sep = "")

str(insurance) # 2182 obs, 7 variables
sum(insurance$Insured) # 2383170 years
sum(insurance$Claims) # 113171 claims
# Check for missing values
colSums(is.na(insurance)) # no missing values

# Preparation ----------------------------------------
# Covariate & Response Scales:
# [Responses]
# - Claims (Frequency): Count variable (Poisson GLM with canonical link (log)). 
# - average_claim_size (Severity): Continuous, strictly positive (Gamma GLM).
# [Exposure / Offset]
# - Insured: Continuous (policy-years). Modeled as log(Insured) due to underlying link function. 
# [Predictors]
# - Kilometres: Ordinal categorical (1 to 5, distance bands).
# - Bonus: Ordinal categorical (1 to 7, years claim-free).
# - Zone: Nominal categorical (1 to 7, geographic regions).
# - Make: Nominal categorical (1 to 9, car models).

# Define average claim size (for gamma regression) as 
insurance <- insurance %>%
  mutate(
    average_claim_size = ifelse(Claims > 0, Payment / Claims, NA)
  )

# Factorization to avoid linear relationship along the categories
insurance$Zone       <- as.factor(insurance$Zone)
insurance$Make       <- as.factor(insurance$Make)

# Factorize ordinal covariates
# Reason: The data consists of policies corresponding to 2383170 years and 113171 claims. As result of the large data aggregation, the approximate confidence limits are very narrow.Thus, a priori the disadvantage of higher complexity (in contrast to metric transformation) due to estimating more parameters (and this higher BIC penalty) is accepted to develop a more precise model

insurance$Kilometres <- as.factor(insurance$Kilometres)
insurance$Bonus      <- as.factor(insurance$Bonus)

levels(insurance$Kilometres) # 1 - 5
levels(insurance$Zone) # 1 - 7
levels(insurance$Bonus) # 1- 7
levels(insurance$Make) # 1 - 9

summary(insurance)
view(insurance)

# Create subset strictly for servity (Gamma) model
insurance_severity <- insurance %>% filter(Claims > 0)

# Pre-Checks: Design Matrix Rank for Algorithmic Stability
X_candidate_freq <- model.matrix(~ Kilometres + Zone + Bonus + Make, data = insurance)
X_candidate_sev <- model.matrix(~ Kilometres + Zone + Bonus + Make, data = insurance_severity)

c(qr(X_candidate_freq)$rank == ncol(X_candidate_freq), qr(X_candidate_sev)$rank == ncol(X_candidate_sev)) # TRUE TRUE
# Interpretation: Both design matrix have full rank

# Explatory Data Analysis

# Helper function pre-check of sparse categories
check_sparsity <- function(
    df,
    group_var,
    var_name,
    min_claims = 100, # portfolio contains 113.000 Claims
    min_exposure_share = 0.01
) {
  
  df %>%
    group_by({{ group_var }}) %>%
    summarise(
      Rows = n(),
      Total_Exposure = sum(Insured, na.rm = TRUE),
      Total_Claims = sum(Claims, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Exposure_Share = Total_Exposure / sum(Total_Exposure),
      Claim_Share = Total_Claims / sum(Total_Claims),
      Claim_Frequency = Total_Claims / Total_Exposure,
      Claim_to_Exposure_Ratio = Claim_Share / Exposure_Share,
      
      # Sparsity warning
      Sparsity_Warning = case_when(
        Total_Exposure == 0 ~ 
          "Critical: 0 exposure (undefined frequency)",
        
        Total_Claims == 0 ~ 
          "Critical: 0 claims (empirical log-rate undefined)",
        
        Total_Claims < min_claims ~ 
          paste0("Warning: < ", min_claims, " claims (high variance)"),
        
        Exposure_Share < min_exposure_share ~ 
          paste0("Warning: < ", 100 * min_exposure_share, "% of total exposure"),
        
        TRUE ~ "OK"
      )
    ) %>%
    # Format for clean README output
    mutate(
      Variable = var_name,
      Category = as.character({{ group_var }}),
      Exposure_Share = sprintf("%.2f %%", 100 * Exposure_Share),
      Claim_Share = sprintf("%.2f %%", 100 * Claim_Share),
      Claim_Frequency = round(Claim_Frequency, 4),
      Claim_to_Exposure_Ratio = round(Claim_to_Exposure_Ratio, 2)
    ) %>%
    # Final selection: Claim_Frequency is now included
    select(
      Variable,
      Category,
      Rows,
      Total_Exposure,
      Total_Claims,
      Exposure_Share,
      Claim_Share,
      Claim_Frequency, 
      Claim_to_Exposure_Ratio,
      Sparsity_Warning
    )
}

# Run sparsity check for all covariates
sparsity_km    <- check_sparsity(insurance, Kilometres, "Kilometres")
sparsity_bonus <- check_sparsity(insurance, Bonus, "Bonus")
sparsity_zone  <- check_sparsity(insurance, Zone, "Zone")
sparsity_make  <- check_sparsity(insurance, Make, "Make")

# Combine all into one master table for the README
sparsity_master <- bind_rows(sparsity_km, sparsity_bonus, sparsity_zone, sparsity_make)
view(sparsity_master)
# Export as CSV (can be easily converted to a Markdown table for the README)
write.csv(sparsity_master, "output/tables/eda_sparsity_check.csv", row.names = FALSE, quote = FALSE, na = "")

# Interpretation of marginal sparsity analysis:
# Only Zone 7 is flagged consisting 0.80 % of total exposure and 0.55 % of total claims.This category contains data corresponding to 19,000 policy years, which is < 1 % of the total policy years aggregated over the whole data set. This record is therefore flagged in the sparsity check.
# All in all: Less than 1% of the total policy years in the portfolio belong to drives from Zone 7 (Gotland)

# Conclusion: Every category across all covariates contains large data aggregation, i.e., no category requires merging in face of variance reduction.


# Plot of Empirical log rates with CIs to identify relationship between categorical covariates and response and investigate which categories can be merged (similar log-rates and overlapping CIs)

# Helper function for empirical log rates regarding poisson model
plot_poisson_eda <- function(df, grouping_var, var_label, alpha = 0.05) {
  df %>%
    group_by({{ grouping_var }}) %>%
    summarise(
      Y_j = sum(Claims),
      t_j = sum(Insured),
      .groups = "drop"
    ) %>%
    mutate(
      # Calculate approximate 100 * (1- \alpha) % CI based on delta method
      rate = Y_j / t_j,
      emp_log_rate = log(rate),
      se_rate = 1 / sqrt(Y_j),
      ci_lower = emp_log_rate - qnorm(1 - alpha / 2) * se_rate,
      ci_upper = emp_log_rate + qnorm(1 - alpha / 2) * se_rate
    ) %>%
    ggplot(aes(x = {{ grouping_var }}, y = emp_log_rate, group = 1)) +
    geom_line(color = "#2c3e50", linetype = "dashed", alpha = 0.6) +
    geom_point(color = "#2c3e50", size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "#2c3e50") +
    theme_light() +
    labs(x = var_label, y = "Empirical Log-Rate") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(breaks = seq(-4, -1, by = 0.2))
}

# Helper function for Poisson EDA Tables
generate_poisson_table <- function(df, grouping_var, var_name, alpha = 0.05) {
  df %>%
    group_by({{ grouping_var }}) %>%
    summarise(
      Y_j = sum(Claims),
      t_j = sum(Insured),
      .groups = "drop"
    ) %>%
    mutate(
      rate = Y_j / t_j,
      se_rate = 1 /sqrt( Y_j),
      
      `Empir. Log-Rate` = round(log(rate), 3),
      `L` = round(mp_log_rate - qnorm(1 - alpha / 2) * se_rate, 3),
      `U` = round(emp_log_rate + qnorm(1 - alpha / 2) * se_rate, 3)
    ) %>%
    select({{ grouping_var }}, `Empir. Log-Rate`, Y_j, t_j, L, U) %>%
    mutate(Variable = var_name, .before = 1)
}

# Generate Poisson Plots
p_freq_km    <- plot_poisson_eda(insurance, Kilometres, "Kilometres (Distance Band)")
p_freq_bonus <- plot_poisson_eda(insurance, Bonus, "Bonus (Claim-Free Years)")
p_freq_zone  <- plot_poisson_eda(insurance, Zone, "Zone (Geographical Region)")
p_freq_make  <- plot_poisson_eda(insurance, Make, "Make (Car Model)")


ggsave("output/figures/eda_freq_kilometres.png", plot = p_freq_km, width = 6, height = 4, dpi = 300)
ggsave("output/figures/eda_freq_bonus.png", plot = p_freq_bonus, width = 6, height = 4, dpi = 300)
ggsave("output/figures/eda_freq_zone.png", plot = p_freq_zone, width = 6, height = 4, dpi = 300)
ggsave("output/figures/eda_freq_make.png", plot = p_freq_make, width = 6, height = 4, dpi = 300)


# Generate Poisson Master Table
tab_freq_km    <- generate_poisson_table(insurance, Kilometres, "Kilometres")
tab_freq_bonus <- generate_poisson_table(insurance, Bonus, "Bonus")
tab_freq_zone  <- generate_poisson_table(insurance, Zone, "Zone")
tab_freq_make  <- generate_poisson_table(insurance, Make, "Make")

eda_freq_master <- bind_rows(tab_freq_km, tab_freq_bonus, tab_freq_zone, tab_freq_make)
write.csv(eda_freq_master, "output/tables/eda_frequency_summary.csv", row.names = FALSE, na = "")

# Interpretation of Empirical Log-Rates by Log-Rate-Delta and CI separation
# Note: Approximate confidence limits are very narrow which is the result of the large data aggregation.

# 1. Primary Risk Driver (Highest Delta, no overlapping CIs):
# - 'Bonus': Marginal effect size Delta = 1.252 (Extremes: Level 1: -2.129 vs. Level 7: -3.381).The plot shows a strictly monotonic downward trend with steep decline in categories 1-2 and 6-7. Higher claim-free years correspond to lower empirical claim frequencies. None of the 7 levels exhibit CI overlap, indicating distinct risk profiles across the categories. No merging is required.

# 2. Secondary Risk Drivers (Moderate Delta, Nominal Structure):
# - 'Zone'(nominal): Marginal effect size Delta = 0.782 (Extremes: Zone 1: -2.645 vs. Zone 7: -3.427). Urban areas (Zone 1) show highest risk, while rural/island areas (Zone 7) show the lowest, but no linear trend in between. All CIs are strictly separated. Category 7 has largest CI width of 0.158 caused by marginal data sparsity as discussed before.
# - 'Make'(nominal): Marginal effect size Delta = 0.618 (Extremes: Make 5: -2.854 vs. Make 4: -3.472).Car Model 5 shows highest empirical claim frequency, while model 4 shows the lowest. Category 7 and 8 show similar risk profile with overlapping CIs.

# 3. Weak Risk Driver (Lowest Delta, neighboring CIs):
# - 'Kilometres': Marginal effect size Delta = 0.431 (Extremes: Level 5: -2.760 vs. Level 1: -3.191). The plot exhibits a strictly monotonic upward trend (more driving = higher risk) with plateau in categories 2, 3 and 4, which show similar empirical claim frequencies with neighboring CIs (2_U = -3.007 = 3_L; 3_U = 2.982 \sim 4_L = -2.975 )

# Note : Merging Zone
# Zones 4, 6, and 7 represent rural areas. Since the empirical log-rate plot reveals the structural distinct risk profile of these bins. In particular,the CIs of Zone 4 (rural south) and Zone 7 (Gotland (South)) are disjoint, reflecting the differences in mainland transit vs. isolated island dynamics.

# Data merging based on empirical log rates

insurance <- insurance %>%
  mutate(
    # Kilometres 2, 3, 4  plateau
    Kilometres_merged = ifelse(Kilometres %in% c("2", "3", "4"), "2-4", as.character(Kilometres)),
    Kilometres_merged = factor(Kilometres_merged, levels = c("1", "2-4", "5")),
    
    # Make 7, 8 similar log rates and CIs
    Make_merged = ifelse(Make %in% c("7", "8"), "7-8", as.character(Make)),
    Make_merged = factor(Make_merged, levels = c("1", "2", "3", "4", "5", "6", "7-8", "9"))
  )

# Plot of empirical Log-Rates:
p_freq_km_merged    <- plot_poisson_eda(insurance, Kilometres_merged, "Kilometres (Distance Band)")
# Interpretation: Unchanged Delta = 0.431, plot shows strictly linear upward trend with strictly disjoint CIs.
p_freq_make_merged  <- plot_poisson_eda(insurance, Make_merged, "Make (Car Model)")
# Interpretation: Unchanged Delta = 0.618 and extremes, no homogeneous pattern and no categories with similar empirical logit rate nor with neighboring CIs.


# Interaction Plots
# Helper function for Poisson Interaction plots (Empirical Log-Rate)
plot_poisson_interaction <- function(df, x_var, group_var,alpha = 0.05, x_label, legend_label) {
  
  # Grouping and calculating metrics
  agg_data <- df %>%
    group_by({{ x_var }}, {{ group_var }}) %>%
    summarise(
      Y_j = sum(Claims),
      t_j = sum(Insured),
      .groups = "drop"
    ) %>%
    # Avoid errors (log(0) or div/0)
    filter(t_j > 0, Y_j > 0) %>%
    mutate(
      # Calculate CI and log rate
      rate = Y_j / t_j,
      emp_log_rate = log(rate),
      se_rate = 1 / sqrt( Y_j ),
      ci_lower = emp_log_rate - qnorm(1 - alpha / 2) * se_rate,
      ci_upper = emp_log_rate + qnorm(1 - alpha / 2) * se_rate
    )
  
  # Doge lines for better overview
  pd <- position_dodge(width = 0.3)
  
  # Generate Plot
  ggplot(agg_data, aes(x = {{ x_var }}, y = emp_log_rate, color = {{ group_var }}, group = {{ group_var }})) +
    geom_point(size = 3, position = pd) +
    geom_line(linewidth = 1, linetype = "dashed", alpha = 0.7, position = pd) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = pd) +
    theme_light() +
    labs(
      x = x_label,
      y = "Empirical Log-Rate",
      color = legend_label
    ) +
    theme(
      legend.position = "right",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}


# 1) Interaction Bonus and Zone: In crowded cities, such as Stockholm, the traffic density is higher than in rural areas, such as Gotland. 
# Q: Does the geographic region correlates with the bonus level, i,e, the time since the last claim was reported?

p_inter_Zone_Bonus <- plot_poisson_interaction(insurance, Zone, Bonus, "Zone", "Bonus")
# Interpretation: Broadly parallel zig-zag trend w.r.t. all Bonus classes in Categories 1-6 of Zone. Non-linearities in Zone 1-2 across all all Bonus bins visible, caused by varying slopes of the linear interpolation leading to changing risk gaps accompanied by narrow and approximately distinct CIs. From Zone 3 onwards, the CIs are strongly overlapping / absorbing and increase in its widths, reaching their maximum in Zone 7 caused by data sparsity, masking visual interaction effects. 
# All in all: Visual indication of local interaction effect between Bonus and Zone bin 1-2 (urban zones).


# Zoom in:
# Urban areas
insurance_urban <- insurance %>% filter(Zone %in% c("1", "2")) %>% mutate(droplevels(Zone))
p_inter_Zone_urban_Bonus <- plot_poisson_interaction(insurance_urban, Zone, Bonus, "Zone urban", "Bonus")
# Interpretation: The linear interpolation looks nearly parallel. The CIs of Bonus 5,6 in Zone 1 strongly overlap and the CIs of Bonus bins 4,5 and 5,6 overlap in Zone 2, all other bins have narrow CI witdhs and are distinct.

# Conclusion (Zone vs. Bonus): The zoomed visualization does not support the hypothesis of a local interaction. There is no significant visual indication of presence of a local or global interaction effect between Bonus and Zone.


# 2) Interaction Bonus and Kilometres: Driving large distances can be exhausting and intuitively leads to higher risk of having accidents. Therefore, does the average distance traveled by car influences the amount of time, since the last claim was reported?
# Q: Does the distance traveled correlates with the amount of claim-free years?
p_inter_Kilometres_Bonus <- plot_poisson_interaction(insurance, Kilometres, Bonus, "Kilometres", "Bonus")
p_inter_Kilometres_merged_Bonus <- plot_poisson_interaction(insurance, Kilometres_merged, Bonus, "Kilometres", "Bonus ")

# Comparison of both plots:
comparison_inter_Kilometres_Bonus <- p_inter_Kilometres_Bonus + p_inter_Kilometres_merged_Bonus + theme(axis.title.y = element_blank()) +
  plot_layout(guides = "collect")

# Interpretation: Both plots show broadly parallel trend across most categories, whereby the parallelism in the merged plot is more clear. In both plots, the narrow CIs in groups 4-6 of Bonus are neighboring across all categories, but are not absorbing, except local instances in the mid-distance range. In the unmerged plot, the linearly interpolated lines reveal non-parallelism caused by difference in slope ( e.g., between Bonus class 1 and 6 ). In contrast, the slopes in the merged plot behave approximately similarly.
# All in all: After taking narrow approximate CIs into account, non-linearities are notable and visually indicate the presence of a global interaction effect regarding the umerged version of Kilometres. On the other hand, the merged plot does not visually indicate a significant interaction effect by smoothing out the irregularities.


# 3) Interaction Zone and Kilometres: In rural areas, shopping facilities / workplaces are rare leading to periodically traveling larger distances, while living in crowded regions, shorter distances are traveled, since essential destinations are packed closely together.
# Q: Does the distance traveled correlates with the zone where the insured lives?

p_inter_Kilometres_Zone <- plot_poisson_interaction(insurance, Kilometres, Zone, "Kilometres", "Zone")
p_inter_Kilometres_merged_Zone <- plot_poisson_interaction(insurance, Kilometres_merged, Zone, "Kilometres", "Zone")

# Comparison of both plots:
comparison_inter_Kilometres_Zone <- p_inter_Kilometres_Zone + p_inter_Kilometres_merged_Zone + theme(axis.title.y = element_blank()) +
  plot_layout(guides = "collect")

# Interpretation: Both plots show broadly parallel trend - Zone 3, 5 and 4,6 intersect between Kilometres bins 3-5, accompanied by overlapping CIs. Marginal sparsity of bin 7 of zone leads to combinatorial sparsity across all kilometres bins leading to absorbing CIs regarding Zone 4, 6 in Kilometres bin 1 and 5. In the unmerged plot, this behavior occurs as well in bins 3, 4 of Kilometres.
# ALl in all: No significant visual indication of an strong interaction effect between Kilometres and Zone - differences in slope can be caused by statistical noise.


# 4) Interaction Bonus and Make: Car types that beginners usually drive or fast sport cars cause intuitively more accidents than car types for families, since the latter drive (intuitively) more defensiv than the former.
# Q: Does the Car type (Make) correlates with the amount of claim-free years? 

p_inter_Make_Bonus <- plot_poisson_interaction(insurance, Make, Bonus, "Make", "Bonus")
p_inter_Make_merged_Bonus <- plot_poisson_interaction(insurance, Make_merged, Bonus, "Make", "Bonus")

# Comparison of both plots:
comparison_inter_Make_Bonus <- p_inter_Make_Bonus + p_inter_Make_merged_Bonus + theme(axis.title.y = element_blank()) +
  plot_layout(guides = "collect")

# Interpretation: Both plots show a broadly parallel and volatile zig-zag trend - in particular, Bonus bins 3-6 reveals neighboring / overlapping CIs across all Make bins. Moreover, Bonus categories 2,3 cross in the unmerged plot between Make categories 7,8 accompanied by overlapping CIs. In the merged plot, there is no crossing regarding these categories, the CI widths are weakened, but remain overlapping.In addition, several weak intersections of bonus classes 4,5,6 between make bins and bins 7-9 are visible, but are accompanied by overlapping / absorbing CIs. In particular, the intersection between categories 5,6 of Bonus in the range of categories 7/8, 9 in the merged plot remains and shows simillary overlapping CIs.
# All in all: Despite the existence of slope differences and corssings, no significant visual indication of an strong interaction effect between Make and Bonus, since these local irregularities are accompanied by heavily overlapping / absorbing CIs.


# Table of Empirical Log Rates for Interactions
generate_poisson_interaction_table <- function(df, var1, var2, name_var1, name_var2, alpha = 0.05) {
  df %>%
    group_by({{ var1 }}, {{ var2 }}) %>%
    summarise(
      Y_j = sum(Claims),
      t_j = sum(Insured),
      .groups = "drop"
    ) %>%
    # Avoid errors (log(0) or div/0) caused by empty cells
    filter(t_j > 0, Y_j > 0) %>%
    # Calculate CIs
    mutate(
      rate = Y_j / t_j,
      se_rate = 1 /sqrt( Y_j),
      
      `Empir. Log-Rate` = round(log(rate), 3),
      `L` = round(emp_log_rate - qnorm(1 - alpha / 2) * se_rate, 3),
      `U` = round(emp_log_rate + qnorm(1 - alpha / 2) * se_rate, 3)
    ) %>%
    select({{ var1 }}, {{ var2 }}, `Empir. Log-Rate`, Y_j, t_j, L, U) %>%
    # Adjustment of format
    mutate(t_j = round(t_j, 2)) %>%
    # Rows
    pivot_longer(
      cols = c(`Empir. Log-Rate`, Y_j, t_j, L, U), 
      names_to = "Metric", 
      values_to = "Value"
    ) %>%
    # Columns
    pivot_wider(
      names_from = {{ var2 }}, 
      values_from = Value
    ) %>%
    # Row order and refinement
    mutate(Metric = factor(Metric, levels = c("Empir. Log-Rate", "Y_j", "t_j", "L", "U"))) %>%
    arrange(Metric, {{ var1 }}) %>%
    rename(Category_Var1 = {{ var1 }}) %>%
    mutate(Interaction = paste0(name_var1, " vs ", name_var2), .before = 1)
}


# Generate interaction tables
tab_inter_bonus_zone <- generate_poisson_interaction_table(
  insurance, Bonus, Zone, "Bonus Class", "Zone"
)
tab_inter_bonus_km <- generate_poisson_interaction_table(
  insurance, Bonus, Kilometres, "Bonus Class", "Kilometres"
)
tab_inter_bonus_km_merged <- generate_poisson_interaction_table(
  insurance, Bonus, Kilometres_merged, "Bonus Class", "Kilometres (Merged)"
)
tab_inter_zone_km <- generate_poisson_interaction_table(
  insurance, Zone, Kilometres, "Zone", "Kilometres"
)
tab_inter_zone_km_merged <- generate_poisson_interaction_table(
  insurance, Zone, Kilometres_merged, "Zone", "Kilometres (Merged)"
)
tab_inter_bonus_make <- generate_poisson_interaction_table(
  insurance, Bonus, Make, "Bonus Class", "Make"
)
tab_inter_bonus_make_merged <- generate_poisson_interaction_table(
  insurance, Bonus, Make_merged, "Bonus Class", "Make (Merged)"
)

# Combine into master table
eda_interaction_master <- bind_rows(
  tab_inter_bonus_zone,
  tab_inter_bonus_km,
  tab_inter_bonus_km_merged,
  tab_inter_zone_km,
  tab_inter_zone_km_merged,
  tab_inter_bonus_make,
  tab_inter_bonus_make_merged
)

# Export table
write.csv(
  eda_interaction_master, 
  "output/tables/eda_interaction_summary.csv", 
  row.names = FALSE, 
  na = "",
  quote = FALSE
)


# 6. Model Specification -------------------------------------------------------
# Main Goal: High predictivity with understandable results. 
# Method: AIC as primary model-selection model and BIC as sensitivity criterion

n_indiv <- nrow(insurance) # 2182

# Null model (intercept and offset only)
model_null <- glm(
  Claims ~ 1 + offset(log(Insured)), 
  data = insurance, 
  family = poisson(link = "log")
)

# Functional Form Selection: Kilometres
mod_km_raw    <- glm(Claims ~ Kilometres + offset(log(Insured)), data = insurance, family = poisson(link = "log"))
mod_km_merged <- glm(Claims ~ Kilometres_merged + offset(log(Insured)), data = insurance, family = poisson(link = "log"))

AIC(mod_km_raw, mod_km_merged, k = log(n_indiv)) # 40270.62, 40285.86 as BIC -> Decision for raw version of Kilometres
AIC(mod_km_raw, mod_km_merged, k = 2) # 40242.18, 40268.80 as AIC -> Decision for raw version of Kilometres
# Interpretation: The raw version of Kilometres outperforms the merged version in BIC and AIC.
# We proceed with the raw version of Kilometres

# Functional Form Selection: Make
mod_make_raw    <- glm(Claims ~ Make + offset(log(Insured)), data = insurance, family = poisson(link = "log"))
mod_make_merged <- glm(Claims ~ Make_merged + offset(log(Insured)), data = insurance, family = poisson(link = "log"))

AIC(mod_make_raw, mod_make_merged, k = log(n_indiv)) # 40660.25, 40652.77 as BIC (low, but significant)-> Decision for raw version of Make
AIC(mod_make_raw, mod_make_merged, k = 2) # 40609.05, 40607.27 as AIC (very low delta)-> Decision for raw version of Make
# Interpretation: The merged version of Make outperforms the raw version in BIC, i.e., the additional parameters yield no significant increase in goodness-of-fit when taking the further complexity into account. The result of AIC is not significant by the very low delta of \sim 1.78, i.e., the predictive power of both models is approximately equal.
# The categories of Make represent distinct car models, whereby the detailed car - classification is not given. Regarding all 9 categories of Make allows a more precise tariff structure in practice and thus a more precise risk classification.
# We proceed with the raw version of Make


# Full model
model_full <- glm(Claims ~ Kilometres + Make + Zone + Bonus + offset(log(Insured)), data = insurance, family = poisson(link = "log"))


## Model selection Forward selection (step()-function) via BIC and AIC ---------------------
# BIC (forward- and backward- selection)
model_full_stepwise_forward_bic <- step(
  object = model_null, 
  direction = "forward", 
  scope = formula(model_full), 
  k = log(n_indiv),
  trace = 1
)
summary(model_full_stepwise_forward_bic) # Claims ~ Bonus + Zone + Kilometres + Make + offset(log(Insured)

model_full_stepwise_backward_bic <- step(
  object = model_full, 
  direction = "backward",
  k = log(n_indiv),
  trace = 1
)
summary(model_full_stepwise_backward_bic) # Claims ~ Kilometres + Make + Zone + Bonus + offset(log(Insured))

# Conclusion: Both methods yield Kilometres, Make, Zone and Bonus as covariates


# AIC (forward- and backward- selection)
model_full_stepwise_forward_aic <- step(
  object = model_null, 
  direction = "forward", 
  scope = formula(model_full), 
  trace = 1
)
summary(model_full_stepwise_forward_aic) # Claims ~ Bonus + Zone + Kilometres + Make + offset(log(Insured))

model_full_stepwise_backward_aic <- step(
  object = model_full, 
  direction = "backward",
  trace = 1
)
summary(model_full_stepwise_backward_aic) # Claims ~ Kilometres + Make + Zone + Bonus + offset(log(Insured))
# Conclusion: Kilometres, Make, Zone and Bonus as covariates

# Conclusion (overall): Forward and backward AIC / BIC yield exact the same covariates. While selection with AIC mainly focuses on prediction, BIC additionally evaluates the model complexity. This means the covariates Kilometres, Make, Zone and Bonus have adequate predictive power and removing one of them would not suffice for a worth complexity - precision trade-off in BIC score. 

model_final_without_interaction <-  model_full


# Interaction Effects

# Define interaction models by adding single interaction term to model_final_without_interaction

model_inter_bonus_zone <- update(model_final_without_interaction, . ~ . + Bonus:Zone)
model_inter_bonus_km <- update(model_final_without_interaction , . ~ . + Bonus:Kilometres)
model_inter_zone_km <- update(model_final_without_interaction, . ~ . + Zone:Kilometres)
model_inter_bonus_make <- update(model_final_without_interaction, . ~ . + Bonus:Make)

# Q: Does any specific interaction term significantly reduce the residual deviance?

anova(model_final_without_interaction, model_inter_bonus_zone, test = "Chisq") # A: yes, by p-value of 8.787e-10
anova(model_final_without_interaction, model_inter_bonus_km, test = "Chisq") # A: yes, by p-value of 2.2e-16
anova(model_final_without_interaction, model_inter_zone_km, test = "Chisq") # A: yes, by p-value of 2.957e-10
anova(model_final_without_interaction, model_inter_bonus_make, test = "Chisq") # A: yes, by p-value of 2.2e-16


## Interaction model selection Forward selection (step()-function) via BIC and AIC

model_full_inter <- update(model_final_without_interaction, . ~ . + Bonus:Zone + Bonus:Kilometres + Zone:Kilometres + Bonus:Make)

# BIC (forward- and backward- selection)
model_full_inter_stepwise_forward_bic <- step(
  object = model_final_without_interaction, 
  direction = "forward", 
  scope = formula(model_full_inter), 
  k = log(n_indiv),
  trace = 1
)
summary(model_full_inter_stepwise_forward_bic) # Claims ~ Kilometres + Make + Zone + Bonus + Kilometres:Bonus + offset(log(Insured)

model_full_inter_stepwise_backward_bic <- step(
  object = model_full_inter, 
  direction = "backward",
  k = log(n_indiv),
  trace = 1
)
summary(model_full_inter_stepwise_backward_bic) # Claims ~ Kilometres + Make + Zone + Bonus + Kilometres:Bonus + offset(log(Insured)

# Conclusion: Both methods yield the same result: Interaction term Kilometres:Bonus as additionally adequate main effect on the response.


# AIC (forward- and backward- selection)
model_full_inter_stepwise_forward_aic <- step(
  object = model_final_without_interaction, 
  direction = "forward", 
  scope = formula(model_full_inter), 
  trace = 1
)
summary(model_full_inter_stepwise_forward_aic) # Claims ~ Kilometres + Make + Zone + Bonus + Kilometres:Bonus + Make:Bonus + Zone:Bonus + Kilometres:Zone + offset(log(Insured))

model_full_inter_stepwise_backward_aic <- step(
  object = model_full_inter, 
  direction = "backward",
  trace = 1
)
summary(model_full_inter_stepwise_backward_aic) # Claims ~ Kilometres + Make + Zone + Bonus + Zone:Bonus + Kilometres:Bonus + Kilometres:Zone + Make:Bonus + offset(log(Insured))

# Conclusion: In addition to the selection based on BIC score, the step function additionally yield Zone:Bonus, Kilometres:Zone and Make:Bonus as adequate predictive for the response Claims.

# Conclusion (overall): Forward and backward AIC / BIC yield yield nested models 'model_full_inter_stepwise_backward_bic' \subseteq 'model_full_inter_stepwise_backward_aic'. In the EDA section, no significant visual indication of the existence of interaction effect was given. 
#Therefore, Likelihood-Ratio Tests are performed in the next step to formally test if additional parameter yield significant improve in fit.

# Q: Does the interaction term Kilometres:Bonus significantly reduce the residual deviance?
anova(model_final_without_interaction, model_full_inter_stepwise_backward_bic, test = "Chisq") # A: yes, by p-value of 2.2e-16

# Q: Does the additionally, the interaction terms Zone:Bonus, Kilometres:Zone and Make:Bonus significantly reduce the residual deviance?
anova(model_full_inter_stepwise_backward_bic, model_full_inter_stepwise_backward_aic, test = "Chisq") # A: yes, by p-value of 2.2e-16

# Compare BIC and AIC of all three nested models:

AIC(model_final_without_interaction, model_full_inter_stepwise_backward_bic, model_full_inter_stepwise_backward_aic, k = log(n_indiv)) # BIC 10796.20, 10784.97, 11192.85 -> \Delta_1 = 11.23, \Delta_2 = -407.88
# Model with interaction term Kilometres:Bonus outperforms the other models in BIC score significantly, i.e., adding the regarded interaction term as additive main effect is worth the precision - complexity tradeoff. According to the test, this is not true for the model containing in addition the covariates Zone:Bonus, Kilometres:Zone and Make:Bonus (i.e., all regarded interaction effects)

AIC(model_final_without_interaction, model_full_inter_stepwise_backward_bic, model_full_inter_stepwise_backward_aic, k = 2) # AIC 10654.00, 10506.25, 10299.83 -> \Delta_1 = 147.75, \Delta_2 = 206.64, 
# Model containing all interaction terms has significantly lowest AIC among the three nested models, i.e., highest predictive power.


# Conclusion: Including all interaction terms, i.e., following the AIC result, stays in contrast to the practical EDA, where significant interaction effects could not be distinguished from statistical noise of the effects Zone:Bonus, Kilometres:Zone and Make:Bonus. Note that the underlying data is aggregated into distinct covariate profiles, hence each cell consists of one observation. Therefore, the full interaction model is saturated, since #parameters = #cells. Consequently, the AIC result leads to overfitting. Moreover, the only explainable interaction effect in the EDA was Kilometres:Bonus, which aligns with the BIC result, hence this model is selected.

model_main <-  model_full_inter_stepwise_backward_bic
summary(model_main)

# Model Diagnostic ----------------------------

# Model integrity check (MLE property)
# In Poisson GLM with log-link and intercept, the score equation yields that sum of fitted values = sum of obs. values
sum_observed_claims <- sum(insurance$Claims)
sum_fitted_claims <- sum(fitted(model_main))

as.integer(sum_fitted_claims) == as.integer(sum_observed_claims) # TRUE


# Residual Deviance Test (Goodness of fit, comparison with saturated model): 
# Reject H_0: Model assumptions for specified GLM are satisfied vs. H_1: not H_0 at level of alpha

# Degrees of freedom
df_res <- df.residual(model_main)
# Test statistic
D_stat <- deviance(model_main)
dev_residuals <- resid(model_main, type = "deviance")
D_raw_stat <- sum(dev_residuals^2)

# Integrity check
tol <- 1.5e-8
abs(D_stat - D_raw_stat) < tol # TRUE

p_val_Deviance <- pchisq(D_stat, df = df_res, lower.tail = FALSE) # 1.996e-19 < 5%
# Interpretation: Reject the hypothesis that the poisson model assumptions for the regarded model are satisfied at level of 5%. 


# Rule of thumb: D > n - p (Overdispersion)?
D_stat / df_res # 1.30
# Observation: Overdispersion might be present, but need more precise checks


# Pearson Chi-Square Test
# Reject H_0: Specified Poisson model is adequate, i.e., equidispersion holds vs H_1: Not H_0
# Under H_0, the pearson \chi^2 statistic follows asymptotically (i.e., for large n) a \chi^2_{n-p} distribution.
pearson_residuals <- resid(model_main, type = "pearson") # pearson_residuals^2 = squared deviation / theoretical variance
X2_pearson <- sum(pearson_residuals^2)

# Calculate the one-sided p-value
p_val_pearson <- pchisq(X2_pearson, df = df_res, lower.tail = FALSE) # 3.69e-22 < 5%
# Conclusion: Reject the Null Hypothesis, that the specified model is adequate at significance level of 5%. According to the result of this hypothesis test, over- or underdispersion might be present. Due to the rule of thumb, overdispersion is more likely. 

# Overdispersion test, based on score / likelihood ratio test for GLMs (ref: https://www.math.cit.tum.de/fileadmin/w00ccg/math/Forschung/forschungsgruppen/statistics/academics/lec7.pdf, pp. 21-22)
# Dean (1992) score test for overdispersion
# Model: Additive random effect on log scale: \theta_\star  = \theta + Z_i, \theta = Offset_i + x_i \cdot \beta with  Z_i iid, E[Z_i] = 0 and Var(Z_i) = \tau \in \mathbb{R}_{\ge 0}.
# Reject H0: tau = 0 vs H1: tau > 0 at level \alpha, if T_s > \chi^2_{1, 1- \alpha}. 
# T_s is the score statistic and follows an asymptotic \chi_1^2 distribution under H_0.
# Use \chi_1^2 \sim N(0,1)^2: T_s = V^2, whereby V is the standardized score.

test_dean_overdispersion <- function(poisson_model, alpha = 0.05) {
  y_obs <- poisson_model$y 
  mu_hat <- fitted(poisson_model)
  
  # Calculate test statistic for additive random effect w.r.t. poisson model
  numerator <- 0.5 * sum((y_obs - mu_hat)^2 - mu_hat)
  denominator <- sqrt(0.5 * sum(mu_hat^2))
  # Standardized score
  V <- numerator / denominator
  T_s <- V^2
  
  # Asymptotic standard normal distribution:
  p_val <- pnorm(V, lower.tail = FALSE) # V has to be positive for overdispersion, one sided p - value
  z_critical <- qnorm(1 - alpha)
  
  # Asymptotic chi_1^2 distribution:
  p_val_chisq <- pchisq(T_s, lower.tail = FALSE, df = 1) # one sided test
  chisq_critical <- qchisq(1 - alpha, df = 1)
  
  # Output
  return(data.frame(
    Test_Statistic = T_s,
    Standardized_Score = V,
    P_Value_Normal = p_val,
    Overdispersion_Significant_Normal = p_val < alpha,
    Z_critical_Normal = z_critical,
    P_Value_Chisq = p_val_chisq,
    Overdispersion_Significant_Chisq = p_val_chisq < alpha & V > 0,
    Z_critical_Chisq = chisq_critical
  ))
}

test_dean_overdispersion(model_main)
# Conclusion: Reject H0 at 5% significance level. The test provides statistical evidence for presence of overdispersion in the poisson model.

# Handling Overdispersion
# Theoretical:
# It holds: \theta_\star  = \theta + Z_i \iff \mu_i^\star \coloneqq \exp(\theta_\star) = \mu_i * Z_i^\prime, whereby Z_i^\prime = \exp(Z_i) and \mu_i = \exp(\theta_i), i.e., Z_i^\prime is multiplicative effect on mean scale.
# Assumption: Z_i^\prime are iid Gamma(\nu, \nu) distributed, i.e., mean 1 and variance 1 / \nu .
# This implies (via integrating the conditional Poisson distribution over Gamma mixing distribution):
# Y_i \sim nb(\mu_i, \nu), in particular E[Y_i] = \mu_i and Var[Y_i] = \mu_i (1 + \mu_i / \nu)
# Note: Since \nu is unkown, this adjusted specification does not follow a GLM.

# Practical:
#glm.nb estimates \nu and regression coefficients via alternating MLE through Scores Function of NB-distribution.
model_nb <- MASS::glm.nb(
  formula = formula(model_main),
  data = insurance,
  link = log,
  x = TRUE
)

# Model comparison: Poisson vs NB (nested models)

# Instead of using deviance(model_main) - deviance(model_nb) for test statistic, use 2(LL_full - LL_reduced), which is mathematically equivalent, since D = -2(LL_Model - LL_Saturated).
# Reason: the R functions glm and glm.nb calculate Deviance not equally due to difference in parameters w.r.t. underlying family distribution (Poisson(\mu_i) vs NB(\mu_i, \nu))

log_like_poisson <- as.numeric(logLik(model_main))
log_like_NB <- as.numeric(logLik(model_nb))

D_stat_poisson <- deviance(model_main)
D_stat_NB <- deviance(model_nb)

log_like_poisson_sat <- (log_like_poisson + D_stat_poisson) / 2
log_like_NB_sat <- (log_like_NB + D_stat_NB) / 2

tol <- 1.5e-8
abs(log_like_poisson_sat - log_like_NB_sat) < tol # FALSE, result is 214.98


# Likelihood Ratio Test / Partial Deviance Test
# Reject H_0: Dispersion parameter 1 / \nu = 0 vs H_1: not H_0 at level \alpha
# Under H_0, the test statistic LRT_stat is asymptotically \chi^2_1 distributed.
LRT_stat <- 2 * (log_like_NB - log_like_poisson)

# p- value: df = 1, since NB has one additional parameter \nu; lower.tail = FALSE, since one-sided test
p_val_LRT <- pchisq(LRT_stat, df = 1, lower.tail = FALSE) # 3.25e-37 < 5%

# Conclusion: According to the hypothesis test, reject that the dispersion parameter is 0. This indicates that the NB - model provides a more adequate fit.

# Re-Calibration of Model Selection under Negative Binomial family ---------- < To Do>
# Re Do EDA with NB Model - note: CI differ, but same link function - short interpretation
# For this: Adjust helper function with additional method s.t. NB CIs are calculated.

# The initial model selection via step()-function (forward / backward, AIC / BIC, Raw vs Merged) was performed under the assumption that the response follows a poisson family distribution. But such a model suffers from overdispersion, which influences the test statistics and favors complex models to absorb the discrepancy between theoretical variance and standard deviation. Hence, the model selection is re-executed under the assumption that the response has NB distribution.

# Null model (intercept only)
model_null_NB <- MASS::glm.nb(Claims ~ 1 + offset(log(Insured)), data = insurance, link = log, x = T)

# Functional Form Check for NB

# Kilometres
nb_km_raw    <- update(model_null_NB, . ~ . + Kilometres)
nb_km_merged <- update(model_null_NB, . ~ . + Kilometres_merged)



# Make
nb_make_raw    <- update(model_null_NB, . ~ . + Make)
nb_make_merged <- update(model_null_NB, . ~ . + Make_merged)



























































# ==============================================================================
# 4.2 SEVERITY EDA (Gamma Regression)
# Objective: Visualize marginal effects on average claim size by plotting 
# empirical log-means and empirical pointwise 95% CIs.
# Note: Use 'insurance_severity' dataset (Claims > 0).
# ==============================================================================

# Helper function for Gamma EDA Plots
plot_gamma_eda <- function(df, grouping_var, var_label) {
  df %>%
    group_by({{ grouping_var }}) %>%
    summarise(
      P_j = sum(Payment),
      C_j = sum(Claims),
      # Empirical standard error of the log response for pointwise CI
      se_log_mean = sd(log(average_claim_size), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      emp_log_mean = log(P_j / C_j),
      ci_lower = emp_log_mean - qnorm(0.975) * se_log_mean,
      ci_upper = emp_log_mean + qnorm(0.975) * se_log_mean
    ) %>%
    ggplot(aes(x = {{ grouping_var }}, y = emp_log_mean)) +
    geom_point(color = "#c0392b", size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "#c0392b") +
    theme_light() +
    labs(x = var_label, y = "Empirical Log-Mean (Severity)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Helper function for Gamma EDA Tables
generate_gamma_table <- function(df, grouping_var, var_name) {
  df %>%
    group_by({{ grouping_var }}) %>%
    summarise(
      P_j = sum(Payment),
      C_j = sum(Claims),
      se = sd(log(average_claim_size), na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      `Empir. Log-Mean` = round(log(P_j / C_j), 3),
      `L` = round(`Empir. Log-Mean` - qnorm(0.975) * se, 3),
      `U` = round(`Empir. Log-Mean` + qnorm(0.975) * se, 3)
    ) %>%
    select({{ grouping_var }}, `Empir. Log-Mean`, P_j, C_j, L, U) %>%
    mutate(Variable = var_name, .before = 1)
}

# Generate Gamma Plots
p_sev_km    <- plot_gamma_eda(insurance_severity, Kilometres, "Kilometres (Distance Band)")
p_sev_bonus <- plot_gamma_eda(insurance_severity, Bonus, "Bonus (Claim-Free Years)")
p_sev_zone  <- plot_gamma_eda(insurance_severity, Zone, "Zone (Geographical Region)")
p_sev_make  <- plot_gamma_eda(insurance_severity, Make, "Make (Car Model)")

ggsave("output/figures/eda_sev_kilometres.png", plot = p_sev_km, width = 6, height = 4, dpi = 300)
ggsave("output/figures/eda_sev_bonus.png", plot = p_sev_bonus, width = 6, height = 4, dpi = 300)
ggsave("output/figures/eda_sev_zone.png", plot = p_sev_zone, width = 6, height = 4, dpi = 300)
ggsave("output/figures/eda_sev_make.png", plot = p_sev_make, width = 6, height = 4, dpi = 300)

# Generate Gamma Master Table
tab_sev_km    <- generate_gamma_table(insurance_severity, Kilometres, "Kilometres")
tab_sev_bonus <- generate_gamma_table(insurance_severity, Bonus, "Bonus")
tab_sev_zone  <- generate_gamma_table(insurance_severity, Zone, "Zone")
tab_sev_make  <- generate_gamma_table(insurance_severity, Make, "Make")

eda_sev_master <- bind_rows(tab_sev_km, tab_sev_bonus, tab_sev_zone, tab_sev_make)
write.csv(eda_sev_master, "output/tables/eda_severity_summary.csv", row.names = FALSE, na = "")

# Interpretation Placeholder (Gamma):
# - Evaluate if variables driving frequency (e.g., Bonus) have the same, opposite, or no effect on severity.
# - Merge categories structurally similarly to frequency analysis if CI widths are large due to sparse claim counts (C_j).

# 5. Data Refinement: Category Merging based on EDA ---------------------------
# [Hier fügst du dann deine mutate()-Befehle ein, um Level zusammenzufassen, 
#  basierend auf den Outputs der Tabellen und CIs]





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