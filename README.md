# Frequency-Severity Modeling for Insurance Claims

## Overview
This repository demonstrates the end-to-end development of a Frequency-Severity model to analyze the impact of mileage, region, bonus class, and vehicle type on claim frequency and severity.

## Mathematical Framework
The modeling approach is split into two generalized linear models (GLMs):

### 1. Claim Frequency (Poisson Regression)
To model the number of claims $N_i$ for a given exposure (insured duration) $t_i$, a Poisson distribution with a log-link function and exposure offset is utilized:
$$N_i \sim \text{Poisson}(\mu_i)$$
$$\log(\mu_i) = \log(t_i) + \beta_0 + \mathbf{x}_i^T \boldsymbol{\beta}$$

### 2. Claim Severity (Gamma Regression)
The average claim size $Z_i = \frac{\text{Payment}}{\text{Claims}}$ is modeled strictly for positive claims using a Gamma distribution with a log-link function:
$$Z_i \sim \text{Gamma}(\alpha, \lambda_i)$$
$$\log(E[Z_i]) = \gamma_0 + \mathbf{x}_i^T \boldsymbol{\gamma}$$

## Project Structure
- `scripts/`: Contains the three main R scripts for preprocessing, Poisson frequency modeling, and Gamma severity modeling.
- `data/`: Contains the processed datasets (`.rds` files).
- `output/`: Contains model objects and diagnostic visualizations.

---
*Developed as part of a statistical modeling project.*
