# Frequency-Severity Modeling for Insurance Claims

## Overview
This repository shows an actuarial **Frequency-Severity model** to analyze the impact of mileage, region, bonus class, and vehicle type on insurance claim frequencies and claim sizes.


## Variable description
\textbf{Source}: Swedish Committee on the Analysis of Risk Premium in Motor Insurance (compiled for the year 1977).

\textbf{Description:} The data give details of third party motor insurance claims in Sweden for the year 1977. In Sweden, all motor insurance companies apply identical risk arguments to classify customers, allowing their portfolios and claims statistics to be combined. The data were compiled by the Committee to analyze the real influence of these risk arguments on claims and to compare this empirical structure with the actual tariff. The number of claims in each category can be treated as Poisson distributed to a good approximation.

The variables \texttt{Kilometres}, \texttt{Zone}, \texttt{Bonus}, and \texttt{Make} represent the a priori risk arguments used for customer classification. The variable \texttt{Insured} is the measure of the risk exposure of the insurance companies.  

More details on the variables  

\texttt{Kilometres} = Kilometres travelled per year  

1 = < 1000  

2 = 1000-15000  

3 = 15000-20000  

4 = 20000-25000

5 = > 25000

\texttt{Zone} = Geographical zone

1 = Stockholm, Göteborg, Malmö with surroundings

2 = Other large cities with surroundings

3 = Smaller cities with surroundings in southern Sweden

4 = Rural areas in southern Sweden

5 = Smaller cities with surroundings in northern Sweden

6 = Rural areas in northern Sweden

7 = Gotland

\texttt{Bonus} = No claims bonus

Equal to the number of years, plus one, since the last claim.

\texttt{Make} = Car model

1-8 = Eight different common car models. Make 4 represents the Volkswagen 1200 (discontinued shortly after 1977). The other makes remain unidentified to prevent potential impacts on the sales of those cars.

9 = All other models combined.  

\texttt{Insured} = Number of insured in policy-years  

\texttt{Claims} = Number of claims  

\texttt{Payment} = Total value of payments in Skr (Swedish Krona)

#----------------------------------
## Mathematical Framework
The modeling approach is split into two generalized linear models (GLMs):

### 1. Claim Frequency (Poisson Regression)
To model the number of claims $N_i$ for a given exposure (insured duration) $t_i$, a Poisson distribution with a log-link function and exposure offset is utilized:

> $$N_i \sim \text{Poisson}(\mu_i)$$
> $$\log(\mu_i) = \log(t_i) + \beta_0 + \mathbf{x}_i^T \boldsymbol{\beta}$$

### 2. Claim Severity (Gamma Regression)
The average claim size $Z_i = \frac{\text{Payment}}{\text{Claims}}$ is modeled strictly for positive claims using a Gamma distribution with a log-link function:

> $$Z_i \sim \text{Gamma}(\alpha, \lambda_i)$$
> $$\log(E[Z_i]) = \gamma_0 + \mathbf{x}_i^T \boldsymbol{\gamma}$$

## Project Structure
- `scripts/`: Contains the three main R scripts for preprocessing, Poisson frequency modeling, and Gamma severity modeling.
- `data/`: Contains the processed datasets (`.rds` files).
- `output/`: Contains model objects and diagnostic visualizations.

---
*Developed as part of a statistical modeling project.*
