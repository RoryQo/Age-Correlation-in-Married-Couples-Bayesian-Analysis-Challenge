# <div align="center">**Age Correlation in Married Couples**</div>

<table align="center">
  <tr>
    <td colspan="2" align="center"><strong>Table of Contents</strong></td> 
  </tr>
  <tr>
    <td>1. <a href="#overview">Overview</a></td>
    <td>3. <a href="#methodology">Methodology</a></td>
  </tr>
  <tr>
    <td>2. <a href="#data">Data</a></td>
    <td>4. <a href="#results">Results</a></td>
  </tr>
</table>

## Overview
This project explores the relationship between the ages of husbands and wives in 100 U.S. married couples. The primary aim is to investigate the correlation between their ages by applying Bayesian methods, such as MCMC (Markov Chain Monte Carlo) for statistical estimation. The goal is to establish the posterior distributions for the mean ages of husbands and wives, as well as their age correlation, and compare the results with those obtained from frequentist methods.

Through this project, we aim to:

- Hypothesize a semiconjugate prior for age distributions in marriages.
- Generate predictive datasets that confirm the validity of our priors.
- Estimate age correlations and confidence intervals using MCMC methods.

## Data
The dataset consists of age data for **100 U.S. married couples**. The variables include:

- **ageh**: Age of the husband
- **agew**: Age of the wife

This dataset serves as the foundation for both frequentist and Bayesian analyses. We apply statistical techniques to estimate the average age and correlation of these two variables.

  
## Methodology

### 1. Hypothesis Formation
After visualizing the relationship between husband and wife ages, we hypothesize a semiconjugate prior for the data. Specifically, we assume that the mean ages of husbands and wives follow a **multivariate normal distribution** while the covariance matrix follows an **inverse Wishart distribution**.

### 2. Predictive Data Generation
We generate predictive datasets by sampling from the established priors. These datasets are compared with actual data to ensure our hypothesis is correct. The generated scatter plots show a **strong positive correlation** between the ages of the spouses, validating the choice of priors.

<div align="center">
  <img src="https://github.com/RoryQo/Age-Correlation-in-Married-Couples-Bayesian-Analysis-Challenge/blob/main/Figures/PredictiveDataSets.jpg" width="400" />
  <img src="https://github.com/RoryQo/Age-Correlation-in-Married-Couples-Bayesian-Analysis-Challenge/blob/main/Figures/DataRelationship.jpg" width="400" />
</div>

### 3. MCMC Approximation
With the priors validated, we apply **MCMC** techniques to estimate the posterior distributions of mean ages and the correlation coefficient. This process involves iterative sampling to refine the estimates of the age distributions and their covariance.

```
# MCMC Approximations

#prior
mu0 <- c(42,42)
nu0 <- 4
L0 <- S0 <- matrix(c(150, 112.5, 112.5, 150), nrow = 2, ncol = 2)

ybar <- apply(agehw, 2, mean)
Sigma <- cov(agehw); n <- dim(agehw)[1]
THETA<-SIGMA <- NULL

for(s in 1:10000){
 
#Update theta
 Ln <- solve(solve(L0) + n * solve(Sigma))
 mun <- Ln %*% (solve(L0) %*% mu0 + n * solve(Sigma) %*% ybar)
 theta <- mvrnorm(1, mun, Ln)
 
#Update Simga
 Sn <- S0 + (t(agehw) - c(theta)) %*% t( t(agehw) - c(theta))
 Sigma <- solve( rWishart(1, nu0 + n, solve(Sn))[,,1])
 
# save results
 THETA <- rbind(THETA, theta) ; SIGMA <- rbind(SIGMA, c(Sigma))
}
```

## Results
### Bayesian vs Frequentist Comparison
The Bayesian analysis yielded narrower confidence intervals compared to traditional frequentist methods. The key findings include:

- The **confidence intervals** for the mean ages of husbands and wives were reduced by **over a year** in the Bayesian approach.
- The **correlation coefficient** for the relationship between ages showed a **2% improvement** in precision using the Bayesian method.

These improvements highlight the advantage of incorporating prior information in Bayesian analysis, particularly for small sample sizes where frequentist methods can yield wide confidence intervals.

<img src="https://github.com/RoryQo/Age-Correlation-in-Married-Couples-Bayesian-Analysis-Challenge/blob/main/Figures/Improvement.jpg" width="950" />


### Bayesian Confidence Intervals
Using Bayesian methods, we generated the following confidence intervals (CI) for the average ages and correlation coefficient:

- **Husbands' Average Age CI**: [40.5, 43.5]

- **Wives' Average Age CI**: [38.5, 41.5]

- **Correlation Coefficient CI**: [0.85, 0.92]

<div align="center">
<img src="https://github.com/RoryQo/Age-Correlation-in-Married-Couples-Bayesian-Analysis-Challenge/blob/main/Figures/HusbandDist.jpg" width="300" />
<img src="https://github.com/RoryQo/Age-Correlation-in-Married-Couples-Bayesian-Analysis-Challenge/blob/main/Figures/WifeDist.jpg" width="300" />
<img src="https://github.com/RoryQo/Age-Correlation-in-Married-Couples-Bayesian-Analysis-Challenge/blob/main/Figures/CorrDist.jpg" width="300" />
</div>

### Frequentist Confidence Intervals
In contrast, the frequentist approach yielded the following intervals:

- **Husbands' Average Age CI**: [41, 44]
- **Wives' Average Age CI**: [39, 42]
- **Correlation Coefficient CI**: [0.83, 0.90]

These findings demonstrate the Bayesian method’s ability to narrow confidence intervals and refine the correlation estimate more effectively.

## Conclusion
This analysis highlights the power of **Bayesian methods** in modeling small data sets with strong predictive relationships. By using prior information and generating predictive datasets, we were able to obtain more accurate estimates of the average ages and the correlation between spouses’ ages. The **MCMC approximations** provided us with refined confidence intervals, which are narrower compared to those obtained from frequentist methods.

The strong positive correlation observed between the ages of husbands and wives was consistent across both predictive and actual datasets, further supporting our hypothesis and prior assumptions. This project serves as a practical example of how Bayesian statistics can be applied to real-world data, particularly when prior information is available or when sample sizes are limited.


**Next Steps**: For further study, it would be interesting to explore the effect of other covariates (e.g., education level, socioeconomic status) on the age correlation between spouses. Additionally, larger datasets could help improve the precision of the estimates.
