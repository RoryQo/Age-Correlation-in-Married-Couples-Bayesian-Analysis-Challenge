# <div align="center">**Age Correlation in Married Couples**</div>

## Table of Contents 
1. [Overview](#overview) 
2. [Methodology](#methodology) 
3. [Data](#data)
4. [Results](#results)

## Overview
This project analyzes the relationship between the ages of husbands and wives using data from 100 sampled U.S. married couples. The goal is to hypothesize a semiconjugate prior, generate a predictive dataset to confirm this hypothesis and utilize Markov Chain Monte Carlo (MCMC) approximation to estimate the mean ages and the correlation coefficient.

## Methodology
1. **Hypothesis Formation**: After visualizing the relationship between the ages of husbands and wives, we hypothesize a semiconjugate prior. Specifically, we model the variation using an inverse Wishart distribution, while the mean ages of spouses follow a multivariate normal distribution.

2. **Predictive Data Generation**: We generate predictive datasets by taking random samples from the established distributions. These datasets help us confirm the validity of our priors by comparing simulated scatter plots to the actual data. Both follow a strong positive correlation

<img src="https://raw.githubusercontent.com/RoryQo/R-Age-Correlation-in-Married-Couples-Bayesian-Analysis-Mini-Project/main/PredictiveDataSets.jpg" width="400" />


4. **MCMC Approximation**: Once the distributions are confirmed, we apply MCMC methods to estimate the mean correlation between the ages of husbands and wives. This allows us to create new confidence intervals for average ages and correlations.

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

## Data
The dataset comprises the ages of spouses from 100 sampled U.S. married couples:
- **ageh**: Age of the husband
- **agew**: Age of the wife

## Results
The limited sample size and lack of prior information about the variables led to a relatively wide confidence interval using standard frequentist approaches. In contrast, the Bayesian approach reduced the confidence interval for the average ages of husbands and wives by over a year. The confidence interval for the correlation coefficient was also reduced by 2% compared to the frequentist method (at the same alpha level).

These findings highlight the advantages of Bayesian methods over traditional frequentist approaches, particularly in providing narrower confidence intervals and better estimations of correlations. The ability to leverage prior information significantly enhances the reliability of our analysis regarding age dynamics in marriages.
