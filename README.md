# Comparative-Analysis-of-Correlation-Tests-in-R
This GitHub project empirically compares the power of Pearson, Spearman, and Kendall correlation tests under different bivariate distributions
The analysis demonstrates:

Pearson's test is more powerful for bivariate normal distributions
Nonparametric tests can outperform Pearson for certain non-normal distributions
The project includes extensive R code with:

Simulation functions, Power calculations, Visualization,and Statistical analysis

Key Findings Demonstration

For the Bivariate Normal Distribution, Pearson’s test demonstrates superior power, correctly rejecting the null hypothesis when correlation exists, whereas nonparametric tests exhibit lower power in this ideal Pearson case. However, for Non-Normal Distributions, nonparametric tests prove more robust. In the case of heavy-tailed distributions, such as the bivariate t-distribution, nonparametric methods maintain higher power compared to Pearson. When analyzing quadratic relationships, Pearson’s test fails to detect dependence, while rank-based methods successfully capture the relationship. Similarly, for the Frank copula, which exhibits nonlinear dependence, rank-based methods outperform Pearson, highlighting their effectiveness in capturing non-monotonic associations.
