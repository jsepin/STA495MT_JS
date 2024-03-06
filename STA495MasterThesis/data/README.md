# data:
---

The data sets in this directory contain simulated data sets that we need to discuss the different simulation approaches to create design matrices with a given collinearity constellation:

* *est_rmv.rds* : Resulting estimands when simulating over multivariate normal distribution
* *est_rmv_normal.rds* : Resulting estimands when simulating over multivariate normal distribution with subsequent transformation to uniform distribution.
* *est_scalefac.rds* : Resulting estimands when simulating over the scale factor method
* *cross_section.rds* : Contains the data set to illustrates how the differently constructed look like for two different condition numbers
* *df_min_trouble.rds* : Contains the data set to with the highest collinearity magnitude we allow. This data set is used to determine the number of how many time we simulate each experimental condition.

Furthermore, it also provides the decided experimental conditions:

* *boston_parameters.rds* 

