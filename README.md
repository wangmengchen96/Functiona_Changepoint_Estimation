# Functional_Changepoint_Estimation

Source code and data for the Bayesian hierarchical model, which is proposed in the paper **Asynchronous Changepoint Estimation for Spatially Correlated Functional Time Series**.

## Data

These are under the folder, "covid_example".

`covidData.RData` is the matrix for the daily COVID-19 cases for all counties in Illinois between 01/01/2021 and 04/05/2021 from the Illinois Department of Public Health (https://www.dph.illinois.gov/covid19/data-portal). This is the dataset after removing the "Unknown" age group.

`CountyName.RData` is a vector with a length of 102. It contains the name of all counties in Illinois, which are the column names in `covidData.RData` correspondingly.

`countyCenter.RData` saves the latitude and longitude of county centers.

`ages.RData` are the ages groups for each row of `covidData.RData`.

`dates.RData` are the dates for each row of `covidData.RData`.

## Code

`Functional_Changepoint_Estimation_Simulation.R` contains the code of the data generation for simulation, the estimation results from two previous methods, and the implementation of the proposed Bayesian hierarchical model. Stepsize may need to be tuned for other datasets. In the code with 50 locations and 50 time points, it takes around 3 hours to run 20,000 iterations.

`covid_data_preprocess.R` in Folder "covid_example" shows how we scale the spatial domain and functional data, getting the rejection region and the initial values of MCMC. Model building is the same as that in the simulation code. Stepsize needs to be tuned to get the ideal acceptance rate.

Versions of software and loaded packages:
* `R 4.1.0`
* `mvtnorm_1.1-1`
* `tmvtnorm_1.4-10`
* `sde_2.0.15`
* `fChange_0.2.1`
* `invgamma_1.1`
* `reshape2_1.4.4`
* `ggplot2_3.3.5`
* `scpt_0.1`
* `splines_4.1.0`
* `matrixcalc_1.0-3`
