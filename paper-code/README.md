This file contains the code for reproducing results in the paper:

Oh-Ran Kwon and Hui Zou (2024) **Enhanced Response Envelope via Envelope Regularization**, *Journal of the American Statistical Association*, 1-22. DOI: [10.1080/01621459.2024.2368844](https://doi.org/10.1080/01621459.2024.2368844)

## Structure

- `simulation.R` contains code for simulation.
- `real-data.R` contains code for real data analysis.

## Running the code

- `simulation.R`: Three arguments — n (sample size), p (number of predictors), and rho (value controlling the covariance of predictors) — should be passed to the R script from the command line. For example, the following command lines
```
Rscript simulation.R 50 5 0.0 > “sim_n50_p5_rho0.Rout”
Rscript simulation.R 50 55 0.8 > “sim_n50_p55_rho0.8.Rout”
```
output the results for the cases of n=50, p=5, rho=0 and n=50, p=55, rho=0.8, respectively. The outputs are saved in the corresponding `.Rout` files. 

- `real-data.R`: Run the following command line to execute the script.
```
R --vanilla < "real-data.R" > "real-data.Rout"
```
The output file for this is `real-data.Rout`. 
