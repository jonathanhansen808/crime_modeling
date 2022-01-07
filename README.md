# crime_modeling
This project models crime rates of 88 counties in North Carolina from 1981 to 1987 using Bayesian methodologies with R and STAN. 

The Data is from the Ecdat R package Crime: Crime in North Carolina (https://rdrr.io/cran/Ecdat/man/Crime.html). 

The explanatory variables selected for analysis were probability of arrest, probability of a prison sentence, the population density of a county (100s of people),
and the ratio of a county's police force density relative to the North Carolina average. 

The response variable was a county's crime rate. 

All variables were transformed logarithmically. Additionally, year 1987 was left out as test data to see if the model was accurate.  

Two models were fit, with the first one being a bayesian log-log regression that pools all 88 counties crime rates as a baseline county for the intercept. The other model was a hierarchical log-log regression model that gives each county its own random intercept. The latter resulting in a better fit, as a result of letting each county have its own unique intercept for crime rate. 

All priors were based on initial visualizations and prior predictive checks were carried out to ensure that the choice of prior(s) were reasonable. 
After the models were fit, posterior predictive checks using the last year (1987) were carried out and the hierarchical model captured more county's true crime rates, than the pooled model, with 95% prediction intervals. 

The anlysis (R) and model (STAN) are included above. 




