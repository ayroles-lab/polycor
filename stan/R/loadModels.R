uniRandSlopes <- cmdstan_model(here::here("stan/models/univariateRandomSlopes.stan"))
uniRandSlopesiid <- cmdstan_model(here::here("stan/models/univariateiidRandomSlopes.stan"))
biRandSlopes <- cmdstan_model(here::here("stan/models/bivariateRandomSlopes.stan"))
biRandSlopesiid <- cmdstan_model(here::here("stan/models/bivariateiidRandomSlopes.stan"))
