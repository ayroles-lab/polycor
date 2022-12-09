uniRandSlopes <- cmdstan_model(here::here("stan/models/univariateRandomSlopes.stan"))
uniRandSlopesiid <- cmdstan_model(here::here("stan/models/univariateiidRandomSlopes.stan"))
uniRandSlopes <- cmdstan_model(here::here("stan/models/bivariateRandomSlopes.stan"))
uniRandSlopesiid <- cmdstan_model(here::here("stan/models/bivariateiidRandomSlopes.stan"))
