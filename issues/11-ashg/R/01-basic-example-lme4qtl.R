### inc
library(tidyverse)

stopifnot(packageVersion("lme4qtl") >= "0.2.2")
library(lme4qtl)

### data
data(dat40, package = "lme4qtl") # -> dat40, kin2

dat <- as_tibble(dat40) %>%
  mutate(
    trait1 = as.numeric(scale(trait1)), trait2 = as.numeric(scale(trait2)),
    id = as.character(ID), id2 = ID, rid = ID) 
G <- kin2 

### output
tab <- tibble()

### model `add`
# cov(y2) = g1^2 + e1^2
mod <- relmatLmer(trait2 ~ (1|id2), data = dat, relmat = list(id2 = G))
add <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id2") %$% prop

tab <- bind_rows(tab, tibble(model = "add", add = add))

### model `bivar` 
# cov([y1, y2]) = [g1^2, rho_g g1 g2; rho_g g1 g2, g2^2] + 
#   [e1^2, rho_e e1 e2; rho_e e1 e2, e2^2] + 1e-10 e^2
bdat <- gather(dat, tname, tvalue, trait1, trait2)
mod <- relmatLmer(tvalue ~ -1 + tname + (0 + tname | id) + (0 + tname | rid), 
  data = bdat, relmat = list(id = G), weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

cov <- VarProp(mod) %>% as.data.frame %>% 
  filter(grp == "id" & var1 == "tnametrait1" & var2 == "tnametrait2") %$% prop

tab <- bind_rows(tab, tibble(model = "bivar", cov = cov))

### model `slope_add`
mod <- relmatLmer(trait2 ~ (1|id2) + (0 + trait1|id) + (0 + trait1|rid), data = dat, 
  relmat = list(id2 = G, id = G), calc.derivs = FALSE)

add <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id2") %$% prop
icov <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
          
tab <- bind_rows(tab, tibble(model = "slope_add", add = add, icov = cov))
print(tab)

### print results
print(tab)
