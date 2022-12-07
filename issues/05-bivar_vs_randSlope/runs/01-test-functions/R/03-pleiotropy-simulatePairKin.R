### inc 
library(tidyverse)
library(magrittr)

library(Matrix)
library(MASS)

library(nlme)
library(lme4)

# change to your path here, e.g. based on `nodemame`
path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
src <- file.path(path, "relmer.R")
source(src)

src <- file.path(path, "simFunctions.R")
source(src)

### par
testing <- FALSE

h2 <- 0.75
nrep <- 4
vals_n <- c(100, 500, 1e3) # / nrep
vals_rhog <- c(0.1, 0.5, 0.9)
vals_sim <- 1:5

if(testing) {
  nrep <- 2
  vals_n <- c(10, 50)
}

### run: no rep.
grid <- expand.grid(n = as.integer(vals_n), rhog = vals_rhog, sim = vals_sim) %>%
  as_tibble

tab <- lapply(seq(nrow(grid)), function(i)
{
  n <- grid$n[i]
  rhog <- grid$rhog[i]
  sim <- grid$sim[i]
  
  cat(" -", i, "/", nrow(grid), "\n")
  print(grid[i, ])
  
  simdat <- pleiotropy.simulatePairKin(n_ind = n, n_perInd = nrep, 
    h2 = h2, overlap = rhog)

  dat <- simdat$pheno %>% 
    mutate(id = as.character(id), rid = id)
  G <- simdat$G


  dat_slope <- mutate(dat, trait1 = as.numeric(scale(trait1)))
    
  vals_models <- c("bivar", "slope", "slope_add", "add")
    tab <- lapply(seq(length(vals_models)), function(k) {
      model <- vals_models[k]
      
      tab <- switch(model,
        "bivar" = {
          bdat <- gather(dat, tname, tvalue, trait1, trait2)
          mod <- relmer(tvalue ~ -1 + tname + obs + (0 + tname | id) + (0 + tname | rid), 
            data = bdat, relmat = list(id = G), calc.derivs = FALSE)
          
          vf <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id")
          cov <- filter(vf, grp == "id" & var1 == "tnametrait1" & var2 == "tnametrait2") %$% vcov
          cov <- cov / sum(vf$vcov)
          
          data_frame(cov = cov, add = NA)
        },
        "slope" = {
          mod <- relmer(trait2 ~ obs + (0 + trait1|id), data = dat_slope, relmat = list(id = G))
          cov <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop

          data_frame(cov = cov, add = NA)
        },
        "slope_add" = {
          mod <- relmer(trait2 ~ obs + (1|rid) + (0 + trait1|id), data = dat_slope, 
            relmat = list(id = G, rid = G))
        
          cov <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
          add <- VarProp(mod) %>% as.data.frame %>% filter(grp == "rid") %$% prop
          
          data_frame(cov = cov, add = add)
        },
        "add" = {
          mod <- relmer(trait2 ~ obs + (1|id), data = dat_slope, relmat = list(id = G))
          add <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop

          data_frame(cov = NA, add = add)        
        },
        stop("error in switch: unknown `model`"))
        
        mutate(tab, model = model)
     }) %>% bind_rows

  tab %>% mutate(n = n, nrep = nrep, sim = sim,
    h2 = h2, rhog = rhog, 
    h2_cov = (rhog * h2)/2 
    # divided by 2, because covariance is shared equally between two traits
    # (the h2 is the same for two traits)
    
 )
}) %>% bind_rows

### plots
n2str <- function(n) paste("N:", n, "x", nrep)

ptab <- mutate(tab,
  n = factor(n2str(n), levels = n2str(vals_n)))

p1 <- ggplot(ptab, aes(as.factor(rhog), add, color = model)) + geom_boxplot() + facet_wrap(~ n) + theme_minimal() + theme(legend.position = "top")


p2 <- ggplot(ptab, aes(as.factor(rhog), cov, color = model)) + geom_boxplot() + facet_wrap(~ n) + geom_hline(aes(yintercept = h2_cov), linetype = 3) + theme_minimal() + theme(legend.position = "top")

