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

h2_cov <- 0.4
nrep <- 2
vals_n <- c(100, 500, 1e3) / nrep
vals_h2_add <- c(0, 0.2, 0.4)
vals_sim <- 1:5

if(testing) {
  nrep <- 2
  vals_n <- c(50, 100)
}

### run: no rep.
grid <- expand.grid(n = as.integer(vals_n), h2_add = vals_h2_add, sim = vals_sim) %>%
  as_tibble

tab <- lapply(seq(nrow(grid)), function(i)
{
  n <- grid$n[i]
  h2_add <- grid$h2_add[i]
  sim <- grid$sim[i]
  
  cat(" -", i, "/", nrow(grid), "\n")
  print(grid[i, ])
  
  simdat <- genCov.simulatePairKin(n_ind = n, n_perInd = nrep, 
    h2.cov = h2_cov, h2.add = h2_add)
  
  dat <- simdat$pheno %>% 
    mutate(id = as.character(id), rid = id, trait1 = as.numeric(scale(trait1)))
  G <- simdat$G

  vals_simtrait <- c("slope", "slope_add", "add")
  vals_models <- c("bivar", "slope", "slope_add", "add")
  tab <- lapply(seq(length(vals_simtrait)), function(j) {
    simtrait <- vals_simtrait[j]
    
    dat <- switch(simtrait,
      "slope" = dat,
      "slope_add" = within(dat, trait2 <- trait2_withAdd),
      "add" = within(dat, trait2 <- trait2_onlyAdd),
      stop("error in switch: unknown 'simtrait'"))
    
    tab <- lapply(seq(length(vals_models)), function(k) {
      model <- vals_models[k]
      
      tab <- switch(model,
        "bivar" = {
          bdat <- gather(dat, tname, tvalue, trait1, trait2)
          mod <- relmer(tvalue ~ -1 + tname + obs + (0 + tname | id) + (0 + tname | rid), 
            data = bdat, relmat = list(id = G), calc.derivs = FALSE)
          
          cov <- VarProp(mod) %>% as.data.frame %>% 
            filter(grp == "id" & var1 == "tnametrait1" & var2 == "tnametrait2") %$% prop
          
          data_frame(cov = cov, add = NA)
        },
        "slope" = {
          mod <- relmer(trait2 ~ obs + (0 + trait1|id), data = dat, relmat = list(id = G))
          cov <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop

          data_frame(cov = cov, add = NA)
        },
        "slope_add" = {
          mod <- relmer(trait2 ~ obs + (1|rid) + (0 + trait1|id), data = dat, 
            relmat = list(id = G, rid = G))
        
          cov <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop
          add <- VarProp(mod) %>% as.data.frame %>% filter(grp == "rid") %$% prop
          
          data_frame(cov = cov, add = add)
        },
        "add" = {
          mod <- relmer(trait2 ~ obs + (1|id), data = dat, relmat = list(id = G))
          add <- VarProp(mod) %>% as.data.frame %>% filter(grp == "id") %$% prop

          data_frame(cov = NA, add = add)        
        },
        stop("error in switch: unknown `model`"))
        
        mutate(tab, model = model)
     }) %>% bind_rows
     
     mutate(tab, simtrait = simtrait)
   }) %>% bind_rows
  
  tab %>% mutate(n = n, nrep = nrep, sim = sim,
    h2_add = h2_add, h2_cov = h2_cov)
}) %>% bind_rows

### plots
n2str <- function(n) paste("N:", n, "x", nrep)

ptab <- mutate(tab,
  n = factor(n2str(n), levels = n2str(vals_n)),
  simtrait = paste("sim:", simtrait))

p1 <- ggplot(ptab, aes(as.factor(h2_add), add, color = model)) + geom_boxplot() + facet_grid(simtrait ~ n) + geom_hline(aes(yintercept = h2_add), linetype = 3) + theme_minimal() + theme(legend.position = "top")

p2 <- ggplot(ptab, aes(as.factor(h2_add), cov, color = model)) + geom_boxplot() + facet_grid(simtrait ~ n) + geom_hline(aes(yintercept = h2_cov), linetype = 3) + theme_minimal() + theme(legend.position = "top")

