### inc 
library(tidyverse)
library(cowplot)

### settings
theme_set(theme_cowplot(font_size = 10))

### par
env <- readRDS("rdata/env_rhoe.rds")

envrep <- readRDS("rdata/env_rep_rhoe.rds")
levs <- unique(sort(envrep$n))
r <- unique(envrep$nrep)
envrep <- mutate(envrep, n = factor(n, levs, paste(levs, "x", r)))

gen <- readRDS("rdata/gen_rhog.rds")

genrep <- readRDS("rdata/gen_rep_rhog.rds") %>%
  bind_rows(readRDS("rdata/gen_rep_rhog-2000.rds")) %>%
  bind_rows(readRDS("rdata/gen_rep_rhog-4000.rds"))
  
levs <- unique(sort(genrep$n))
r <- unique(genrep$nrep)
genrep <- mutate(genrep, n = factor(n, levs, paste(levs, "x", r)))

### plots
update_plot <- function(p, ylab = "") 
{
  p + geom_boxplot() + facet_wrap(~ n, nrow = 1) + 
    geom_abline(linetype = 3) + coord_equal() + 
    labs(x = "", y = ylab)
}

p1 <- ggplot(env, aes(rhoe, rhoe_hat, group = as.factor(rhoe))) %>% update_plot("rhoe: Env")

p2 <- ggplot(envrep, aes(rhoe, rhoe_hat, group = as.factor(rhoe))) %>% update_plot("rhoe: Env + Rep")

p3 <- ggplot(gen, aes(rhog, rhog_hat, group = as.factor(rhog))) %>% update_plot("rhog: Env + Gen")

p4 <- (ggplot(genrep, aes(rhog, rhog_hat, group = as.factor(rhog))) + geom_hline(yintercept = 0, linetype = 3) + xlim(c(-1, 1)) + ylim(c(-1, 1))) %>% update_plot("rhog: Env + Gen + Rep")

### plot 1
plot_grid(p1, p2, p3, ncol = 1)
#ggsave("rho-sim-v1.png", dpi = 150)

### plot 2
p4
#ggsave("rho-sim-4-v1.png", dpi = 150)

