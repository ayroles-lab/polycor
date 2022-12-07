### inc
load_all("~/git/hemostat/sumstats/")

library(tidyverse)

### par
thr_pval <- 1e-3

### variables
# results from `ss_ukbb_search_meta("fat")`
fat <- "100004_irnt"
vitd <- "100021_irnt"

### get assoc
assoc_fat <- ss_ukbb_assoc(fat, pval = thr_pval)

assoc_vitd <- ss_ukbb_assoc(vitd, pval = thr_pval)

### plots
title <- paste0("Fat, ", fat, ", n = ", assoc_fat$n_complete_samples[1])  
p_fat <- ggplot(assoc_fat, aes(pos, -log10(pval))) + geom_point() + 
  facet_wrap(~ chr, scale = "free_x") + 
  geom_hline(yintercept = -log10(5e-8), linetype = 2, color = "red") + 
  labs(title = title) +
  theme_void()

title <- paste0("Vitamin D, ", vitd, ", n = ", assoc_vitd$n_complete_samples[1])  
p_vitd <- ggplot(assoc_fat, aes(pos, -log10(pval))) + geom_point() + 
  facet_wrap(~ chr, scale = "free_x") + 
  geom_hline(yintercept = -log10(5e-8), linetype = 2, color = "red") + 
  labs(title = title) +
  theme_void()
  
