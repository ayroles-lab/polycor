### inc
library(MASS)

library(Matrix)
library(lme4)
library(lme4qtl)

library(magrittr)
library(tidyverse)

path <- switch(Sys.info()[['nodename']],
  "~/git/hemostat/polycor/" )
  
file.path(path, "relmer.R") %>% source

### data
data(raz2005, package = "rmcorr")

dat <- as_tibble(raz2005) %>% 
  mutate(Obs = as.factor(1:n()), Obs2 = Obs, Participant2 = Participant,
  Volume = as.numeric(scale(Volume)),
  Age = as.numeric(scale(Age)))

### explore data
p1 <- ggplot(dat, aes(Age, Volume, group = Participant)) + geom_line()

gr <- group_by(dat, Participant) %>% 
  summarize(Slope = diff(Volume) / diff(Age)) %>% 
  ungroup %>% arrange(Slope)


dat <- left_join(dat,
  select(gr, Participant, Slope), by = "Participant") %>%
  mutate(Group = cut(Slope, breaks = c(-3.1, 0, 2.2), labels = c("Neg", "Pos")),
    Group = droplevels(Group))
    
p2 <- ggplot(dat, aes(Age, Volume, group = Participant, color = Group)) + geom_line()

### random-slope model
m1 <- relmer(Volume ~ Age + (1 | Participant), data = dat, REML = FALSE)

m2 <- relmer(Volume ~ Age + (1 | Participant) + (0 + Age | Participant), data = dat, REML = FALSE)

m3 <- relmer(Volume ~ Age + (1 | Participant) + (0 + Age | Group), data = dat, REML = FALSE)

### bi-variate model
bdat <- gather(dat, tname, tvalue, Age, Volume) %>%
  mutate(tname_gr = paste(tname, Group, sep = "_"))

b1 <- relmer(tvalue ~ tname - 1 + (0 + tname | Participant) + (0 + tname | Obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

b2 <- relmer(tvalue ~ tname - 1 + (0 + tname_gr | Participant) + (0 + tname_gr | Obs), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

b3 <- relmer(tvalue ~ tname - 1 + (0 + dummy(tname_gr, c("Age_Neg", "Volume_Neg")) | Participant) + (0 + dummy(tname_gr, c("Age_Pos", "Volume_Pos")) | Participant2) + (0 + dummy(tname_gr, c("Age_Neg", "Volume_Neg")) | Obs) + (0 + dummy(tname_gr, c("Age_Pos", "Volume_Pos")) | Obs2), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

b4 <- relmer(tvalue ~ tname - 1 + (0 + tname | Participant) + (0 + dummy(tname_gr, c("Age_Neg", "Volume_Neg")) | Obs) + (0 + dummy(tname_gr, c("Age_Pos", "Volume_Pos")) | Obs2), data = bdat, weights = rep(1e10, nrow(bdat)), calc.derivs = FALSE)

### anova on bi-variate models
anova(b2, b3) # no difference

anova(b1, b2) # significant difference

anova(b4, b2) # shouldn't be any diffrence: Group was constructed based on Slope rather than Baseline
