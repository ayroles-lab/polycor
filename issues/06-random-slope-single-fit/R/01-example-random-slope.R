### inc
library(lme4)
library(lme4qtl)

### data
data(dat40)
dat40 <- within(dat40, ID2 <- ID)

m0 <- relmatLmer(trait1 ~ AGE + (1|ID), dat40, relmat = list(ID = kin2))

m1 <- lmer(trait1 ~ AGE + (1 + AGE||FAMID), dat40)
m2 <- relmatLmer(trait1 ~ AGE + (1 + AGE||ID), dat40, relmat = list(ID = kin2))

m3 <- relmatLmer(trait1 ~ AGE + (1|ID) + (0 + AGE|ID), dat40, relmat = list(ID = kin2))
m4 <- relmatLmer(trait1 ~ AGE + (1|ID) + (0 + AGE|ID2), dat40, relmat = list(ID = kin2, ID2 = kin2))

m5 <- relmatLmer(trait1 ~ AGE + (0 + AGE|ID), dat40, relmat = list(ID = kin2))
