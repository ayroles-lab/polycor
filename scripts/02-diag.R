library(ggdag)
library(tidyverse)

coords <- tribble(
  ~name,      ~x,  ~y,
  "t1",       2,   3,
  "t2",       3,   3,
  "G1_slope", 2.5, 4,
  "holder",   2.5, 2.9,
  "G1_add",   1.5,   1.5,
  "G2_add",   3.5,   1.5,
  "G12_pleio", 2.5,   1.5,
)

dagify(
  t2 ~ t1,
  holder ~ G1_slope,
  t1 ~ G1_add, t2 ~ G2_add,
  t1 ~ G12_pleio, t2 ~ G12_pleio,
  coords = coords
) %>% 
  ggplot(aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_dag_point(data = function(x) filter(x, name != "holder"), col = "white") +
    geom_dag_edges() + 
    geom_dag_text(data = function(x) filter(x, name != "holder"), col = "black", size = 5) +
    theme_dag()
    
# save
# ggsave("scripts/figures/01-diag-G1-slope-and-pleio.png")
