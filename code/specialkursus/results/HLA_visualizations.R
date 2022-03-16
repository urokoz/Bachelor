library(tidyverse)

df <- read.table("../results/donor/_raw/ragweed_HLA_count.txt", 
                 sep="\t")
as_tibble(df)

ggplot(df, mapping = aes(x)) + 
  geom_bar()