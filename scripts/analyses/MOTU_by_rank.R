# ---- 0. Load packages ----
library(tidyverse)
# __________________________



# ---- 1. Calculate the number of assigned MOTUs at each taxonomic rank ----
## 1-1. Import data
MOTU_table_final <- read.csv('../../data/MOTU_table_final.csv')


## 1-2. Taxonomic ranks
ranks <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')


## 1-3. Calculate on each rank
MOTU_by_rank <- MOTU_table_final %>%
  select(all_of(ranks)) %>%
  pivot_longer(cols = everything(),
               names_to = 'rank',
               values_to = 'taxon') %>%
  mutate(assigned = if_else(is.na(taxon) | taxon == '', 'Unassigned MOTUs', 'Assigned MOTUs')) %>%
  count(rank, assigned) %>%
  group_by(rank) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()
# __________________________________________________________________________



# ---- 2. Plot the proportions of MOTUs assigned ----
MOTU_by_rank_plot <- ggplot(MOTU_by_rank,
                            aes(x = rank, y = prop, fill = assigned)) +
  geom_col(position = 'stack') +
  geom_text(data = subset(MOTU_by_rank,
                          assigned == 'Assigned MOTUs'),
            aes(y = 0.95, label = n),
            size = 3) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = 'Taxonomic rank',
       y = '% of MOTUs',
       fill = NULL) +
  scale_x_discrete(limits = tax_ranks) +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10, angle = 45, hjust = 1),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        legend.position = 'bottom')
MOTU_by_rank_plot

#### Save in pdf
ggsave(MOTU_by_rank_plot,
       file = '../../outputs/MOTU_by_rank_plot.pdf',
       height = 8,
       width = 12,
       unit = 'cm')
# __________________________________________________________________________