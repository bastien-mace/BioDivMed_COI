# ---- 0. Load packages ----
library(tidyverse)
# __________________________


# ---- 1. Calculate the number of ASVs per MOTU ----
## 1-1. Import data
ASV_table_final <- read.csv('../../data/ASV_table_final.csv')


## 1-2. Calculate ASV richness for each MOTU
ASV_per_MOTU <- ASV_table_final %>%
  count(MOTU, name = 'n_ASV')


## 1-3. Separate in different classes
ASV_per_MOTU <- ASV_per_MOTU %>%
  mutate(class = case_when(n_ASV < 16  ~ as.character(n_ASV),
                           n_ASV > 15 & n_ASV < 21 ~ '[16–20]',
                           n_ASV > 20 & n_ASV < 26 ~ '[21–25]',
                           n_ASV > 25 & n_ASV < 31 ~ '[26–30]',
                           n_ASV > 30 & n_ASV < 41 ~ '[31–40]',
                           n_ASV > 40 & n_ASV < 51 ~ '[41–50]',
                           n_ASV > 50 & n_ASV < 101 ~ '[51–100]',
                           n_ASV > 100 ~ '> 100')) %>% 
  mutate(class = factor(class, levels = c(as.character(1:15), '[16–20]', '[21–25]', '[26–30]', '[31–40]', '[41–50]', '[51–100]', '> 100')))


## 1-4. Count the number of MOTU in each class
class_summary <- ASV_per_MOTU %>%
  count(class, name = 'n_MOTU')
# __________________________________________________



# ---- 2. Plot the Number of MOTUs vs ASV richness ----
ASV_richness_per_MOTU_plot <- ggplot(class_summary,
                                     aes(x = class, y = n_MOTU)) +
  geom_col(width = 0.8) +
  xlab('ASV richness') + ylab('Number of MOTUs') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10, vjust = c(rep(1, 15), 1, 0.8, 1, 0.8, 1, 0.8, 1)),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15))
ASV_richness_per_MOTU_plot

#### Save in pdf
ggsave(ASV_richness_per_MOTU_plot,
       file = '../../outputs/ASV_richness_per_MOTU_plot.pdf',
       height = 10,
       width = 20,
       unit = 'cm')
# _____________________________________________________