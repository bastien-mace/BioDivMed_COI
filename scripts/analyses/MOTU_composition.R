# ---- 0. Load packages ----
library(tidyverse)
library(RColorBrewer)
# __________________________



# ---- 1. Import data ----
## 1-1. Import MOTU table
MOTU_table_final <- read.csv('../../data/MOTU_table_final.csv') %>% 
  select(-c(pid, ltg_taxid, ltg_name, ltg_rank, ltg_rank_index, domain, class, order, family, genus, species, COUNT, seq, NIS))


## 1-2. Import extended metadata
metadata <- read.csv2('../../data/site_metadata.csv') %>%
  pivot_longer(c(spygen_code_1, spygen_code_2),
               values_to = 'sample') %>% 
  select(pool, Name, Habitat, sample)
metadata$Name <- gsub('Saint-Cyprien', 'Saint-Cyprien\n(West)', metadata$Name)
metadata$Name <- gsub('Agde', 'Agde\n(Central)', metadata$Name)
metadata$Name <- gsub('Le Grau-du-Roi', 'Le Grau-du-Roi\n(East)', metadata$Name)


## 1-3. Display kingdoms and phyla detected in the dataset
kingdoms_phyla <- MOTU_table_final %>%
  select(kingdom, phylum) %>%
  unique() %>%
  arrange(kingdom, phylum)
kingdoms_phyla

### NOTE: Some MOTUs are completely unassigned
###       Some Metazoa MOTUs have unassigned phyla
###       Fungi MOTUs have unassigned phyla

### Transform assignations in MOTU table
MOTU_table_kingdoms_phyla <- MOTU_table_final %>% 
  mutate(kingdom = ifelse(is.na(kingdom), 'Unassigned Eukaryota', kingdom),
         phylum = ifelse(is.na(phylum) & kingdom == 'Unassigned Eukaryota', 'Unassigned Eukaryota', phylum),
         phylum = ifelse(is.na(phylum) & kingdom == 'Metazoa', 'Unassigned Metazoa', phylum),
         phylum = ifelse(is.na(phylum) & kingdom == 'Fungi', 'Unassigned Fungi', phylum))
# ________________________



# ---- 2. Phyla composition of dominant MOTUs per site ----
## 2-1. Calculate proportions of dominant phyla in every site
dominant_MOTU_composition_table <- MOTU_table_kingdoms_phyla %>%
  mutate(across(starts_with('SPY'), ~ as.integer(.x > 0))) %>% # transform abundance into presence/absence
  pivot_longer(cols = starts_with('SPY'),
               names_to = 'sample',
               values_to = 'presence') %>%
  filter(presence == 1) %>%
  left_join(metadata %>%
              distinct(sample, pool), by = 'sample') %>%
  count(pool, phylum, name = 'n_MOTU') %>% # calculate the number of MOTUs assigned to each phylum in every site
  group_by(pool) %>%
  mutate(prop = n_MOTU / sum(n_MOTU)) %>% # calculate the proportions of each phylum in every site
  ungroup() %>%
  group_by(phylum) %>%
  mutate(mean_prop = mean(prop), # calculate the mean proportion of phylum across sites
         dominant_phylum = if_else(mean_prop < 0.01, 'Other Eukaryota', phylum)) %>% # only consider dominant MOTUs (mean proportion > 1%)
  ungroup() %>%
  group_by(pool, dominant_phylum) %>%
  summarise(prop = sum(prop), .groups = 'drop') %>% # sum proportions to have one 'Other Eukaryotes' per site
  left_join(metadata %>% distinct(pool, Name, Habitat), by = 'pool') %>%
  mutate(Habitat = factor(Habitat,
                          levels = c('Seaport', 'Rock', 'Sand')),
         Name = factor(Name,
                       levels = c('Saint-Cyprien\n(West)', 'Agde\n(Central)', 'Le Grau-du-Roi\n(East)')),
         pool = factor(pool, levels = unique(pool[order(Name, Habitat)])))


## 2-2. Calculate the total number of MOTUs per site
total_MOTU_table <- MOTU_table_kingdoms_phyla %>%
  mutate(across(starts_with('SPY'), ~ as.integer(.x > 0))) %>%
  pivot_longer(cols = starts_with('SPY'),
               names_to = 'sample',
               values_to = 'presence') %>%
  filter(presence == 1) %>%
  left_join(metadata %>%
              distinct(sample, pool, Name), by = 'sample') %>%
  distinct(pool, Name, id) %>%
  count(pool, Name, name = 'n_MOTU') %>%
  mutate(pool = factor(pool, levels = levels(dominant_MOTU_composition_table$pool)),
         Name = factor(Name, levels = levels(dominant_MOTU_composition_table$Name)))


## 2-3. Display stacked barplot
### Order dominant phyla
unique(dominant_MOTU_composition_table$dominant_phylum)
dominant_MOTU_composition_table$dominant_phylum <- factor(dominant_MOTU_composition_table$dominant_phylum,
                                                          levels = c('Annelida', 'Arthropoda', 'Chordata', 'Cnidaria', 'Mollusca', 'Nemertea', 'Unassigned Metazoa',
                                                                     'Bacillariophyta', 'Ochrophyta',  'Oomycota', 
                                                                     'Chlorophyta', 'Rhodophyta',
                                                                     'Unassigned Eukaryota', 'Other Eukaryota'))

### Associated colors
pal <- c('#654522', '#e25822', '#604e97', '#e68fac', '#882d17', '#f38400', 'grey90',
         '#008856', '#b3446c', '#0067a5',
         '#8db600', '#be0032',
         'grey70', 'grey30')
names(pal) <- levels(dominant_MOTU_composition_table$dominant_phylum)

### Plot
dominant_MOTUs_plot <- ggplot(dominant_MOTU_composition_table,
                              aes(pool, prop, fill = dominant_phylum)) +
  geom_col(color = 'black', width = 0.85) +
  scale_fill_manual(values = pal) +
  scale_x_discrete(labels = levels(dominant_MOTU_composition_table$Habitat), 
                   position = 'bottom') +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.02))) +
  facet_wrap(~ Name, scales = 'free_x', strip.position = 'top') +
  geom_text(data = total_MOTU_table,
            aes(x = pool, y = 1.03, label = n_MOTU),
            inherit.aes = FALSE, size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(face = 'bold', colour = 'black', size = 7.5),
        axis.text.y = element_text(face = 'bold', colour = 'black'),
        axis.title = element_blank(),
        legend.text = element_text(colour = 'black'),
        legend.title = element_blank(),
        strip.text = element_text(face = 'bold', colour = 'black', size = 10),
        strip.background = element_blank(),
        panel.spacing = unit(1, 'lines'),
        panel.grid.major.x  = element_blank(),
        panel.grid.major.y = element_line(colour = 'black'),
        panel.grid.minor.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(colour = 'black'))
dominant_MOTUs_plot

### Save in pdf
ggsave(dominant_MOTUs_plot,
       file = '../../outputs/dominant_MOTUs_plot.pdf',
       height = 10,
       width = 20,
       unit = 'cm')
# _________________________________________________________



# ---- 3. Phyla composition of non-dominant MOTUs per site ----
## 3-1. Calculate proportions of non-dominant phyla in every site
other_MOTU_composition_table <- MOTU_table_kingdoms_phyla %>%
  mutate(across(starts_with('SPY'), ~ as.integer(.x > 0))) %>% # transform abundance into presence/absence
  pivot_longer(cols = starts_with('SPY'),
               names_to = 'sample',
               values_to = 'presence') %>%
  filter(presence == 1) %>%
  left_join(metadata %>%
              distinct(sample, pool), by = 'sample') %>%
  filter(!phylum %in% levels(dominant_MOTU_composition_table$dominant_phylum)) %>%  # only consider non-dominant MOTUs
  count(pool, phylum, name = 'n_MOTU') %>% # calculate the number of MOTUs assigned to each non-dominant phylum in every site
  group_by(pool) %>%
  mutate(prop = n_MOTU / sum(n_MOTU)) %>% # calculate the proportions of each non-dominant phylum in every site
  ungroup() %>%
  left_join(metadata %>% distinct(pool, Name, Habitat), by = 'pool') %>%
  mutate(Habitat = factor(Habitat,
                          levels = c('Seaport', 'Rock', 'Sand')),
         Name = factor(Name,
                       levels = c('Saint-Cyprien\n(West)', 'Agde\n(Central)', 'Le Grau-du-Roi\n(East)')),
         pool = factor(pool, levels = unique(pool[order(Name, Habitat)])))


## 3-2. Calculate the total number of MOTUs assigned to non-dominant phyla per site
total_other_MOTU_table <- MOTU_table_kingdoms_phyla %>%
  mutate(across(starts_with('SPY'), ~ as.integer(.x > 0))) %>%
  pivot_longer(cols = starts_with('SPY'),
               names_to = 'sample',
               values_to = 'presence') %>%
  filter(presence == 1) %>%
  left_join(metadata %>%
              distinct(sample, pool, Name), by = 'sample') %>%
  filter(!phylum %in% levels(dominant_MOTU_composition_table$dominant_phylum)) %>%  # only consider non-dominant MOTUs
  distinct(pool, Name, id) %>%
  count(pool, Name, name = 'n_MOTU') %>%
  mutate(pool = factor(pool, levels = levels(dominant_MOTU_composition_table$pool)),
         Name = factor(Name, levels = levels(dominant_MOTU_composition_table$Name)))


## 3-3. Display stacked barplot
### Order dominant phyla
unique(other_MOTU_composition_table$phylum)
other_MOTU_composition_table$phylum <- factor(other_MOTU_composition_table$phylum,
                                              levels = c('Bryozoa', 'Chaetognatha', 'Echinodermata', 'Phoronida', 'Platyhelminthes', 'Porifera', 'Rotifera',
                                                         'Cercozoa', 'Cryptophyta', 'Haptophyta', 'Myzozoa',
                                                         'Amoebozoa',
                                                         'Unassigned Fungi'))

### Associated colors
pal <- c(brewer.pal(12, 'Set3'), 'white')
names(pal) <- levels(other_MOTU_composition_table$phylum)

### Plot
other_MOTUs_plot <- ggplot(other_MOTU_composition_table,
                           aes(pool, prop, fill = phylum)) +
  geom_col(color = 'black', width = 0.55) +
  scale_fill_manual(values = pal) +
  scale_x_discrete(labels = levels(other_MOTU_composition_table$Habitat),
                   position = 'top') +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0.02, 0))) +
  facet_wrap(~ Name, scales = 'free_x', strip.position = 'bottom') +
  geom_text(data = total_other_MOTU_table,
            aes(x = pool, y = -0.03, label = n_MOTU),
            inherit.aes = FALSE, size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(face = 'bold', colour = 'black', size = 7.5),
        axis.text.y = element_text(face = 'bold', colour = 'black'),
        axis.title = element_blank(),
        legend.text = element_text(colour = 'black'),
        legend.title = element_blank(),
        strip.text = element_text(face = 'bold', colour = 'black', size = 10),
        strip.background = element_rect(color = NA),
        panel.spacing = unit(1, 'lines'),
        panel.grid.major.x  = element_blank(),
        panel.grid.major.y = element_line(colour = 'black'),
        panel.grid.minor.y  = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_line(colour = 'black'))
other_MOTUs_plot

### Save in pdf
ggsave(other_MOTUs_plot,
       file = '../../outputs/other_MOTUs_plot.pdf',
       height = 9,
       width = 20,
       unit = 'cm')
# _____________________________________________________________