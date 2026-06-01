# ---- 0. Load packages ----
library(tidyverse)
# __________________________



# ---- 1. Prepare data ----
## 1-1. Site metadata
metadata <- read.csv2('../../data/site_metadata.csv') %>%
  pivot_longer(cols = c(spygen_code_1, spygen_code_2),
               values_to = 'sample')


## 1-2. MOTU presence/absence matrix
### Restrict MOTU table to samples' columns
MOTU_table_final <- read.csv('../../data/MOTU_table_final.csv')
MOTU_pa_matrix <- MOTU_table_final %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

### Transform into presence-absence_matrix
MOTU_pa_matrix[MOTU_pa_matrix > 0] <- 1
colSums(MOTU_pa_matrix)

### Pool replicates
final_MOTU_matrix <- t(MOTU_pa_matrix) %>%  # transpose matrix to get samples in rows and MOTUs in columns
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  left_join(metadata %>% select(sample, pool, Name),
            by = 'sample')%>%
  select(-sample) %>%
  group_by(pool, Name) %>%
  summarise(across(where(is.numeric), max), .groups = 'drop') %>%
  column_to_rownames('pool') %>%
  select(-Name)


## 1-3. ASV presence/absence matrix
### Restrict ASV table to samples' columns
ASV_table_final <- read.csv('../../data/ASV_table_final.csv')
ASV_pa_matrix <- ASV_table_final %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

### Transform into presence-absence_matrix
ASV_pa_matrix[ASV_pa_matrix > 0] <- 1
colSums(ASV_pa_matrix)

### Pool replicates
final_ASV_matrix <- t(ASV_pa_matrix) %>%  # transpose matrix to get samples in rows and MOTUs in columns
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  left_join(metadata %>% select(sample, pool, Name),
            by = 'sample')%>%
  select(-sample) %>%
  group_by(pool, Name) %>%
  summarise(across(where(is.numeric), max), .groups = 'drop') %>%
  column_to_rownames('pool') %>%
  select(-Name)


## 1-4. Total read count
read_count_matrix <- MOTU_table_final %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

### Pool replicates
tot_read_count <- t(read_count_matrix) %>%  # transpose matrix to get samples in rows and MOTUs in columns
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  left_join(metadata %>% select(sample, pool, Name),
            by = 'sample')%>%
  select(-sample) %>%
  group_by(pool, Name) %>%
  summarise(across(where(is.numeric), sum), .groups = 'drop') %>%
  column_to_rownames('pool') %>%
  select(-Name) %>% 
  rowSums() %>% 
  unlist()


## 1-5. Reorder metadata
metadata_final <- metadata %>% 
  select(pool, Habitat, Area) %>% 
  unique() %>% 
  mutate(pool = factor(pool, levels = rownames(final_MOTU_matrix))) %>%
  arrange(pool)
# _________________________



# ---- 2. Calculate MOTU and ASV richness ----
## 2-1. MOTU richness
MOTU_richness_df <- data.frame(pool = rownames(final_MOTU_matrix),
                               Richness = rowSums(final_MOTU_matrix),
                               Count = tot_read_count) %>%
  inner_join(metadata_final) %>%
  select(Richness, pool, Count, Habitat, Area) %>% 
  mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
         Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))


## 2-2. ASV richness
ASV_richness_df <- data.frame(pool = rownames(final_ASV_matrix),
                              Richness = rowSums(final_ASV_matrix),
                              Count = tot_read_count) %>%
  inner_join(metadata_final) %>%
  select(Richness, pool, Count, Habitat, Area) %>% 
  mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
         Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))
# ____________________________________________



# ---- 3. Compute correlations ----
## 3-1. MOTU richness vs Total read count
cor.test(MOTU_richness_df$Count, MOTU_richness_df$Richness)

### NOTE: Negative correlation significant ***
###       => No sampling bias
###       => Sites with lower MOTU richness may host higher biomass (e.g., seaports)

MOTU_richness_vs_Total_read_count_plot <- ggplot(MOTU_richness_df,
                                                 aes(x = Count,
                                                     y = Richness,
                                                     color = Habitat)) +
  geom_point(size = 2) +
  scale_color_manual(values = wes_palette('Moonrise2')) +
  guides(color = guide_legend(title.position = 'top',
                              title.hjust = 0.5, nrow = 1)) +
  xlab('Total read count') + ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10, angle = 45, vjust = 0.5),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15),
        legend.background = element_rect(color = 'black', fill = 'white'),
        legend.position = 'top',
        legend.title.align = 0.5,
        legend.title = element_text(face = 'bold'))
MOTU_richness_vs_Total_read_count_plot

### Save in pdf
ggsave(MOTU_richness_vs_Total_read_count_plot,
       file = '../../outputs/MOTU_richness_vs_Total_read_count_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 3-2. ASV richness vs Total read count
cor.test(ASV_richness_df$Count, ASV_richness_df$Richness)

### NOTE: Negative correlation not significant
###       => No sampling bias

ASV_richness_vs_Total_read_count_plot <- ggplot(ASV_richness_df,
                                                aes(x = Count,
                                                    y = Richness,
                                                    color = Habitat)) +
  geom_point(size = 2) +
  scale_color_manual(values = wes_palette('Moonrise2')) +
  guides(color = guide_legend(title.position = 'top',
                              title.hjust = 0.5, nrow = 1)) +
  xlab('Total read count') + ylab('ASV richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10, angle = 45, vjust = 0.5),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15),
        legend.background = element_rect(color = 'black', fill = 'white'),
        legend.position = 'top',
        legend.title.align = 0.5,
        legend.title = element_text(face = 'bold'))
ASV_richness_vs_Total_read_count_plot

### Save in pdf
ggsave(ASV_richness_vs_Total_read_count_plot,
       file = '../../outputs/ASV_richness_vs_Total_read_count_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

## 3-3. MOTU richness vs ASV_richness
ASV_MOTU_richness_df <- inner_join(ASV_richness_df %>%
                                     rename(ASV_richness = Richness),
                                   MOTU_richness_df %>%
                                     rename(MOTU_richness = Richness))
cor.test(ASV_MOTU_richness_df$MOTU_richness, ASV_MOTU_richness_df$ASV_richness)

### NOTE: Positive correlation significant **
###       => Expected pattern: ASV richness increases with MOTU richness

MOTU_richness_vs_ASV_richness_plot <- ggplot(ASV_MOTU_richness_df,
                                             aes(x = MOTU_richness,
                                                 y = ASV_richness,
                                                 color = Habitat)) +
  geom_point(size = 2) +
  scale_color_manual(values = wes_palette('Moonrise2')) +
  guides(color = guide_legend(title.position = 'top',
                              title.hjust = 0.5, nrow = 1)) +
  xlab('MOTU richness') + ylab('ASV richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15),
        legend.background = element_rect(color = 'black', fill = 'white'),
        legend.position = 'top',
        legend.title.align = 0.5,
        legend.title = element_text(face = 'bold'))
MOTU_richness_vs_ASV_richness_plot

### Save in pdf
ggsave(MOTU_richness_vs_ASV_richness_plot,
       file = '../../outputs/MOTU_richness_vs_ASV_richness_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 3-4. Number of reads per MOTU vs Number of ASVs per MOTU across all samples
reads_per_MOTU <- MOTU_table_final %>% 
  select(id, COUNT)
ASVs_per_MOTU <- ASV_table_final %>% 
  count(MOTU, name = 'n_ASV')
reads_and_ASVs_per_MOTU_df <- inner_join(reads_per_MOTU, ASVs_per_MOTU, by = c('id' = 'MOTU'))
cor.test(log(reads_and_ASVs_per_MOTU_df$COUNT), log(reads_and_ASVs_per_MOTU_df$n_ASV), method = 'spearman')

### NOTE: Positive correlation significant ***
###       => MOTU ASV richness increases with MOTU read count
###       => Abundant MOTUs are more genetically diverse
###       => Could be interpreted like some ASVs are artifacts when MOTU has a high read count, but the correlation was negative and not significant between read count and ASV richness

reads_and_ASVs_per_MOTU_plot <- ggplot(reads_and_ASVs_per_MOTU_df,
                                       aes(x = log(COUNT),
                                           y = log(n_ASV))) +
  geom_point(size = 2) +
  xlab('log(MOTU read count)') + ylab('log(Number of ASVs per MOTU)') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15))
reads_and_ASVs_per_MOTU_plot

### Save in pdf
ggsave(reads_and_ASVs_per_MOTU_plot,
       file = '../../outputs/reads_and_ASVs_per_MOTU_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')
# _________________________________