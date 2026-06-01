# ---- 0. Load packages ----
library(tidyverse)
library(nlme) # linear mixed-models
library(corMLPE) # create correlation structure for MLPE models
library(parameters) # compute standardized effect sizes
library(wesanderson)
# __________________________



# ---- 1. Import data ----
## 1-1. MOTU table
MOTU_table_subset <- read.csv('../../data/MOTU_table_subset.csv')


## 1-2. Load Jost's D distance matrices
load('../../data/beta_GD.Rdata')
MOTU_D_matrices <- D$MOTU_D_matrices


## 1-3. Site metadata
metadata <- read.csv2('../../data/site_metadata.csv') %>%
  pivot_longer(cols = c(spygen_code_1, spygen_code_2),
               values_to = 'sample')

### Pool the problematic sample
samples_to_pool <- c('SPY232133', 'SPY232146')
pooled_name <- 'SPY232133_SPY232146'

### Modify metadata
metadata_new <- metadata[!metadata$sample %in% samples_to_pool, ]
pooled_meta <- metadata[metadata$sample == samples_to_pool[1], ]
pooled_meta$sample <- pooled_name
metadata_new <- rbind(metadata_new, pooled_meta) %>% 
  mutate(sample = factor(sample, levels = rownames(as.matrix(MOTU_D_matrices[[1]])))) %>%
  arrange(sample)


## 1-4. Geographic distance matrix
### Load the matrix of marine geographic distances
load('../../data/geo_dist.Rdata')

### Convert into long format
geo_dist_long <- as.data.frame(as.matrix(geo_dist))
geo_dist_long$sample_A <- rownames(geo_dist_long)
geo_dist_long <- pivot_longer(geo_dist_long, 
                              cols = -sample_A, 
                              names_to = 'sample_B', 
                              values_to = 'geo_dist')
# ________________________



# ---- 2. Compute maximum-likelihood population-effects (MLPE) linear mixed-effects models ----
## 2-1. Create an empty dataframe
MLPE_results_df <- data.frame(MOTU = as.character(),
                              perc = as.numeric(),
                              Habitat_beta = as.numeric(),
                              Habitat_p = as.numeric(),
                              Geo_beta = as.numeric(),
                              Geo_p = as.numeric())


## 2-2. Loop across all MOTUs
for (i in 1:length(MOTU_D_matrices)){
  ### Retrieve MOTU Jost's D matrix
  jost_dist <- MOTU_D_matrices[[i]]
  
  ### Convert into long format
  jost_dist_long <- as.data.frame(as.matrix(jost_dist))
  jost_dist_long$sample_A <- rownames(jost_dist_long)
  jost_dist_long <- pivot_longer(jost_dist_long, 
                                 cols = -sample_A, 
                                 names_to = 'sample_B', 
                                 values_to = 'jost')
  
  ### Add geographic distance and habitat data
  data_long <- left_join(jost_dist_long, geo_dist_long, by = c('sample_A', 'sample_B')) %>% 
    left_join(metadata_new %>% select(sample, Habitat, pool), 
              by = c('sample_A' = 'sample')) %>%
    rename(Habitat_A = Habitat,
           site_A = pool) %>%
    left_join(metadata_new %>% select(sample, Habitat, pool), 
              by = c('sample_B' = 'sample')) %>%
    rename(Habitat_B = Habitat,
           site_B = pool) %>%
    mutate(habitat_num_A = ifelse(Habitat_A == 'Seaport', 1, 0),
           habitat_num_B = ifelse(Habitat_B == 'Seaport', 1, 0),
           habitat_dist = abs(habitat_num_A - habitat_num_B)) %>% # Distance of 1 for Seaport vs Sand/Rock, 0 otherwise
    filter(site_A != site_B) %>%
    drop_na()
  
  ### Run MLPE (skip errors)
  MLPE_model <- try(lme(jost ~ geo_dist + habitat_dist,
                        random = ~ 1 | site_A,
                        correlation = corMLPE(form = ~ sample_A + sample_B),
                        data = data_long,
                        method = 'ML'),
                    silent = TRUE)
  if(inherits(MLPE_model, 'try-error')) next
  
  ### Record model estimates
  summary <- summary(MLPE_model)
  
  ### Compute standardized effect sizes
  betas <- standardize_parameters(MLPE_model)
  
  ### Record the results
  MLPE_results_df <- rbind(MLPE_results_df,
                           data.frame(MOTU = names(MOTU_D_matrices[i]),
                                      perc = sum(!is.na(MOTU_D_matrices[[i]]))/136,
                                      Habitat_beta = betas$Std_Coefficient[3],
                                      Habitat_p = summary$tTable[3, 5],
                                      Geo_beta = betas$Std_Coefficient[2],
                                      Geo_p = summary$tTable[2, 5]))
  
}


## 2-3. Compute MLPE on mean Jost's D matrix
### Import the average Jost's D matrix
mean_D_matrix <- D$mean_D_matrix

### Convert into long format
jost_dist_long <- as.data.frame(as.matrix(mean_D_matrix))
jost_dist_long$sample_A <- rownames(jost_dist_long)
jost_dist_long <- pivot_longer(jost_dist_long, 
                               cols = -sample_A, 
                               names_to = 'sample_B', 
                               values_to = 'jost')

### Add geographic distance and habitat data
data_long <- left_join(jost_dist_long, geo_dist_long, by = c('sample_A', 'sample_B')) %>% 
  left_join(metadata_new %>% select(sample, Habitat, pool), 
            by = c('sample_A' = 'sample')) %>%
  rename(Habitat_A = Habitat,
         site_A = pool) %>%
  left_join(metadata_new %>% select(sample, Habitat, pool), 
            by = c('sample_B' = 'sample')) %>%
  rename(Habitat_B = Habitat,
         site_B = pool) %>%
  mutate(habitat_num_A = ifelse(Habitat_A == 'Seaport', 1, 0),
         habitat_num_B = ifelse(Habitat_B == 'Seaport', 1, 0),
         habitat_dist = abs(habitat_num_A - habitat_num_B)) %>% # Distance of 1 for Seaport vs Sand/Rock, 0 elsewhere
  filter(site_A != site_B) %>%
  drop_na()

### Run MLPE
MLPE_model <- lme(jost ~ geo_dist + habitat_dist,
                  random = ~ 1 | site_A,
                  correlation = corMLPE(form = ~ sample_A + sample_B),
                  data = data_long,
                  method = 'ML')

### Record model estimates
summary <- summary(MLPE_model)

### Compute standardized effect sizes
betas <- standardize_parameters(MLPE_model)

### Merge the results in a matrix
MLPE_model_mean_df <- data.frame(Habitat_beta = betas$Std_Coefficient[3],
                                 Habitat_p = summary$tTable[3, 5],
                                 Geo_beta = betas$Std_Coefficient[2],
                                 Geo_p = summary$tTable[2, 5])
# _____________________________________________________________________________________________



# ---- 3. Plot the Habitat/Geography standardized effect sizes ----
## 3-1. Transform the results into a long dataframe
MLPE_results_long_df <- MLPE_results_df %>%
  filter(perc > 0.1) %>% # remove MOTUs with very few comparisons
  drop_na() %>% 
  left_join(MOTU_table_subset %>% 
              select(id, ltg_name, ltg_rank, kingdom),
            by = c('MOTU' = 'id')) %>% 
  mutate(kingdom = ifelse(is.na(kingdom), 'Unassigned Eukaryota', kingdom),
         Habitat_p_adj = p.adjust(Habitat_p, method = 'fdr'), # adjust p-values
         Geo_p_adj = p.adjust(Geo_p, method = 'fdr'), # adjust p-values
         signif = case_when(Habitat_p_adj < 0.05 & Geo_p_adj < 0.05 ~ 'Habitat + Geography',
                            Habitat_p_adj < 0.05 & Geo_p_adj >= 0.05 ~ 'Habitat only',
                            Habitat_p_adj >= 0.05 & Geo_p_adj < 0.05 ~ 'Geography only',
                            TRUE ~ 'None')) %>%
  mutate(signif = factor(signif, levels = c('Habitat + Geography', 'Habitat only', 'Geography only', 'None')))

## 3-2. Plot Habitat/Geography standardized effect sizes scatterplot
INTRA_MOTU_BETA_MLPE_plot <- ggplot() +
  geom_point(data = MLPE_model_mean_df,
            aes(x = Habitat_beta, y = Geo_beta),
            size = 5, shape = 21, color = 'black', fill = '#ff0000') +
  geom_text(data = MLPE_model_mean_df,
           aes(x = Habitat_beta, y = Geo_beta,
               label = paste0("Mean Jost's D\n(", round(Habitat_beta, 2), ', ', round(Geo_beta, 2), ')')),
           size = 4, vjust = -0.5, fontface = 'bold') +
  geom_point(data = MLPE_results_long_df,
             aes(x = Habitat_beta, y = Geo_beta,
                 fill = signif),
             size = 2, shape = 21, color = 'black') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey50') +
  geom_abline(slope = 1, intercept = 0, linetype = 'dotted', color = 'black') +
  scale_fill_manual(values = c(wes_palette('Darjeeling1')[1:3], 'grey90')) +
  guides(fill = guide_legend(title.position = 'top',
                             title.hjust = 0.5, nrow = 1)) +
  labs(x = 'Habitat effect (SES)',
       y = 'Geography effect (SES)',
       fill = 'Significance (FDR corrected)') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15),
        legend.background = element_rect(color = 'black', fill = 'white'),
        legend.position = 'bottom',
        legend.title = element_text(face = 'bold'))
INTRA_MOTU_BETA_MLPE_plot

### Save in pdf
ggsave(INTRA_MOTU_BETA_MLPE_plot,
       file = '../../outputs/INTRA_MOTU_BETA_MLPE_plot.pdf',
       height = 15,
       width = 15,
       unit = 'cm')
# _________________________________________________________________



# ---- 4. Plot the proportion of MOTUs driven by each component ----
## 4-1. Create a dataframe of kingdom proportions
kingdom_proportions_df <- MLPE_results_long_df %>%
  group_by(kingdom, signif) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(kingdom) %>%
  mutate(prop = n / sum(n)) %>% 
  ungroup() %>% 
  mutate(kingdom = factor(kingdom, levels = c('Metazoa', 'Chromista', 'Plantae', 'Unassigned Eukaryota')),
         signif = factor(signif, levels = c('Geography only', 'Habitat + Geography', 'Habitat only', 'None')))


## 4-2. Count the total number of MOTUs in each kingdom
kingdom_total <- kingdom_proportions_df %>% 
  select(kingdom, n) %>% 
  group_by(kingdom) %>% 
  mutate(sum = sum(n)) %>% 
  ungroup() %>% 
  select(kingdom, sum) %>% 
  unique()


## 4-3. Plot kingdom proportions
INTRA_MOTU_BETA_MLPE_by_kingdoms_plot <- ggplot(kingdom_proportions_df,
                                                aes(x = kingdom, y = prop, fill = signif)) +
  geom_col(width = 0.5, color = 'black') +
  scale_y_continuous(labels = scales::percent_format()) +
  geom_text(data = kingdom_total,
            aes(x = kingdom, y = 1.03, label = sum),
           inherit.aes = FALSE, size = 3) +
  scale_fill_manual(values = c(wes_palette('Darjeeling1')[c(3, 1, 2)], 'grey90')) +
  guides(fill = guide_legend(title.position = 'top',
                             title.hjust = 0.5, nrow = 1)) +
  labs(x = 'Kingdom',
       y = 'Proportion',
       fill = 'Significance (FDR corrected)') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15),
        legend.background = element_rect(color = 'black', fill = 'white'),
        legend.position = 'bottom',
        legend.title = element_text(face = 'bold'))
INTRA_MOTU_BETA_MLPE_by_kingdoms_plot

### Save in pdf
ggsave(INTRA_MOTU_BETA_MLPE_by_kingdoms_plot,
       file = '../../outputs/INTRA_MOTU_BETA_MLPE_by_kingdoms_plot.pdf',
       height = 15,
       width = 15,
       unit = 'cm')
# __________________________________________________________________