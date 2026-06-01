# ---- 0. Load packages ----
library(tidyverse)
library(wesanderson)
library(lme4) # linear mixed-models
library(car) # compute ANOVA
library(emmeans) # Tukey's HSD test
library(MuMIn) # obtain R²
library(ggsignif)
# __________________________



# ---- 1. Compute MOTU richness on all MOTUs ----
## 1-1. Import data
### 1-1-1. Presence/Absence matrix
#### Restrict MOTU table to samples' columns
MOTU_table_final <- read.csv('../../data/MOTU_table_final.csv')
pa_matrix <- MOTU_table_final %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

#### Transform into presence-absence_matrix
pa_matrix[pa_matrix > 0] <- 1
colSums(pa_matrix)
final_matrix <- t(pa_matrix) # transpose matrix to get samples in rows and MOTUs in columns

### 1-1-2. Site metadata
metadata <- read.csv2('../../data/site_metadata.csv') %>%
  pivot_longer(cols = c(spygen_code_1, spygen_code_2),
               values_to = 'sample') %>% 
  mutate(sample = factor(sample, levels = rownames(final_matrix))) %>%
  arrange(sample)

### 1-1-3. Pool the problematic sample
samples_to_pool <- c('SPY232133', 'SPY232146')
pooled_sample <- colSums(final_matrix[samples_to_pool, , drop = FALSE]) > 0
pooled_sample <- as.numeric(pooled_sample)
pooled_name <- 'SPY232133_SPY232146'

#### Remove the samples pooled and add the corresponding pool
final_matrix_new <- final_matrix[!rownames(final_matrix) %in% samples_to_pool, ]
final_matrix_new <- rbind(final_matrix_new, pooled_sample)
rownames(final_matrix_new)[nrow(final_matrix_new)] <- pooled_name

#### Modify metadata
metadata_new <- metadata[!metadata$sample %in% samples_to_pool, ]
pooled_meta <- metadata[metadata$sample == samples_to_pool[1], ]
pooled_meta$sample <- pooled_name
metadata_new <- rbind(metadata_new, pooled_meta)


## 1-2. Compute MOTU richness
richness_df <- data.frame(sample = rownames(final_matrix_new),
                            Richness = rowSums(final_matrix_new)) %>%
  inner_join(metadata_new) %>%
  select(Richness, pool, Habitat, Area) %>% 
  mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
         Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))


## 1-3. Data exploration
### 1-3-1. Descriptive statistics
str(richness_df)
summary(richness_df$Richness) # Y summary (Richness)
summary(richness_df$Habitat) # X summary (Habitat)
summary(richness_df$Area) # X summary (Area)

### 1-3-2. Distribution of Y
par(mfrow = c(2,2))

#### Boxplot
boxplot(richness_df$Richness, col = 'blue', ylab = 'Richness')

#### Cleveland plot
dotchart(richness_df$Richness, pch = 16, col = 'blue', xlab = 'Richness')

#### Histogram
hist(richness_df$Richness, col = 'blue', xlab = 'Richness', main = '')

#### Quantile-Quantile plot
qqnorm(richness_df$Richness, pch = 16, col = 'blue', xlab = '')
qqline(richness_df$Richness, col = 'red')

### NOTE: Y does not follow a Gaussian distribution (left-skewed)

dev.off()

### 1-3-3. Potential relationship Y/X
#### Habitat
boxplot(richness_df$Richness ~ richness_df$Habitat,
        varwidth = TRUE,
        ylab = 'MOTU richness',
        xlab = 'Habitat',
        col = wes_palette('Moonrise2')[1:3])
#### NOTE: Richness looks lower in seaports

#### Area
boxplot(richness_df$Richness ~ richness_df$Area,
        varwidth = TRUE,
        ylab = 'MOTU richness',
        xlab = 'Area',
        col = wes_palette('AsteroidCity1')[c(1, 3, 4)])
#### NOTE: No difference between areas is observed


## 1-4. Fit the model
LMM_model <- lmer(Richness ~ Area + Habitat + (1|pool),
                  data = richness_df)

### 1-4-1. Check normality of residuals
par(mfrow = c(1, 2))

#### Histogram
hist(residuals(LMM_model), col = 'blue', xlab = 'Residuals', main = 'Check Normality')

#### Quantile-Quantile plot
qqnorm(residuals(LMM_model), pch = 16, col = 'blue', xlab = '')
qqline(residuals(LMM_model), col = 'red')

#### Shapiro test
shapiro.test(residuals(LMM_model))
#### NOTE: Residuals follow normality

### 1-4-2. Check for homogeneity of variances (= Homoscedasticity)
par(mfrow = c(1, 3))

#### Residuals vs Fitted
plot(residuals(LMM_model) ~ fitted(LMM_model), col = 'blue', pch = 16, ylab = 'Residuals', xlab = 'Fitted')
abline(h = 0)

#### Residuals vs Area
boxplot(residuals(LMM_model) ~ richness_df$Area, varwidth = TRUE, ylab = 'Residuals', xlab = 'Area')
abline(h = 0)

#### Residuals vs Habitat
boxplot(residuals(LMM_model) ~ richness_df$Habitat, varwidth = TRUE, ylab = 'Residuals', xlab = 'Habitat')
abline(h = 0)
##### NOTE: There is no evident pattern in Residuals vs Fitted

### 1-4-3. Perform ANOVA
Anova(LMM_model)
#### NOTE: Effect of 'Area' is  not significant
#### NOTE: Effect of 'Habitat' is significant **

### 1-4-4. Perform post-hoc Tukey's HSD
emm <- emmeans(LMM_model, ~ Habitat)
pairs(emm, adjust = 'tukey')
#### NOTE: Difference in richness is significant * between 'Seaport' and 'Sand'
####       Difference in richness is significant ** between 'Seaport' and 'Rock'
####       Difference in richness is not significant between 'Rock' and 'Sand'

### 1-4-5. Obtain the R²
r.squaredGLMM(LMM_model)
#### NOTE: Marginal R² = 78.92%
####       Conditional R² = 88.27%


## 1-5. Plot MOTU richness vs Habitat
all_MOTU_ALPHA_Habitat_plot <- ggplot(richness_df,
                                      aes(x = Habitat,
                                          y = Richness,
                                          fill = Habitat)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  geom_signif(comparisons = list(c('Sand', 'Seaport'),
                                 c('Rock', 'Seaport')),
              annotations = c('*', '**'),
              y_position = c(500, 470),
              tip_length = 0, size = 1, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = wes_palette('Moonrise2')) +
  ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = 'grey95')) 
all_MOTU_ALPHA_Habitat_plot

### Save in pdf
ggsave(all_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/all_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 1-6. Plot MOTU richness vs Area
all_MOTU_ALPHA_Area_plot <- ggplot(richness_df,
                                   aes(x = Area,
                                       y = Richness,
                                       fill = Area)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
  ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        panel.background = element_rect(fill = 'grey95')) 
all_MOTU_ALPHA_Area_plot

### Save in pdf
ggsave(all_MOTU_ALPHA_Area_plot,
       file = '../../outputs/all_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')
# _______________________________________________



# ---- 2. Compute MOTU richness by kingdoms ----
## 2-1. Select the kingdoms
MOTUs_by_kingdom <- MOTU_table_final %>%
  filter(!is.na(kingdom)) %>%
  group_by(kingdom) %>%
  summarise(n_MOTUs = n(), .groups = 'drop')
MOTUs_by_kingdom

### Restrict to kingdoms with >= 10 MOTUs
kept_kingdom <- MOTUs_by_kingdom %>% 
  filter(n_MOTUs >= 10) %>% 
  select(kingdom) %>% 
  unlist()
kept_kingdom


## 2-2. Compute MOTU richness per kingdom
kingdom_MOTU_richness <- function(selected_kingdom){
  
  ### Restrict MOTU table to samples' columns
  pa_matrix <- MOTU_table_final %>%
    filter(kingdom == selected_kingdom) %>% 
    select(id, starts_with('SPY')) %>% 
    column_to_rownames('id')
    
    ### Transform into presence-absence_matrix
    pa_matrix[pa_matrix > 0] <- 1
    colSums(pa_matrix)
    final_matrix <- t(pa_matrix)
    
    ### Handle the problematic sample
    pooled_sample <- colSums(final_matrix[samples_to_pool, , drop = FALSE]) > 0
    pooled_sample <- as.numeric(pooled_sample)
    final_matrix_new <- final_matrix[!rownames(final_matrix) %in% samples_to_pool, ]
    final_matrix_new <- rbind(final_matrix_new,
                              setNames(as.data.frame(t(pooled_sample)),
                                       colnames(final_matrix)))
    rownames(final_matrix_new)[nrow(final_matrix_new)] <- pooled_name
    
    ### Compute MOTU richness
    richness_df <- data.frame(sample = rownames(final_matrix_new),
                              Richness = rowSums(final_matrix_new)) %>%
      inner_join(metadata_new) %>%
      select(Richness, pool, Habitat, Area) %>% 
      mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
             Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))
    
    return(richness_df)
}

### Metazoa
metazoa_richness_df <- kingdom_MOTU_richness('Metazoa')

### Chromista
chromista_richness_df <- kingdom_MOTU_richness('Chromista')

### Plantae
plantae_richness_df <- kingdom_MOTU_richness('Plantae')


## 2-3. Display descriptive statistics and variables distributions
data_explo <- function(richness_df){
  
  ### Descriptive statistics
  str(richness_df)
  summary(richness_df$Richness) # Y summary (Richness)
  summary(richness_df$Habitat) # X summary (Habitat)
  summary(richness_df$Area) # X summary (Area)
  
  ### Distribution of Y
  par(mfrow = c(3, 2))
  
  ### Boxplot
  boxplot(richness_df$Richness, col = 'blue', ylab = 'Richness')
  
  ### Cleveland plot
  dotchart(richness_df$Richness, pch = 16, col = 'blue', xlab = 'Richness')
  
  ### Histogram
  hist(richness_df$Richness, col = 'blue', xlab = 'Richness', main = '')
  
  ### Quantile-Quantile plot
  qqnorm(richness_df$Richness, pch = 16, col = 'blue', xlab = '')
  qqline(richness_df$Richness, col = 'red')
  
  ### Potential relationship Y/X
  #### Habitat
  boxplot(richness_df$Richness ~ richness_df$Habitat,
          varwidth = TRUE,
          ylab = 'MOTU richness',
          xlab = 'Habitat',
          col = wes_palette('Moonrise2')[1:3])
  
  #### Area
  boxplot(richness_df$Richness ~ richness_df$Area,
          varwidth = TRUE,
          ylab = 'MOTU richness',
          xlab = 'Area',
          col = wes_palette('AsteroidCity1')[c(1, 3, 4)])
}
### Metazoa
data_explo(metazoa_richness_df)
#### NOTE: Y follows a Gaussian distribution
#### NOTE: Richness looks lower in seaports
#### NOTE: No difference between areas is observed

### Chromista
data_explo(chromista_richness_df)
#### NOTE: Y follows a Gaussian distribution
#### NOTE: Richness looks lower in seaports
#### NOTE: No difference between areas is observed

### Plantae
data_explo(plantae_richness_df)
#### NOTE: Y follows a Gaussian distribution
#### NOTE: Richness looks lower in seaports
#### NOTE: No difference between areas is observed


## 2-4. Fit linear mixed models
model <- function(richness_df){
  
  ### Fit the model
  LMM_model <- lmer(Richness ~ Area + Habitat + (1|pool),
                    data = richness_df)
  
  ### Check normality of residuals
  par(mfrow = c(1, 2))
  
  #### Histogram
  hist(residuals(LMM_model), col = 'blue', xlab = 'Residuals', main = 'Check Normality')
  
  #### Quantile-Quantile plot
  qqnorm(residuals(LMM_model), pch = 16, col = 'blue', xlab = '')
  qqline(residuals(LMM_model), col = 'red')
  
  #### Shapiro test
  print(shapiro.test(residuals(LMM_model)))
  
  ### Check for homogeneity of variances (= Homoscedasticity)
  par(mfrow = c(1, 3))
  
  #### Residuals vs Fitted
  plot(residuals(LMM_model) ~ fitted(LMM_model), col = 'blue', pch = 16, ylab = 'Residuals', xlab = 'Fitted')
  abline(h = 0)
  
  #### Residuals vs Area
  boxplot(residuals(LMM_model) ~ richness_df$Area, varwidth = TRUE, ylab = 'Residuals', xlab = 'Area')
  abline(h = 0)
  
  #### Residuals vs Habitat
  boxplot(residuals(LMM_model) ~ richness_df$Habitat, varwidth = TRUE, ylab = 'Residuals', xlab = 'Habitat')
  abline(h = 0)
  
  ### Perform ANOVA
  print(Anova(LMM_model))
  
  ### Perform post-hoc Tukey's HSD
  emm <- emmeans(LMM_model, ~ Habitat)
  print(pairs(emm, adjust = 'tukey'))
  
  ### Obtain the R²
  r.squaredGLMM(LMM_model)
}

### Metazoa
model(metazoa_richness_df)
#### NOTE: Residuals follow normality
#### NOTE: There is no evident pattern in Residuals vs Fitted
#### NOTE: Effect of 'Area' is  not significant
#### NOTE: Effect of 'Habitat' is significant *
#### NOTE: Difference in richness is significant * between 'Seaport' and 'Sand'
####       Difference in richness is significant * between 'Seaport' and 'Rock'
####       Difference in richness is not significant between 'Rock' and 'Sand'
#### NOTE: Marginal R² = 77.05%
####       Conditional R² = 90.22%

### Chromista
model(chromista_richness_df)
#### NOTE: Residuals follow normality
#### NOTE: There is no evident pattern in Residuals vs Fitted
#### NOTE: Effect of 'Area' is  not significant
#### NOTE: Effect of 'Habitat' is not significant
#### NOTE: Marginal R² = 57.20%
####       Conditional R² = 90.66%

### Plantae
model(plantae_richness_df)
#### NOTE: Residuals follow normality
#### NOTE: There is no evident pattern in Residuals vs Fitted
#### NOTE: Effect of 'Area' is  not significant
#### NOTE: Effect of 'Habitat' is significant *
#### NOTE: Difference in richness is significant * between 'Seaport' and 'Sand'
####       Difference in richness is significant * between 'Seaport' and 'Rock'
####       Difference in richness is not significant between 'Rock' and 'Sand'
#### NOTE: Marginal R² = 77.05%
####       Conditional R² = 90.22%


## 2-5. Plot MOTU richness vs Habitat for each kingdom
### 2-5-1. Metazoa
metazoa_MOTU_ALPHA_Habitat_plot <- ggplot(metazoa_richness_df,
                                          aes(x = Habitat,
                                              y = Richness,
                                              fill = Habitat)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  geom_signif(comparisons = list(c('Sand', 'Seaport'),
                                 c('Rock', 'Seaport')),
              annotations = c('*', '*'),
              y_position = c(215, 205),
              tip_length = 0, size = 1, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = wes_palette('Moonrise2')) +
  ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank())
metazoa_MOTU_ALPHA_Habitat_plot

#### Save in pdf
ggsave(metazoa_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/metazoa_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### 2-5-2. Chromista
chromista_MOTU_ALPHA_Habitat_plot <- ggplot(chromista_richness_df,
                                            aes(x = Habitat,
                                                y = Richness,
                                                fill = Habitat)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  scale_fill_manual(values = wes_palette('Moonrise2')) +
  ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank())
chromista_MOTU_ALPHA_Habitat_plot

#### Save in pdf
ggsave(chromista_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/chromista_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### 2-5-3. Plantae
plantae_MOTU_ALPHA_Habitat_plot <- ggplot(plantae_richness_df,
                                          aes(x = Habitat,
                                              y = Richness,
                                              fill = Habitat)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  geom_signif(comparisons = list(c('Sand', 'Seaport'),
                                 c('Rock', 'Seaport')),
              annotations = c('*', '*'),
              y_position = c(33, 31),
              tip_length = 0, size = 1, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = wes_palette('Moonrise2')) +
  ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank())
plantae_MOTU_ALPHA_Habitat_plot

#### Save in pdf
ggsave(plantae_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/plantae_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 2-6. Plot MOTU richness vs Area for each kingdom
### 2-6-1. Metazoa
metazoa_MOTU_ALPHA_Area_plot <- ggplot(metazoa_richness_df,
                                       aes(x = Area,
                                           y = Richness,
                                           fill = Area)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
  ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank())
metazoa_MOTU_ALPHA_Area_plot

#### Save in pdf
ggsave(metazoa_MOTU_ALPHA_Area_plot,
       file = '../../outputs/metazoa_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### 2-6-2. Chromista
chromista_MOTU_ALPHA_Area_plot <- ggplot(chromista_richness_df,
                                         aes(x = Area,
                                             y = Richness,
                                             fill = Area)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
  ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank())
chromista_MOTU_ALPHA_Area_plot

#### Save in pdf
ggsave(chromista_MOTU_ALPHA_Area_plot,
       file = '../../outputs/chromista_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### 2-6-3. Plantae
plantae_MOTU_ALPHA_Area_plot <- ggplot(plantae_richness_df,
                                       aes(x = Area,
                                           y = Richness,
                                           fill = Area)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
  ylab('MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank())
plantae_MOTU_ALPHA_Area_plot

#### Save in pdf
ggsave(plantae_MOTU_ALPHA_Area_plot,
       file = '../../outputs/plantae_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')
# ______________________________________________



# ---- 3. Compute MOTU richness by phylum ----
## 3-1. Select the phyla
MOTUs_by_phylum <- MOTU_table_final %>%
  filter(!is.na(phylum)) %>%
  group_by(phylum) %>%
  summarise(n_MOTUs = n(), .groups = 'drop')
MOTUs_by_phylum

### Restrict to phylums with >= 10 MOTUs
kept_phylum <- MOTUs_by_phylum %>% 
  filter(n_MOTUs >= 10) %>% 
  select(phylum) %>% 
  unlist()
kept_phylum


## 3-2. Compute MOTU richness per phylum
phylum_MOTU_richness <- function(selected_phylum){
  
  ### Restrict MOTU table to samples' columns
  pa_matrix <- MOTU_table_final %>%
    filter(phylum == selected_phylum) %>% 
    select(id, starts_with('SPY')) %>% 
    column_to_rownames('id')
  
  ### Transform into presence-absence_matrix
  pa_matrix[pa_matrix > 0] <- 1
  colSums(pa_matrix)
  final_matrix <- t(pa_matrix)
  
  ### Handle the problematic sample
  pooled_sample <- colSums(final_matrix[samples_to_pool, , drop = FALSE]) > 0
  pooled_sample <- as.numeric(pooled_sample)
  final_matrix_new <- final_matrix[!rownames(final_matrix) %in% samples_to_pool, ]
  final_matrix_new <- rbind(final_matrix_new,
                            setNames(as.data.frame(t(pooled_sample)),
                                     colnames(final_matrix)))
  rownames(final_matrix_new)[nrow(final_matrix_new)] <- pooled_name
  
  ### Compute MOTU richness
  richness_df <- data.frame(sample = rownames(final_matrix_new),
                            Richness = rowSums(final_matrix_new)) %>%
    inner_join(metadata_new) %>%
    select(Richness, pool, Habitat, Area) %>% 
    mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
           Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))
  
  return(richness_df)
}

### Annelida
annelida_richness_df <- phylum_MOTU_richness('Annelida')

### Arthropoda
arthropoda_richness_df <- phylum_MOTU_richness('Arthropoda')

### Chordata
chordata_richness_df <- phylum_MOTU_richness('Chordata')

### Cnidaria
cnidaria_richness_df <- phylum_MOTU_richness('Cnidaria')

### Mollusca
mollusca_richness_df <- phylum_MOTU_richness('Mollusca')

### Nemertea
nemertea_richness_df <- phylum_MOTU_richness('Nemertea')

### Bacillariophyta
bacillariophyta_richness_df <- phylum_MOTU_richness('Bacillariophyta')

### Ochrophyta
ochrophyta_richness_df <- phylum_MOTU_richness('Ochrophyta')

### Oomycota
oomycota_richness_df <- phylum_MOTU_richness('Oomycota')

### Rhodophyta
rhodophyta_richness_df <- phylum_MOTU_richness('Rhodophyta')


## 3-3. Plot MOTU richness vs Area for each phylum
richness_vs_Habitat <- function(phylum_richness_df){
  phylum_MOTU_ALPHA_Habitat_plot <- ggplot(phylum_richness_df,
                                        aes(x = Habitat,
                                            y = Richness,
                                            fill = Habitat)) +
    geom_boxplot(width = 0.5, show.legend = F) +
    scale_fill_manual(values = wes_palette('Moonrise2')) +
    ylab('MOTU richness') +
    theme_classic() +
    theme(axis.text.y = element_text(colour = 'black', size = 10),
          axis.text.x = element_text(colour = 'black', size = 10),
          axis.title.y = element_text(colour = 'black', size = 15),
          axis.title.x = element_blank())
  
  return(phylum_MOTU_ALPHA_Habitat_plot)
}

### Annelida
annelida_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(annelida_richness_df)

#### Save in pdf
ggsave(annelida_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/annelida_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Arthropoda
arthropoda_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(arthropoda_richness_df)

#### Save in pdf
ggsave(arthropoda_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/arthropoda_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Chordata
chordata_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(chordata_richness_df)

#### Save in pdf
ggsave(chordata_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/chordata_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Cnidaria
cnidaria_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(cnidaria_richness_df)

#### Save in pdf
ggsave(cnidaria_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/cnidaria_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Mollusca
mollusca_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(mollusca_richness_df)

#### Save in pdf
ggsave(mollusca_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/mollusca_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Nemertea
nemertea_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(nemertea_richness_df)

#### Save in pdf
ggsave(nemertea_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/nemertea_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Bacillariophyta
bacillariophyta_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(bacillariophyta_richness_df)

#### Save in pdf
ggsave(bacillariophyta_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/bacillariophyta_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Ochrophyta
ochrophyta_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(ochrophyta_richness_df)

#### Save in pdf
ggsave(ochrophyta_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/ochrophyta_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Oomycota
oomycota_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(oomycota_richness_df)

#### Save in pdf
ggsave(oomycota_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/oomycota_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Rhodophyta
rhodophyta_MOTU_ALPHA_Habitat_plot <- richness_vs_Habitat(rhodophyta_richness_df)

#### Save in pdf
ggsave(rhodophyta_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/rhodophyta_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 3-4. Plot MOTU richness vs Area for each phylum
richness_vs_Area <- function(phylum_richness_df){
  phylum_MOTU_ALPHA_Area_plot <- ggplot(phylum_richness_df,
                                        aes(x = Area,
                                            y = Richness,
                                            fill = Area)) +
    geom_boxplot(width = 0.5, show.legend = F) +
    scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
    ylab('MOTU richness') +
    theme_classic() +
    theme(axis.text.y = element_text(colour = 'black', size = 10),
          axis.text.x = element_text(colour = 'black', size = 10),
          axis.title.y = element_text(colour = 'black', size = 15),
          axis.title.x = element_blank())
  
  return(phylum_MOTU_ALPHA_Area_plot)
}

### Annelida
annelida_MOTU_ALPHA_Area_plot <- richness_vs_Area(annelida_richness_df)

#### Save in pdf
ggsave(annelida_MOTU_ALPHA_Area_plot,
       file = '../../outputs/annelida_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Arthropoda
arthropoda_MOTU_ALPHA_Area_plot <- richness_vs_Area(arthropoda_richness_df)

#### Save in pdf
ggsave(arthropoda_MOTU_ALPHA_Area_plot,
       file = '../../outputs/arthropoda_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Chordata
chordata_MOTU_ALPHA_Area_plot <- richness_vs_Area(chordata_richness_df)

#### Save in pdf
ggsave(chordata_MOTU_ALPHA_Area_plot,
       file = '../../outputs/chordata_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Cnidaria
cnidaria_MOTU_ALPHA_Area_plot <- richness_vs_Area(cnidaria_richness_df)

#### Save in pdf
ggsave(cnidaria_MOTU_ALPHA_Area_plot,
       file = '../../outputs/cnidaria_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Mollusca
mollusca_MOTU_ALPHA_Area_plot <- richness_vs_Area(mollusca_richness_df)

#### Save in pdf
ggsave(mollusca_MOTU_ALPHA_Area_plot,
       file = '../../outputs/mollusca_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Nemertea
nemertea_MOTU_ALPHA_Area_plot <- richness_vs_Area(nemertea_richness_df)

#### Save in pdf
ggsave(nemertea_MOTU_ALPHA_Area_plot,
       file = '../../outputs/nemertea_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Bacillariophyta
bacillariophyta_MOTU_ALPHA_Area_plot <- richness_vs_Area(bacillariophyta_richness_df)

#### Save in pdf
ggsave(bacillariophyta_MOTU_ALPHA_Area_plot,
       file = '../../outputs/bacillariophyta_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Ochrophyta
ochrophyta_MOTU_ALPHA_Area_plot <- richness_vs_Area(ochrophyta_richness_df)

#### Save in pdf
ggsave(ochrophyta_MOTU_ALPHA_Area_plot,
       file = '../../outputs/ochrophyta_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Oomycota
oomycota_MOTU_ALPHA_Area_plot <- richness_vs_Area(oomycota_richness_df)

#### Save in pdf
ggsave(oomycota_MOTU_ALPHA_Area_plot,
       file = '../../outputs/oomycota_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')

### Rhodophyta
rhodophyta_MOTU_ALPHA_Area_plot <- richness_vs_Area(rhodophyta_richness_df)

#### Save in pdf
ggsave(rhodophyta_MOTU_ALPHA_Area_plot,
       file = '../../outputs/rhodophyta_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')
# ____________________________________________