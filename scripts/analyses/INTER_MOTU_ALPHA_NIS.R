# ---- 0. Load packages ----
library(tidyverse)
library(wesanderson)
library(lme4) # linear mixed-models
library(car) # compute ANOVA
library(emmeans) # Tukey's HSD test
library(MuMIn) # obtain R²
library(ggsignif)
# __________________________



# ---- 1. Compute MOTU richness NIS MOTUs ----
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


## 1-2. Compute relative NIS MOTU richness
### Calculate NIS richness
NIS_ids <- which(MOTU_table_final$NIS == 'NIS')
NIS_richness <- rowSums(final_matrix_new[, NIS_ids])

### Calculate total richness
tot_richness <- rowSums(final_matrix_new)

### Relative richness table
richness_df <- data.frame(sample = rownames(final_matrix_new),
                          Richness = NIS_richness / tot_richness) %>%
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

### NOTE: Y follows a Gaussian distribution

dev.off()

### 1-3-3. Potential relationship Y/X
#### Habitat
boxplot(richness_df$Richness ~ richness_df$Habitat,
        varwidth = TRUE,
        ylab = 'MOTU richness',
        xlab = 'Habitat',
        col = wes_palette('Moonrise2')[1:3])
#### NOTE: Richness looks lower in sand

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
#### NOTE: Effect of 'Habitat' is significant ***

### 1-4-4. Perform post-hoc Tukey's HSD
emm <- emmeans(LMM_model, ~ Habitat)
pairs(emm, adjust = 'tukey')
#### NOTE: Difference in richness is significant * between 'Seaport' and 'Sand'
####       Difference in richness is not significant between 'Seaport' and 'Rock'
####       Difference in richness is not significant between 'Rock' and 'Sand'

### 1-4-5. Obtain the R²
r.squaredGLMM(LMM_model)
#### NOTE: Marginal R² = 66.04%
####       Conditional R² = 70.49%


## 1-5. Plot MOTU richness vs Habitat
NIS_MOTU_ALPHA_Habitat_plot <- ggplot(richness_df,
                                      aes(x = Habitat,
                                          y = Richness,
                                          fill = Habitat)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  geom_signif(comparisons = list(c('Sand', 'Seaport')),
              annotations = c('*'),
              y_position = c(0.055),
              tip_length = 0, size = 1, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = wes_palette('Moonrise2')) +
  ylab('Relative NIS MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank())
NIS_MOTU_ALPHA_Habitat_plot

### Save in pdf
ggsave(NIS_MOTU_ALPHA_Habitat_plot,
       file = '../../outputs/NIS_MOTU_ALPHA_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 1-6. Plot MOTU richness vs Area
NIS_MOTU_ALPHA_Area_plot <- ggplot(richness_df,
                                   aes(x = Area,
                                       y = Richness,
                                       fill = Area)) +
  geom_boxplot(width = 0.5, show.legend = F) +
  scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
  ylab('Relative NIS MOTU richness') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank())
NIS_MOTU_ALPHA_Area_plot

### Save in pdf
ggsave(NIS_MOTU_ALPHA_Area_plot,
       file = '../../outputs/NIS_MOTU_ALPHA_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')
# _______________________________________________