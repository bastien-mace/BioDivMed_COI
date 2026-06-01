# ---- 0. Load packages ----
library(tidyverse)
library(Biostrings) # handle DNA sequences
library(DECIPHER) # align DNA sequences
library(ape) # measure genetic distances
library(wesanderson)
library(lme4) # linear mixed-models
library(glmmTMB) # generalized linear mixed-models
library(performance) # check overdispersion
library(car) # compute ANOVA
library(emmeans) # Tukey's HSD test
library(MuMIn) # obtain R²
library(ggsignif)
# __________________________



# ---- 1. Import data ----
## 1-1. ASV table
ASV_table_final <- read.csv('../../data/ASV_table_final.csv')


## 1-2. Presence/Absence matrix
### Restrict MOTU table to samples' columns
MOTU_table_final <- read.csv('../../data/MOTU_table_final.csv')
pa_matrix <- MOTU_table_final %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

### Transform into presence-absence_matrix
pa_matrix[pa_matrix > 0] <- 1
colSums(pa_matrix)
final_matrix <- t(pa_matrix) # transpose matrix to get samples in rows and MOTUs in columns


## 1-3. Site metadata
metadata <- read.csv2('../../data/site_metadata.csv') %>%
  pivot_longer(cols = c(spygen_code_1, spygen_code_2),
               values_to = 'sample') %>% 
  mutate(sample = factor(sample, levels = rownames(final_matrix))) %>%
  arrange(sample)


## 1-4. Pool the problematic sample
samples_to_pool <- c('SPY232133', 'SPY232146')
pooled_sample <- colSums(final_matrix[samples_to_pool, , drop = FALSE]) > 0
pooled_sample <- as.numeric(pooled_sample)
pooled_name <- 'SPY232133_SPY232146'

### Remove the samples pooled and add the corresponding pool
final_matrix_new <- final_matrix[!rownames(final_matrix) %in% samples_to_pool, ]
final_matrix_new <- rbind(final_matrix_new, pooled_sample)
rownames(final_matrix_new)[nrow(final_matrix_new)] <- pooled_name

### Modify metadata
metadata_new <- metadata[!metadata$sample %in% samples_to_pool, ]
pooled_meta <- metadata[metadata$sample == samples_to_pool[1], ]
pooled_meta$sample <- pooled_name
metadata_new <- rbind(metadata_new, pooled_meta)
# ________________________



# ---- 2. Compute genetic diversity metrics ----
## 2-1. Build an empty list of matrix
### Calculate the number of ASVs per MOTU (for Hd and pi)
nb_ASV_per_MOTU <- ASV_table_final %>%
  filter(MOTU %in% MOTU_table_final$id) %>% 
  count(MOTU, name = 'n_ASV') %>% 
  filter(n_ASV > 1) # only keep polymorphic MOTUs

intra_alpha_metrics <- list(ASV_richness = matrix(NA,
                                                  nrow = dim(final_matrix_new)[2],
                                                  ncol = dim(final_matrix_new)[1],
                                                  dimnames = list(colnames(final_matrix_new),
                                                                  rownames(final_matrix_new))),
                            Hd = matrix(NA,
                                        nrow = dim(nb_ASV_per_MOTU)[1],
                                        ncol = dim(final_matrix_new)[1],
                                        dimnames = list(nb_ASV_per_MOTU$MOTU,
                                                        rownames(final_matrix_new))),
                            pi = matrix(NA,
                                        nrow = dim(nb_ASV_per_MOTU)[1],
                                        ncol = dim(final_matrix_new)[1],
                                        dimnames = list(nb_ASV_per_MOTU$MOTU,
                                                        rownames(final_matrix_new))))


## 2-2. Complete the ASV richness matrix MOTU by MOTU
for (i in 1:dim(final_matrix_new)[2]){
  
  #### Select MOTU
  MOTU_id <- colnames(final_matrix)[i]
  
  #### Extract its ASVs presence/absence table
  ASVs_table <- ASV_table_final %>% 
    filter(MOTU == MOTU_id) %>% 
    select(id, starts_with('SPY')) %>% 
    column_to_rownames('id')
  ASVs_table[ASVs_table > 0] <- 1
  ASVs_table <- t(ASVs_table)
  
  #### Pool the problematic sample
  pooled_sample <- colSums(ASVs_table[samples_to_pool, , drop = FALSE]) > 0
  pooled_sample <- as.numeric(pooled_sample)
  final_ASVs_table <- ASVs_table[!rownames(ASVs_table) %in% samples_to_pool, ]
  
  if(is.null(dim(final_ASVs_table))){ # there is only one ASV in the MOTU
    final_ASVs_table <- c(final_ASVs_table, pooled_sample)
    names(final_ASVs_table)[length(final_ASVs_table)] <- pooled_name
    richness <- final_ASVs_table
    
  }else{ # there are 2 ASVs or more in the MOTU
    final_ASVs_table <- rbind(final_ASVs_table, pooled_sample)
    rownames(final_ASVs_table)[nrow(final_ASVs_table)] <- pooled_name
    richness <- rowSums(final_ASVs_table)
  }
  
  ### Only consider MOTU that are present 
  motu_presence <- final_matrix_new[, i] # MOTU presence/absence across samples
  richness[motu_presence == 0] <- NA  # Replace value by NA if the MOTU is absent
  
  ### Add ASV richness
  intra_alpha_metrics$ASV_richness[i, ] <- richness
}


## 2-3. Compute haplotype diversity (Hd) and nucleotide diversity (π) matrices MOTU by MOTU
for (i in 1:dim(nb_ASV_per_MOTU)[1]){
  
  ### Select MOTU
  MOTU_id <- nb_ASV_per_MOTU$MOTU[i]
  
  ### Extract its ASVs presence/absence table
  ASVs_table <- ASV_table_final %>% 
    filter(MOTU == MOTU_id) %>% 
    select(id, starts_with('SPY')) %>% 
    column_to_rownames('id')
  ASVs_table <- t(ASVs_table)
  
  ### Extract ASVs' sequences
  ASVs_seq <- ASV_table_final %>% 
    filter(MOTU == MOTU_id) %>%
    select(id, seq)
  
  ### Pool the problematic sample
  pooled_sample <- colSums(ASVs_table[samples_to_pool, , drop = FALSE]) # somme des abondances
  final_ASVs_table <- ASVs_table[!rownames(ASVs_table) %in% samples_to_pool, ]
  final_ASVs_table <- rbind(final_ASVs_table, pooled_sample)
  rownames(final_ASVs_table)[nrow(final_ASVs_table)] <- pooled_name
    
  ### Calculate haplotype diversity (Hd) per sample
  Hd <- apply(final_ASVs_table, 1, function(x){
      
    #### Only consider ASVs present in the sample
    x <- x[x > 0]
    n <- sum(x)
    if(length(x) <= 1) return(0) # Hd is 0 if all ASVs are absent in the sample or if there is only one ASV
    p <- x / n # ASVs' frequencies
    Hd_val <- (n / (n - 1)) * (1 - sum(p^2))
    return(Hd_val)
    })
    
  ### Calculate nucleotide diversity (π) per sample
  #### Align ASV sequences
  dna_seqs <- DNAStringSet(ASVs_seq$seq)
  names(dna_seqs) <- ASVs_seq$id
  dna_seqs_aligned <- AlignSeqs(dna_seqs)
  
  #### Compute nucleotide distance between ASVs
  dist_mat <- as.matrix(dist.dna(as.DNAbin(dna_seqs_aligned), model = 'raw'))
  
  pi <- apply(final_ASVs_table, 1, function(x){
    
    #### Only consider ASVs present in the sample
    x <- x[x > 0]
    n <- sum(x)
    if(length(x) <= 1) return(0) # π is 0 if all ASVs are absent in the sample or if there is only one ASV
    p <- x / n # ASVs' frequencies
    ids <- names(x)
    d_sub <- dist_mat[ids, ids] # distance matrix of ASVs present in the sample
    pi_val <- 0
    for(a in 1:length(p)){
      for(b in 1:length(p)){
        if(a != b){
          pi_val <- pi_val + p[a] * p[b] * d_sub[a, b]
        }
      }
    }
    return(pi_val)
    })
  
  ### Only consider MOTU that are present
  motu_col_index <- which(colnames(final_matrix_new) == MOTU_id)
  motu_presence <- final_matrix_new[, motu_col_index] # MOTU presence/absence across samples
  Hd[motu_presence == 0] <- NA   # Replace value by NA if the MOTU is absent
  pi[motu_presence == 0] <- NA
  
  ### Add Hd and π
  intra_alpha_metrics$Hd[i, ] <- Hd
  intra_alpha_metrics$pi[i, ] <- pi
}
# ______________________________________________
 


# ---- 3. ASV richness analysis ----  
## 3-1. Prepare a long dataframe
ASV_richness_long <- intra_alpha_metrics$ASV_richness %>%
  as.data.frame() %>% 
  rownames_to_column('MOTU') %>%
  pivot_longer(cols = -MOTU,
               names_to = 'sample',
               values_to = 'Richness') %>%
  inner_join(metadata_new) %>%
  select(Richness, pool, Habitat, Area, MOTU) %>% 
  drop_na() %>% 
  mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
         Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')),
         pool = factor(pool),
         MOTU = factor(MOTU))


## 3-2. Data exploration
### 3-2-1. Descriptive statistics
str(ASV_richness_long)
summary(ASV_richness_long$Richness) # Y summary (Richness)
summary(ASV_richness_long$Habitat) # X summary (Habitat)
summary(ASV_richness_long$Area) # X summary (Area)

### 3-2-2. Distribution of Y
par(mfrow = c(2,2))

#### Boxplot
boxplot(ASV_richness_long$Richness, col = 'blue', ylab = 'Richness')

#### Cleveland plot
dotchart(ASV_richness_long$Richness, pch = 16, col = 'blue', xlab = 'Richness')
##### NOTE: A log scale looks more appropriate

#### Histogram
hist(log(ASV_richness_long$Richness), col = 'blue', xlab = 'Richness', main = '')

#### Quantile-Quantile plot
qqnorm(log(ASV_richness_long$Richness), pch = 16, col = 'blue', xlab = '')
qqline(log(ASV_richness_long$Richness), col = 'red')
##### NOTE: Y does not follow a Gaussian distribution (right-skewed)

dev.off()

### 3-2-3. Potential relationship Y/X
#### Habitat
boxplot(log(ASV_richness_long$Richness) ~ ASV_richness_long$Habitat,
        varwidth = TRUE,
        ylab = 'ASV richness per MOTU',
        xlab = 'Habitat',
        col = wes_palette('Moonrise2')[1:3])
#### NOTE: Richness looks lower in seaports

#### Area
boxplot(log(ASV_richness_long$Richness) ~ ASV_richness_long$Area,
        varwidth = TRUE,
        ylab = 'ASV richness per MOTU',
        xlab = 'Area',
        col = wes_palette('AsteroidCity1')[c(1, 3, 4)])
#### NOTE: No difference between areas is observed


### 3-3. Fit a linear mixed-model (LMM)
LMM_model <- lmer(log(Richness) ~ Area + Habitat + (1|pool) + (1|MOTU),
                  data = ASV_richness_long)

### 3-3-1. Check normality of residuals
par(mfrow = c(1, 2))

#### Histogram
hist(residuals(LMM_model), col = 'blue', xlab = 'Residuals', main = 'Check Normality')

#### Quantile-Quantile plot
qqnorm(residuals(LMM_model), pch = 16, col = 'blue', xlab = '')
qqline(residuals(LMM_model), col = 'red')

#### Shapiro test
shapiro.test(residuals(LMM_model))
#### NOTE: Residuals does not follow normality

### 3-3-2. Check for homogeneity of variances (= Homoscedasticity)
par(mfrow = c(1, 3))

#### Residuals vs Fitted
plot(residuals(LMM_model) ~ fitted(LMM_model), col = 'blue', pch = 16, ylab = 'Residuals', xlab = 'Fitted')
abline(h = 0)

#### Residuals vs Area
boxplot(residuals(LMM_model) ~ ASV_richness_long$Area, varwidth = TRUE, ylab = 'Residuals', xlab = 'Area')
abline(h = 0)

#### Residuals vs Habitat
boxplot(residuals(LMM_model) ~ ASV_richness_long$Habitat, varwidth = TRUE, ylab = 'Residuals', xlab = 'Habitat')
abline(h = 0)
##### NOTE: There is an evident pattern in Residuals vs Fitted


### 3-4. Fit a generalized linear mixed-model (GLMM)
#### 3-4-1. Poisson GLMM
GLMM_model <- glmmTMB(Richness ~ Area + Habitat + (1|pool) + (1|MOTU),
                      data = ASV_richness_long,
                      family = poisson)

##### Check overdispersion
check_overdispersion(GLMM_model)
##### NOTE: There is a slight overdispersion (1.315)

#### 3-4-2. Negative binomial GLMM
GLMM_model <- glmmTMB(Richness ~ Area + Habitat + (1|pool) + (1|MOTU),
                      data = ASV_richness_long, family = nbinom2)
summary(GLMM_model)

### 3-4-3. Perform ANOVA
Anova(GLMM_model)
#### NOTE: Effect of 'Area' is  not significant
#### NOTE: Effect of 'Habitat' is significant ***

### 3-4-4. Perform post-hoc Tukey's HSD
emm <- emmeans(GLMM_model, ~ Habitat)
pairs(emm, adjust = 'tukey')
#### NOTE: Difference in ASV richness per MOTU is significant *** between 'Seaport' and 'Sand'
####       Difference in ASV richness per MOTU is significant *** between 'Seaport' and 'Rock'
####       Difference in ASV richness per MOTU is not significant between 'Rock' and 'Sand'

### 3-4-5. Obtain the R²
r.squaredGLMM(GLMM_model)
#### NOTE: Marginal R² = 0.77%
####       Conditional R² = 54.47%


## 3-5. Plot ASV richness per MOTU vs Habitat
### 3-5-1. Extract marginal effects
emm_habitat_df <- as.data.frame(emmeans(GLMM_model, ~ Habitat, type = 'response'))

### 3-5-2. Plot the marginal effects
INTRA_MOTU_ASV_richness_Habitat_plot <- ggplot(emm_habitat_df,
                                               aes(x = Habitat,
                                                   y = response,
                                                   fill = Habitat)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, linewidth = 0.75, color = 'black') +
  geom_point(size = 5, shape = 21, show.legend = FALSE) +
  geom_signif(comparisons = list(c('Sand', 'Seaport'),
                                 c('Rock', 'Seaport')),
              annotations = c('***', '***'),
              y_position = c(1.9, 1.86),
              tip_length = 0, size = 1, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = wes_palette('Moonrise2')) +
  ylab('Model-adjusted value') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = 'grey95'))
INTRA_MOTU_ASV_richness_Habitat_plot

### Save in pdf
ggsave(INTRA_MOTU_ASV_richness_Habitat_plot,
       file = '../../outputs/INTRA_MOTU_ASV_richness_Habitat_plot.pdf',
       height = 8,
       width = 6,
       unit = 'cm')


## 3-5. Plot ASV richness per MOTU vs Area
### 3-5-1. Extract marginal effects
emm_area_df <- as.data.frame(emmeans(GLMM_model, ~ Area, type = 'response'))

### 3-5-2. Plot the marginal effects
INTRA_MOTU_ASV_richness_Area_plot <- ggplot(emm_area_df,
                                            aes(x = Area,
                                                y = response,
                                                fill = Area)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, linewidth = 0.75, color = 'black') +
  geom_point(size = 5, shape = 21, show.legend = FALSE) +
  scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
  ylab('Model-adjusted value') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = 'grey95'))
INTRA_MOTU_ASV_richness_Area_plot

### Save in pdf
ggsave(INTRA_MOTU_ASV_richness_Area_plot,
       file = '../../outputs/INTRA_MOTU_ASV_richness_Area_plot.pdf',
       height = 8,
       width = 6,
       unit = 'cm')
# _______________________________________________



# ---- 4. Haplotype diversity (Hd) analysis ----  
## 4-1. Prepare a long dataframe
Hd_long <- intra_alpha_metrics$Hd %>%
  as.data.frame() %>% 
  rownames_to_column('MOTU') %>%
  pivot_longer(cols = -MOTU,
               names_to = 'sample',
               values_to = 'Hd') %>%
  inner_join(metadata_new) %>%
  select(Hd, pool, Habitat, Area, MOTU) %>% 
  drop_na() %>% 
  mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
         Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')),
         pool = factor(pool),
         MOTU = factor(MOTU))


## 4-2. Data exploration
### 4-2-1. Descriptive statistics
str(Hd_long)
summary(Hd_long$Hd) # Y summary (Hd)
summary(Hd_long$Habitat) # X summary (Habitat)
summary(Hd_long$Area) # X summary (Area)

### 4-2-2. Distribution of Y
par(mfrow = c(2,2))

#### Boxplot
boxplot(Hd_long$Hd, col = 'blue', ylab = 'Hd')

#### Cleveland plot
dotchart(Hd_long$Hd, pch = 16, col = 'blue', xlab = 'Hd')

#### Histogram
hist(Hd_long$Hd, col = 'blue', xlab = 'Hd', main = '')

#### Quantile-Quantile plot
qqnorm(Hd_long$Hd, pch = 16, col = 'blue', xlab = '')
qqline(Hd_long$Hd, col = 'red')
##### NOTE: Y does not follow a Gaussian distribution (right-skewed)

dev.off()

### 4-2-3. Potential relationship Y/X
#### Habitat
boxplot(Hd_long$Hd ~ Hd_long$Habitat,
        varwidth = TRUE,
        ylab = 'Hd per MOTU',
        xlab = 'Habitat',
        col = wes_palette('Moonrise2')[1:3])
#### NOTE: No difference between habitats is observed

#### Area
boxplot(Hd_long$Hd ~ Hd_long$Area,
        varwidth = TRUE,
        ylab = 'Hd per MOTU',
        xlab = 'Area',
        col = wes_palette('AsteroidCity1')[c(1, 3, 4)])
#### NOTE: No difference between areas is observed


### 4-3. Fit a linear mixed-model (LMM)
LMM_model <- lmer(Hd ~ Area + Habitat + (1|pool) + (1|MOTU),
                  data = Hd_long)

### 4-3-1. Check normality of residuals
par(mfrow = c(1, 2))

#### Histogram
hist(residuals(LMM_model), col = 'blue', xlab = 'Residuals', main = 'Check Normality')

#### Quantile-Quantile plot
qqnorm(residuals(LMM_model), pch = 16, col = 'blue', xlab = '')
qqline(residuals(LMM_model), col = 'red')

#### Shapiro test
shapiro.test(residuals(LMM_model))
#### NOTE: Residuals follow normality

### 4-3-2. Check for homogeneity of variances (= Homoscedasticity)
par(mfrow = c(1, 3))

#### Residuals vs Fitted
plot(residuals(LMM_model) ~ fitted(LMM_model), col = 'blue', pch = 16, ylab = 'Residuals', xlab = 'Fitted')
abline(h = 0)

#### Residuals vs Area
boxplot(residuals(LMM_model) ~ Hd_long$Area, varwidth = TRUE, ylab = 'Residuals', xlab = 'Area')
abline(h = 0)

#### Residuals vs Habitat
boxplot(residuals(LMM_model) ~ Hd_long$Habitat, varwidth = TRUE, ylab = 'Residuals', xlab = 'Habitat')
abline(h = 0)
##### NOTE: There is an evident pattern in Residuals vs Fitted


### 4-4. Fit a beta generalized linear mixed-model (GLMM)
#### Apply Smithson & Verkuilen transformation (10.1037/1082-989X.11.1.54)
n <- nrow(Hd_long)
Hd_long$Hd_beta <- (Hd_long$Hd * (n - 1) + 0.5) / n
GLMM_model <- glmmTMB(Hd_beta ~ Area + Habitat + (1|pool) + (1|MOTU),
                      data = Hd_long,
                      family = beta_family())
summary(GLMM_model)

### 4-4-1. Perform ANOVA
Anova(GLMM_model)
#### NOTE: Effect of 'Area' is  not significant
#### NOTE: Effect of 'Habitat' is significant ***

### 4-4-2. Perform post-hoc Tukey's HSD
emm <- emmeans(GLMM_model, ~ Habitat)
pairs(emm, adjust = 'tukey')
#### NOTE: Difference in Hd per MOTU is significant *** between 'Seaport' and 'Sand'
####       Difference in Hd per MOTU is significant *** between 'Seaport' and 'Rock'
####       Difference in Hd per MOTU is not significant between 'Rock' and 'Sand'

### 4-4-4. Obtain the R²
r.squaredGLMM(GLMM_model)
#### NOTE: Marginal R² = 3.18%
####       Conditional R² = 35.54%


## 4-5. Plot ASV richness per MOTU vs Habitat
### 4-5-1. Extract marginal effects
emm_habitat_df <- as.data.frame(emmeans(GLMM_model, ~ Habitat, type = 'response'))

### 4-5-2. Plot the marginal effects
INTRA_MOTU_Hd_Habitat_plot <- ggplot(emm_habitat_df,
                                     aes(x = Habitat,
                                         y = response,
                                         fill = Habitat)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, linewidth = 0.75, color = 'black') +
  geom_point(size = 5, shape = 21, show.legend = FALSE) +
  geom_signif(comparisons = list(c('Sand', 'Seaport'),
                                 c('Rock', 'Seaport')),
              annotations = c('***', '***'),
              y_position = c(0.31, 0.3),
              tip_length = 0, size = 1, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = wes_palette('Moonrise2')) +
  ylab('Model-adjusted value') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = 'grey95'))
INTRA_MOTU_Hd_Habitat_plot

### Save in pdf
ggsave(INTRA_MOTU_Hd_Habitat_plot,
       file = '../../outputs/INTRA_MOTU_Hd_Habitat_plot.pdf',
       height = 8,
       width = 6,
       unit = 'cm')


## 4-6. Plot ASV richness per MOTU vs Area
### 4-6-1. Extract marginal effects
emm_area_df <- as.data.frame(emmeans(GLMM_model, ~ Area, type = 'response'))

### 4-6-2. Plot the marginal effects
INTRA_MOTU_Hd_Area_plot <- ggplot(emm_area_df,
                                     aes(x = Area,
                                         y = response,
                                         fill = Area)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, linewidth = 0.75, color = 'black') +
  geom_point(size = 5, shape = 21, show.legend = FALSE) +
  scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
  ylab('Model-adjusted value') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = 'grey95'))
INTRA_MOTU_Hd_Area_plot

### Save in pdf
ggsave(INTRA_MOTU_Hd_Area_plot,
       file = '../../outputs/INTRA_MOTU_Hd_Area_plot.pdf',
       height = 8,
       width = 6,
       unit = 'cm')
# ______________________________________________



# ---- 5. Nucleotide diversity (π) analysis ----  
## 5-1. Prepare a long dataframe
pi_long <- intra_alpha_metrics$pi %>%
  as.data.frame() %>% 
  rownames_to_column('MOTU') %>%
  pivot_longer(cols = -MOTU,
               names_to = 'sample',
               values_to = 'pi') %>%
  inner_join(metadata_new) %>%
  select(pi, pool, Habitat, Area, MOTU) %>% 
  drop_na() %>% 
  mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
         Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')),
         pool = factor(pool),
         MOTU = factor(MOTU))


## 5-2. Data exploration
### 5-2-1. Descriptive statistics
str(pi_long)
summary(pi_long$pi) # Y summary (pi)
summary(pi_long$Habitat) # X summary (Habitat)
summary(pi_long$Area) # X summary (Area)

### 5-2-2. Distribution of Y
par(mfrow = c(2,2))

#### Boxplot
boxplot(pi_long$pi, col = 'blue', ylab = 'pi')

#### Cleveland plot
dotchart(pi_long$pi, pch = 16, col = 'blue', xlab = 'pi')

#### Histogram
hist(pi_long$pi, col = 'blue', xlab = 'pi', main = '')

#### Quantile-Quantile plot
qqnorm(pi_long$pi, pch = 16, col = 'blue', xlab = '')
qqline(pi_long$pi, col = 'red')
##### NOTE: Y does not follow a Gaussian distribution (right-skewed)

dev.off()

### 5-2-3. Potential relationship Y/X
#### Habitat
boxplot(pi_long$pi ~ pi_long$Habitat,
        varwidth = TRUE,
        ylab = 'π per MOTU',
        xlab = 'Habitat',
        col = wes_palette('Moonrise2')[1:3])
#### NOTE: No difference between habitats is observed

#### Area
boxplot(pi_long$pi ~ pi_long$Area,
        varwidth = TRUE,
        ylab = 'π per MOTU',
        xlab = 'Area',
        col = wes_palette('AsteroidCity1')[c(1, 3, 4)])
#### NOTE: No difference between areas is observed


### 5-3. Fit a linear mixed-model (LMM)
LMM_model <- lmer(pi ~ Area + Habitat + (1|pool) + (1|MOTU),
                  data = pi_long)

### 5-3-1. Check normality of residuals
par(mfrow = c(1, 2))

#### Histogram
hist(residuals(LMM_model), col = 'blue', xlab = 'Residuals', main = 'Check Normality')

#### Quantile-Quantile plot
qqnorm(residuals(LMM_model), pch = 16, col = 'blue', xlab = '')
qqline(residuals(LMM_model), col = 'red')

#### Shapiro test
shapiro.test(residuals(LMM_model))
#### NOTE: Residuals follow normality

### 5-3-2. Check for homogeneity of variances (= Homoscedasticity)
par(mfrow = c(1, 3))

#### Residuals vs Fitted
plot(residuals(LMM_model) ~ fitted(LMM_model), col = 'blue', pch = 16, ylab = 'Residuals', xlab = 'Fitted')
abline(h = 0)

#### Residuals vs Area
boxplot(residuals(LMM_model) ~ pi_long$Area, varwidth = TRUE, ylab = 'Residuals', xlab = 'Area')
abline(h = 0)

#### Residuals vs Habitat
boxplot(residuals(LMM_model) ~ pi_long$Habitat, varwidth = TRUE, ylab = 'Residuals', xlab = 'Habitat')
abline(h = 0)
##### NOTE: There is an evident pattern in Residuals vs Fitted


### 5-4. Fit a beta generalized linear mixed-model (GLMM)
#### Apply Smithson & Verkuilen transformation (10.1037/1082-989X.11.1.54)
n <- nrow(pi_long)
pi_long$pi_beta <- (pi_long$pi * (n - 1) + 0.5) / n
GLMM_model <- glmmTMB(pi_beta ~ Area + Habitat + (1|pool) + (1|MOTU),
                      data = pi_long,
                      family = beta_family())
summary(GLMM_model)

### 5-4-1. Perform ANOVA
Anova(GLMM_model)
#### NOTE: Effect of 'Area' is significant *
#### NOTE: Effect of 'Habitat' is significant ***

### 4-4-2. Perform post-hoc Tukey's HSD
emm_habitat <- emmeans(GLMM_model, ~ Habitat)
pairs(emm_habitat, adjust = 'tukey')
#### NOTE: Difference in π per MOTU is significant ** between 'Seaport' and 'Sand'
####       Difference in π per MOTU is significant ** between 'Seaport' and 'Rock'
####       Difference in π per MOTU is not significant between 'Rock' and 'Sand'

emm_area <- emmeans(GLMM_model, ~ Area)
pairs(emm_area, adjust = 'tukey')
#### NOTE: Difference in π per MOTU is significant * between 'West' and 'East'
####       Difference in π per MOTU is significant * between 'West' and 'Central'
####       Difference in π per MOTU is not significant between 'East' and 'Central'

### 4-4-4. Obtain the R²
r.squaredGLMM(GLMM_model)
#### NOTE: Marginal R² = 2.59%
####       Conditional R² = 34.37%


## 4-5. Plot ASV richness per MOTU vs Habitat
### 4-5-1. Extract marginal effects
emm_habitat_df <- as.data.frame(emmeans(GLMM_model, ~ Habitat, type = 'response'))

### 4-5-2. Plot the marginal effects
INTRA_MOTU_pi_Habitat_plot <- ggplot(emm_habitat_df,
                                     aes(x = Habitat,
                                         y = response,
                                         fill = Habitat)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, linewidth = 0.75, color = 'black') +
  geom_point(size = 5, shape = 21, show.legend = FALSE) +
  geom_signif(comparisons = list(c('Sand', 'Seaport'),
                                 c('Rock', 'Seaport')),
              annotations = c('***', '***'),
              y_position = c(0.0032, 0.0031),
              tip_length = 0, size = 1, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = wes_palette('Moonrise2')) +
  ylab('Model-adjusted value') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = 'grey95'))
INTRA_MOTU_pi_Habitat_plot

### Save in pdf
ggsave(INTRA_MOTU_pi_Habitat_plot,
       file = '../../outputs/INTRA_MOTU_pi_Habitat_plot.pdf',
       height = 8,
       width = 6,
       unit = 'cm')


## 5-6. Plot π per MOTU vs Area
### 5-6-1. Extract marginal effects
emm_area_df <- as.data.frame(emmeans(GLMM_model, ~ Area, type = 'response'))

### 5-6-2. Plot the marginal effects
INTRA_MOTU_pi_Area_plot <- ggplot(emm_area_df,
                                  aes(x = Area,
                                      y = response,
                                      fill = Area)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.2, linewidth = 0.75, color = 'black') +
  geom_point(size = 5, shape = 21, show.legend = FALSE) +
  geom_signif(comparisons = list(c('West', 'East'),
                                 c('West', 'Central')),
              annotations = c('*', '*'),
              y_position = c(0.00295, 0.0029),
              tip_length = 0, size = 1, textsize = 5, vjust = 0.5) +
  scale_fill_manual(values = wes_palette('AsteroidCity1')[c(1, 3, 4)]) +
  ylab('Model-adjusted value') +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_blank(),
        legend.position = 'none',
        panel.background = element_rect(fill = 'grey95'))
INTRA_MOTU_pi_Area_plot

### Save in pdf
ggsave(INTRA_MOTU_pi_Area_plot,
       file = '../../outputs/INTRA_MOTU_pi_Area_plot.pdf',
       height = 8,
       width = 6,
       unit = 'cm')
# ______________________________________________