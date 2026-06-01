# ---- 0. Load packages ----
library(tidyverse)
library(betapart) # dissimilarity partitioning
library(vegan) # beta-diversity analyses
library(cluster) # draw ellipsoids
library(wesanderson)
# __________________________



# ---- 1. β-diversity on all MOTUs ----
## 1-1. Import data
### 1-1-1. Presence/Absence matrix (Subset MOTUs)
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


## 1-2. Compute PCoA
## 1-2-1. Measure Jaccard distance between samples
jac_dist <- beta.pair(final_matrix_new, 'jaccard')

#### Save the Jaccard distance matrix
save(jac_dist, file = '../data/beta_SD.Rdata')

### 1-2-2. PCoA ordination
pcoa_jac <- cmdscale(jac_dist$beta.jac, k = dim(final_matrix_new)[1] - 1, eig = T)


## 1-3. Check influence of richness on patterns
### 1-3-1. Compare Jaccard and Raup-Crick distances
raup_dist <- raupcrick(final_matrix_new)
mantel(jac_dist$beta.jac, raup_dist)
#### NOTE: Matrices are correlated: 0.608

### 1-3-2.Compare ordinations
pcoa_raup <- cmdscale(raup_dist, k = dim(final_matrix_new)[1] - 1, eig = T)
protest(pcoa_jac, pcoa_raup)

#### NOTE: Ordinations are highly correlated: 0.8087
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns


## 1-4. Visualize PCoA ordination
### 1-4-1. Extract PCoA data
pcoa_df <- as.data.frame(pcoa_jac$points[, 1:2])
colnames(pcoa_df) <- c('PCoA1', 'PCoA2')
pcoa_df$sample <- rownames(pcoa_df)

#### Add metadata
pcoa_df <- pcoa_df %>%
  left_join(metadata_new %>%
              select(sample, Habitat, Area),
            by = c('sample' = 'sample'))

#### Record the percentage of variance explained by axes 1 and 2
eig_vals <- pcoa_jac$eig
var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])
perc1 <- round(100 * var_explained[1], 1)
perc2 <- round(100 * var_explained[2], 1)

### 1-4-2. Calculate habitat's centroid
#### Extract axes names
axes <- c('PCoA1', 'PCoA2')

#### Extract coordinates
seaport_coords <- subset(pcoa_df, Habitat == 'Seaport', select = axes)
sand_coords <- subset(pcoa_df, Habitat == 'Sand', select = axes)
rock_coords <- subset(pcoa_df, Habitat == 'Rock', select = axes)

#### Calculate mean across columns (= centroid) for each habitat
seaport_centroid <- colMeans(seaport_coords)
sand_centroid <- colMeans(sand_coords)
rock_centroid <- colMeans(rock_coords)

#### Merge habitats
centroids <- rbind(seaport_centroid,
                   sand_centroid,
                   rock_centroid) %>% 
  as.data.frame() %>% 
  mutate(Habitat = c('Seaport', 'Sand', 'Rock'))
colnames(centroids)[1:2] <- c('X1_centroid', 'X2_centroid')

#### Add the centroids to PCoA dataframe
pcoa_df <- pcoa_df %>%
  left_join(centroids, by = 'Habitat') %>% 
  mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
         Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))

### 1-4-3. Draw ellipsoids around habitats in the PCoA
ellipse_df <- pcoa_df %>%
  group_split(Habitat) %>% # we divide the dataframe into a list of dataframe (one per habitat)
  map_dfr(~ {hull <- ellipsoidhull(as.matrix(.x[, c("PCoA1", "PCoA2")])) # compute the ellipsoid hull for each dataframe
  pts <- predict(hull) # record the predicted points of the hull
  tibble(X1 = pts[,1],
         X2 = pts[,2],
         Habitat = unique(.x$Habitat))})

### 1-4-4. Plot PCoA
INTER_all_MOTU_PCoA_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2,
                                                fill = Habitat,
                                                color = Habitat,
                                                shape = Area)) +
  geom_point(size = 3, color = 'black') +
  geom_segment(aes(xend = X1_centroid,
                   yend = X2_centroid)) +
  geom_polygon(data = ellipse_df,
               aes(x = X1, y = X2,
                   fill = Habitat, group = Habitat),
               alpha = 0.4, color = 'black', inherit.aes = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
  scale_color_manual(values = wes_palette('Moonrise2')[1:3]) +
  scale_shape_manual(values = c(21, 22, 25)) +
  labs(x = paste0('PCoA 1 (', perc1, '%)'),
       y = paste0('PCoA 2 (', perc2, '%)'),
       fill = 'Habitat', shape = 'Area') +
  guides(fill = guide_legend(title.position = 'top',
                             title.hjust = 0.5, nrow = 1),
         shape = guide_legend(title.position = 'top',
                              title.hjust = 0.5, nrow = 1)) +
  xlim(c(-0.6, 0.6)) + ylim(c(-0.6, 0.6)) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(size = 15, colour = 'black'),
        panel.background = element_rect(fill = 'grey95'),
        legend.background = element_rect(color = 'black', fill = 'white'),
        legend.position = 'top',
        legend.title.align = 0.5,
        legend.title = element_text(face = 'bold'))
INTER_all_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_all_MOTU_PCoA_plot,
       file = '../../outputs/INTER_all_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')


## 1-5. PERMANOVA
### 1-5-1. Bootstrap PERMANOVAs to account for site replication
set.seed(2026)
nboot <- 999

#### Create an empty list
permanovas <- vector('list', nboot)

#### Perform 999 bootstrap iterations
for (i in seq_len(nboot)){
  
  ##### Randomly select one sample per site
  meta_sub <- metadata_new %>%
    group_by(pool) %>%
    slice_sample(n = 1) %>%
    ungroup()
  samples <- meta_sub %>% 
    pull() %>%
    as.vector()
  
  ##### Subset the distance matrix to these samples
  jac_sub <- as.dist(as.matrix(jac_dist$beta.jac)[samples, samples, drop = FALSE])
  
  ##### Compute PERMANOVA
  fit <- adonis2(jac_sub ~ Habitat + Area,
                 data = meta_sub,
                 permutations = 999,
                 by = 'margin')
  
  ##### Record the results
  permanovas[[i]] <- data.frame(term = rownames(fit),
                                R2 = fit$R2,
                                p = fit$`Pr(>F)`)
}

### 1-5-2. Check R² and p-value distributions
#### 1-5-2-1. Habitat
##### 1-5-2-1-1. R²
###### Extract values
r2_habitat <- sapply(permanovas, function(df){
  filter(df, term == 'Habitat')$R2
})
mean_r2_habitat <- mean(r2_habitat)

###### Save in pdf
pdf('../../outputs/INTER_all_MOTU_PCoA_R2_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_habitat, col = 'skyblue',
     main = 'Habitat R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_habitat, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_habitat, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_habitat, 3)), adj = 0)
dev.off()

##### 1-5-2-1-2. p-value
###### Extract values
p_habitat <- sapply(permanovas, function(df){
  filter(df, term == 'Habitat')$p
})
mean_p_habitat <- mean(p_habitat)
prop_habitat_0.05 <- mean(p_habitat <= 0.05)
prop_habitat_0.01 <- mean(p_habitat <= 0.01)
prop_habitat_0.001 <- mean(p_habitat <= 0.001)

###### Save in pdf
pdf('../../outputs/INTER_all_MOTU_PCoA_p_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(p_habitat, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_habitat, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_habitat, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_habitat_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_habitat_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_habitat_0.001, 3)), adj = 0)
dev.off()

###### NOTE: 35.7% of the variance is explained by the Habitat variable, and the effect is significant **

#### 1-5-2-2. Area
##### 1-5-2-2-1. R²
###### Extract values
r2_area <- sapply(permanovas, function(df){
  filter(df, term == 'Area')$R2
})
mean_r2_area <- mean(r2_area)

###### Save in pdf
pdf('../../outputs/INTER_all_MOTU_PCoA_R2_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_area, col = 'skyblue',
     main = 'Area R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_area, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_area, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_area, 3)), adj = 0)
dev.off()

##### 1-5-2-2-2. p-value
###### Extract values
p_area <- sapply(permanovas, function(df){
  filter(df, term == 'Area')$p
})
mean_p_area <- mean(p_area)
prop_area_0.05 <- mean(p_area <= 0.05)
prop_area_0.01 <- mean(p_area <= 0.01)
prop_area_0.001 <- mean(p_area <= 0.001)

###### Save in pdf
pdf('../../outputs/INTER_all_MOTU_PCoA_p_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(p_area, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Area p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_area, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_area, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_area_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_area_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_area_0.001, 3)), adj = 0)
dev.off()

###### NOTE: 37.4% of the variance is explained by the Area variable, and the effect is significant **


## 1-6. PERMDISP
### 1-6-1. Bootstrap PERMDISPs to account for site replication
set.seed(2026)
nboot <- 999

#### Create an empty list
permdisps <- vector('list', nboot)

#### Perform 999 bootstrap iterations
for (i in seq_len(nboot)){
  
  ##### Randomly select one sample per site
  meta_sub <- metadata_new %>%
    group_by(pool) %>%
    slice_sample(n = 1) %>%
    ungroup()
  samples <- meta_sub %>% 
    pull() %>%
    as.vector()
  
  ##### Subset the distance matrix to these samples
  jac_sub <- as.dist(as.matrix(jac_dist$beta.jac)[samples, samples, drop = FALSE])
  
  ##### Compute PERMDISP
  habitat_fit <- anova(betadisper(jac_sub, meta_sub$Habitat))
  area_fit <- anova(betadisper(jac_sub, meta_sub$Area))
  
  ##### Record the results
  permdisps[[i]] <- data.frame(p_habitat = habitat_fit$`Pr(>F)`[1],
                               p_area = area_fit$`Pr(>F)`[1])
}

### 1-6-2. Check p-value distributions
#### 1-6-2-1. Habitat
p_habitat <- sapply(permdisps, function(df){
  unlist(select(df, p_habitat))
})
mean_p_habitat <- mean(p_habitat)
prop_habitat_0.05 <- mean(p_habitat <= 0.05)
prop_habitat_0.01 <- mean(p_habitat <= 0.01)
prop_habitat_0.001 <- mean(p_habitat <= 0.001)

##### Save in pdf
pdf('../../outputs/INTER_all_MOTU_PCoA_p_Habitat_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_habitat, col = 'skyblue', xlim = c(0, 1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.1, y = max(hist(p_habitat, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_habitat, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_habitat_0.05, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_habitat_0.01, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_habitat_0.001, 3)), adj = 0)
dev.off()
##### NOTE: The effect of Habitat is not significant

#### 1-6-2-2. Area
p_area <- sapply(permdisps, function(df){
  unlist(select(df, p_area))
})
mean_p_area <- mean(p_area)
prop_area_0.05 <- mean(p_area <= 0.05)
prop_area_0.01 <- mean(p_area <= 0.01)
prop_area_0.001 <- mean(p_area <= 0.001)

##### Save in pdf
pdf('../../outputs/INTER_all_MOTU_PCoA_p_Area_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_area, col = 'skyblue', xlim = c(0, 1),
     main = 'Area p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.1, y = max(hist(p_area, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_area, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_area_0.05, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_area_0.01, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_area_0.001, 3)), adj = 0)
dev.off()
##### NOTE: The effect of Area is not significant
# _____________________________________



# ---- 2. β-diversity on MOTUs by kingdoms ----
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


## 2-2. Compute MOTU dissimilarity matrix per kingdom
kingdom_MOTU_jac <- function(selected_kingdom){
  
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
  
  ### Compute MOTU dissimilarity matrix
  jac_dist <- beta.pair(final_matrix_new, 'jaccard')
  
  ### PCoA ordination
  pcoa_jac <- cmdscale(jac_dist$beta.jac, k = dim(final_matrix_new)[1] - 1, eig = T)
  
  ### Check influence of richness on patterns
  #### Compare Jaccard and Raup-Crick distances
  raup_dist <- raupcrick(final_matrix_new)
  print(mantel(jac_dist$beta.jac, raup_dist))
  
  #### Compare ordinations
  pcoa_raup <- cmdscale(raup_dist, k = dim(final_matrix_new)[1] - 1, eig = T)
  print(protest(pcoa_jac, pcoa_raup))
  
  return(list(jac_dist = jac_dist,
              pcoa_jac = pcoa_jac))
}

### Metazoa
metazoa_beta <- kingdom_MOTU_jac('Metazoa')
#### NOTE: Matrices are correlated: 0.6378
#### NOTE: Ordinations are highly correlated: 0.7962
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns

### Chromista
chromista_beta <- kingdom_MOTU_jac('Chromista')
#### NOTE: Matrices are correlated: 0.7371
#### NOTE: Ordinations are highly correlated: 0.8657
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns

### Plantae
plantae_beta <- kingdom_MOTU_jac('Plantae')
#### NOTE: Matrices are correlated: 0.5464
#### NOTE: Ordinations are highly correlated: 0.8327
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns


## 2-3- Visualize PCoA per kingdom
kingdom_pcoa <- function(pcoa_jac){
  
  ### Extract PCoA data
  pcoa_df <- as.data.frame(pcoa_jac$points[, 1:2])
  colnames(pcoa_df) <- c('PCoA1', 'PCoA2')
  pcoa_df$sample <- rownames(pcoa_df)
  
  #### Add metadata
  pcoa_df <- pcoa_df %>%
    left_join(metadata_new %>%
                select(sample, Habitat, Area),
              by = c('sample' = 'sample'))
  
  #### Record the percentage of variance explained by axes 1 and 2
  eig_vals <- pcoa_jac$eig
  var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])
  perc1 <- round(100 * var_explained[1], 1)
  perc2 <- round(100 * var_explained[2], 1)
  
  ### Calculate habitat's centroid
  #### Extract axes names
  axes <- c('PCoA1', 'PCoA2')
  
  #### Extract coordinates
  seaport_coords <- subset(pcoa_df, Habitat == 'Seaport', select = axes)
  sand_coords <- subset(pcoa_df, Habitat == 'Sand', select = axes)
  rock_coords <- subset(pcoa_df, Habitat == 'Rock', select = axes)
  
  #### Calculate mean across columns (= centroid) for each habitat
  seaport_centroid <- colMeans(seaport_coords)
  sand_centroid <- colMeans(sand_coords)
  rock_centroid <- colMeans(rock_coords)
  
  #### Merge habitats
  centroids <- rbind(seaport_centroid,
                     sand_centroid,
                     rock_centroid) %>% 
    as.data.frame() %>% 
    mutate(Habitat = c('Seaport', 'Sand', 'Rock'))
  colnames(centroids)[1:2] <- c('X1_centroid', 'X2_centroid')
  
  #### Add the centroids to PCoA dataframe
  pcoa_df <- pcoa_df %>%
    left_join(centroids, by = 'Habitat') %>% 
    mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
           Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))
  
  ### Draw ellipsoids around habitats in the PCoA
  ellipse_df <- pcoa_df %>%
    group_split(Habitat) %>% # we divide the dataframe into a list of dataframe (one per habitat)
    map_dfr(~ {hull <- ellipsoidhull(as.matrix(.x[, c("PCoA1", "PCoA2")])) # compute the ellipsoid hull for each dataframe
    pts <- predict(hull) # record the predicted points of the hull
    tibble(X1 = pts[,1],
           X2 = pts[,2],
           Habitat = unique(.x$Habitat))})
  
  ### Plot PCoA
  INTER_kingdom_MOTU_PCoA_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2,
                                                      fill = Habitat,
                                                      color = Habitat,
                                                      shape = Area)) +
    geom_point(size = 3, color = 'black') +
    geom_segment(aes(xend = X1_centroid,
                     yend = X2_centroid)) +
    geom_polygon(data = ellipse_df,
                 aes(x = X1, y = X2,
                     fill = Habitat, group = Habitat),
                 alpha = 0.4, color = 'black', inherit.aes = FALSE) +
    scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
    scale_color_manual(values = wes_palette('Moonrise2')[1:3]) +
    scale_shape_manual(values = c(21, 22, 25)) +
    labs(x = paste0('PCoA 1 (', perc1, '%)'),
         y = paste0('PCoA 2 (', perc2, '%)'),
         fill = 'Habitat', shape = 'Area') +
    guides(fill = guide_legend(title.position = 'top',
                               title.hjust = 0.5, nrow = 1),
           shape = guide_legend(title.position = 'top',
                                title.hjust = 0.5, nrow = 1)) +
    xlim(c(-0.6, 0.6)) + ylim(c(-0.6, 0.6)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          axis.text = element_text(colour = 'black', size = 10),
          axis.title = element_text(size = 15, colour = 'black'),
          legend.background = element_rect(color = 'black', fill = 'white'),
          legend.position = 'top',
          legend.title.align = 0.5,
          legend.title = element_text(face = 'bold'))
  
  return(INTER_kingdom_MOTU_PCoA_plot)
}

### Metazoa
INTER_metazoa_MOTU_PCoA_plot <- kingdom_pcoa(metazoa_beta$pcoa_jac)
INTER_metazoa_MOTU_PCoA_plot
#### Save in pdf
ggsave(INTER_metazoa_MOTU_PCoA_plot,
       file = '../../outputs/INTER_metazoa_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Chromista
INTER_chromista_MOTU_PCoA_plot <- kingdom_pcoa(chromista_beta$pcoa_jac)
INTER_chromista_MOTU_PCoA_plot
#### Save in pdf
ggsave(INTER_chromista_MOTU_PCoA_plot,
       file = '../../outputs/INTER_chromista_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Plantae
INTER_plantae_MOTU_PCoA_plot <- kingdom_pcoa(plantae_beta$pcoa_jac)
INTER_plantae_MOTU_PCoA_plot
#### Save in pdf
ggsave(INTER_plantae_MOTU_PCoA_plot,
       file = '../../outputs/INTER_plantae_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')


## 2-4. Compute PERMANOVA per kingdom
kingdom_permanovas <- function(beta_jac){
  
  ### Bootstrap PERMANOVAs to account for site replication
  set.seed(2026)
  nboot <- 999
  
  #### Create an empty list
  permanovas <- vector('list', nboot)
  
  #### Perform 999 bootstrap iterations
  for (i in seq_len(nboot)){
    
    ##### Randomly select one sample per site
    meta_sub <- metadata_new %>%
      group_by(pool) %>%
      slice_sample(n = 1) %>%
      ungroup()
    samples <- meta_sub %>% 
      pull() %>%
      as.vector()
    
    ##### Subset the distance matrix to these samples
    jac_sub <- as.dist(as.matrix(beta_jac)[samples, samples, drop = FALSE])
    
    ##### Compute PERMANOVA
    fit <- adonis2(jac_sub ~ Habitat + Area,
                   data = meta_sub,
                   permutations = 999,
                   by = 'margin')
    
    ##### Record the results
    permanovas[[i]] <- data.frame(term = rownames(fit),
                                  R2 = fit$R2,
                                  p = fit$`Pr(>F)`)
  }
  
  return(permanovas)
}

### Metazoa
metazoa_permanovas <- kingdom_permanovas(metazoa_beta$jac_dist$beta.jac)

#### Check R² and p-value distributions
##### Habitat
###### R²
####### Extract values
r2_habitat_metazoa <- sapply(metazoa_permanovas, function(df){
  filter(df, term == 'Habitat')$R2
})
mean_r2_habitat_metazoa <- mean(r2_habitat_metazoa)

####### Save in pdf
pdf('../../outputs/INTER_metazoa_MOTU_PCoA_R2_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_habitat_metazoa, col = 'skyblue',
     main = 'Habitat R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_habitat_metazoa, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_habitat_metazoa, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_habitat_metazoa, 3)), adj = 0)
dev.off()

###### p-value
####### Extract values
p_habitat_metazoa <- sapply(metazoa_permanovas, function(df){
  filter(df, term == 'Habitat')$p
})
mean_p_habitat_metazoa <- mean(p_habitat_metazoa)
prop_habitat_0.05 <- mean(p_habitat_metazoa <= 0.05)
prop_habitat_0.01 <- mean(p_habitat_metazoa <= 0.01)
prop_habitat_0.001 <- mean(p_habitat_metazoa <= 0.001)

####### Save in pdf
pdf('../../outputs/INTER_metazoa_MOTU_PCoA_p_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(p_habitat_metazoa, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_habitat_metazoa, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_habitat_metazoa, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_metazoa, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_habitat_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_metazoa, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_habitat_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_metazoa, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_habitat_0.001, 3)), adj = 0)
dev.off()

####### NOTE: 37% of the variance is explained by the Habitat variable, and the effect is significant **

##### Area
###### R²
####### Extract values
r2_area_metazoa <- sapply(metazoa_permanovas, function(df){
  filter(df, term == 'Area')$R2
})
mean_r2_area_metazoa <- mean(r2_area_metazoa)

####### Save in pdf
pdf('../../outputs/INTER_metazoa_MOTU_PCoA_R2_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_area_metazoa, col = 'skyblue',
     main = 'Area R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_area_metazoa, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_area_metazoa, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_area_metazoa, 3)), adj = 0)
dev.off()

###### p-value
####### Extract values
p_area_metazoa <- sapply(metazoa_permanovas, function(df){
  filter(df, term == 'Area')$p
})
mean_p_area_metazoa <- mean(p_area_metazoa)
prop_area_0.05 <- mean(p_area_metazoa <= 0.05)
prop_area_0.01 <- mean(p_area_metazoa <= 0.01)
prop_area_0.001 <- mean(p_area_metazoa <= 0.001)

####### Save in pdf
pdf('../../outputs/INTER_metazoa_MOTU_PCoA_p_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(p_area_metazoa, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Area p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_area_metazoa, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_area_metazoa, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_metazoa, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_area_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_metazoa, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_area_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_metazoa, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_area_0.001, 3)), adj = 0)
dev.off()

####### NOTE: 34.8% of the variance is explained by the Area variable, and the effect is significant **

### Chromista
chromista_permanovas <- kingdom_permanovas(chromista_beta$jac_dist$beta.jac)

#### Check R² and p-value distributions
##### Habitat
###### R²
####### Extract values
r2_habitat_chromista <- sapply(chromista_permanovas, function(df){
  filter(df, term == 'Habitat')$R2
})
mean_r2_habitat_chromista <- mean(r2_habitat_chromista)

####### Save in pdf
pdf('../../outputs/INTER_chromista_MOTU_PCoA_R2_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_habitat_chromista, col = 'skyblue',
     main = 'Habitat R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_habitat_chromista, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_habitat_chromista, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_habitat_chromista, 3)), adj = 0)
dev.off()

###### p-value
####### Extract values
p_habitat_chromista <- sapply(chromista_permanovas, function(df){
  filter(df, term == 'Habitat')$p
})
mean_p_habitat_chromista <- mean(p_habitat_chromista)
prop_habitat_0.05 <- mean(p_habitat_chromista <= 0.05)
prop_habitat_0.01 <- mean(p_habitat_chromista <= 0.01)
prop_habitat_0.001 <- mean(p_habitat_chromista <= 0.001)

####### Save in pdf
pdf('../../outputs/INTER_chromista_MOTU_PCoA_p_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(p_habitat_chromista, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_habitat_chromista, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_habitat_chromista, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_chromista, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_habitat_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_chromista, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_habitat_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_chromista, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_habitat_0.001, 3)), adj = 0)
dev.off()

####### NOTE: 26.5% of the variance is explained by the Habitat variable, and the effect is significant *

##### Area
###### R²
####### Extract values
r2_area_chromista <- sapply(chromista_permanovas, function(df){
  filter(df, term == 'Area')$R2
})
mean_r2_area_chromista <- mean(r2_area_chromista)

####### Save in pdf
pdf('../../outputs/INTER_chromista_MOTU_PCoA_R2_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_area_chromista, col = 'skyblue',
     main = 'Area R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_area_chromista, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_area_chromista, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_area_chromista, 3)), adj = 0)
dev.off()

###### p-value
####### Extract values
p_area_chromista <- sapply(chromista_permanovas, function(df){
  filter(df, term == 'Area')$p
})
mean_p_area_chromista <- mean(p_area_chromista)
prop_area_0.05 <- mean(p_area_chromista <= 0.05)
prop_area_0.01 <- mean(p_area_chromista <= 0.01)
prop_area_0.001 <- mean(p_area_chromista <= 0.001)

####### Save in pdf
pdf('../../outputs/INTER_chromista_MOTU_PCoA_p_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(p_area_chromista, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Area p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_area_chromista, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_area_chromista, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_chromista, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_area_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_chromista, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_area_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_chromista, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_area_0.001, 3)), adj = 0)
dev.off()

####### NOTE: 46% of the variance is explained by the Area variable, and the effect is significant ***

### Plantae
plantae_permanovas <- kingdom_permanovas(plantae_beta$jac_dist$beta.jac)

#### Check R² and p-value distributions
##### Habitat
###### R²
####### Extract values
r2_habitat_plantae <- sapply(plantae_permanovas, function(df){
  filter(df, term == 'Habitat')$R2
})
mean_r2_habitat_plantae <- mean(r2_habitat_plantae)

####### Save in pdf
pdf('../../outputs/INTER_plantae_MOTU_PCoA_R2_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_habitat_plantae, col = 'skyblue',
     main = 'Habitat R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_habitat_plantae, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_habitat_plantae, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_habitat_plantae, 3)), adj = 0)
dev.off()

###### p-value
####### Extract values
p_habitat_plantae <- sapply(plantae_permanovas, function(df){
  filter(df, term == 'Habitat')$p
})
mean_p_habitat_plantae <- mean(p_habitat_plantae)
prop_habitat_0.05 <- mean(p_habitat_plantae <= 0.05)
prop_habitat_0.01 <- mean(p_habitat_plantae <= 0.01)
prop_habitat_0.001 <- mean(p_habitat_plantae <= 0.001)

####### Save in pdf
pdf('../../outputs/INTER_plantae_MOTU_PCoA_p_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(p_habitat_plantae, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_habitat_plantae, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_habitat_plantae, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_plantae, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_habitat_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_plantae, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_habitat_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_habitat_plantae, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_habitat_0.001, 3)), adj = 0)
dev.off()

####### NOTE: 45.4% of the variance is explained by the Habitat variable, and the effect is significant **

##### Area
###### R²
####### Extract values
r2_area_plantae <- sapply(plantae_permanovas, function(df){
  filter(df, term == 'Area')$R2
})
mean_r2_area_plantae <- mean(r2_area_plantae)

####### Save in pdf
pdf('../../outputs/INTER_plantae_MOTU_PCoA_R2_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_area_plantae, col = 'skyblue',
     main = 'Area R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_area_plantae, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_area_plantae, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_area_plantae, 3)), adj = 0)
dev.off()

###### p-value
####### Extract values
p_area_plantae <- sapply(plantae_permanovas, function(df){
  filter(df, term == 'Area')$p
})
mean_p_area_plantae <- mean(p_area_plantae)
prop_area_0.05 <- mean(p_area_plantae <= 0.05)
prop_area_0.01 <- mean(p_area_plantae <= 0.01)
prop_area_0.001 <- mean(p_area_plantae <= 0.001)

####### Save in pdf
pdf('../../outputs/INTER_plantae_MOTU_PCoA_p_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(p_area_plantae, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Area p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_area_plantae, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_area_plantae, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_plantae, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_area_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_plantae, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_area_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_area_plantae, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_area_0.001, 3)), adj = 0)
dev.off()

####### NOTE: 33.2% of the variance is explained by the Area variable, and the effect is significant *


## 2-5. Compute PERMDISP per kingdom
kingdom_permdisps <- function(beta_jac){
  
  ### Bootstrap PERMDISPs to account for site replication
  set.seed(2026)
  nboot <- 999
  
  #### Create an empty list
  permdisps <- vector('list', nboot)
  
  #### Perform 999 bootstrap iterations
  for (i in seq_len(nboot)){
    
    ##### Randomly select one sample per site
    meta_sub <- metadata_new %>%
      group_by(pool) %>%
      slice_sample(n = 1) %>%
      ungroup()
    samples <- meta_sub %>% 
      pull() %>%
      as.vector()
    
    ##### Subset the distance matrix to these samples
    jac_sub <- as.dist(as.matrix(beta_jac)[samples, samples, drop = FALSE])
    
    ##### Compute PERMDISP
    habitat_fit <- anova(betadisper(jac_sub, meta_sub$Habitat))
    area_fit <- anova(betadisper(jac_sub, meta_sub$Area))
    
    ##### Record the results
    permdisps[[i]] <- data.frame(p_habitat = habitat_fit$`Pr(>F)`[1],
                                 p_area = area_fit$`Pr(>F)`[1])
  }
  
  return(permdisps)
}

### Metazoa
metazoa_permdisps <- kingdom_permdisps(metazoa_beta$jac_dist$beta.jac)

### Check p-value distributions
#### Habitat
p_habitat_metazoa <- sapply(metazoa_permdisps, function(df){
  unlist(select(df, p_habitat))
})
mean_p_habitat_metazoa <- mean(p_habitat_metazoa)
prop_habitat_0.05 <- mean(p_habitat_metazoa <= 0.05)
prop_habitat_0.01 <- mean(p_habitat_metazoa <= 0.01)
prop_habitat_0.001 <- mean(p_habitat_metazoa <= 0.001)

##### Save in pdf
pdf('../../outputs/INTER_metazoa_MOTU_PCoA_p_Habitat_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_habitat_metazoa, col = 'skyblue', xlim = c(0, 1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.1, y = max(hist(p_habitat_metazoa, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_habitat_metazoa, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_metazoa, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_habitat_0.05, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_metazoa, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_habitat_0.01, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_metazoa, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_habitat_0.001, 3)), adj = 0)
dev.off()
##### NOTE: The effect of Habitat is not significant

#### Area
p_area_metazoa <- sapply(metazoa_permdisps, function(df){
  unlist(select(df, p_area))
})
mean_p_area_metazoa <- mean(p_area_metazoa)
prop_area_0.05 <- mean(p_area_metazoa <= 0.05)
prop_area_0.01 <- mean(p_area_metazoa <= 0.01)
prop_area_0.001 <- mean(p_area_metazoa <= 0.001)

##### Save in pdf
pdf('../../outputs/INTER_metazoa_MOTU_PCoA_p_Area_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_area_metazoa, col = 'skyblue', xlim = c(0, 1),
     main = 'Area p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.1, y = max(hist(p_area_metazoa, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_area_metazoa, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_metazoa, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_area_0.05, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_metazoa, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_area_0.01, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_metazoa, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_area_0.001, 3)), adj = 0)
dev.off()
##### NOTE: The effect of Area is not significant

### Chromista
chromista_permdisps <- kingdom_permdisps(chromista_beta$jac_dist$beta.jac)

### Check p-value distributions
#### Habitat
p_habitat_chromista <- sapply(chromista_permdisps, function(df){
  unlist(select(df, p_habitat))
})
mean_p_habitat_chromista <- mean(p_habitat_chromista)
prop_habitat_0.05 <- mean(p_habitat_chromista <= 0.05)
prop_habitat_0.01 <- mean(p_habitat_chromista <= 0.01)
prop_habitat_0.001 <- mean(p_habitat_chromista <= 0.001)

##### Save in pdf
pdf('../../outputs/INTER_chromista_MOTU_PCoA_p_Habitat_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_habitat_chromista, col = 'skyblue', xlim = c(0, 1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.1, y = max(hist(p_habitat_chromista, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_habitat_chromista, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_chromista, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_habitat_0.05, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_chromista, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_habitat_0.01, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_chromista, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_habitat_0.001, 3)), adj = 0)
dev.off()
##### NOTE: The effect of Habitat is not significant

#### Area
p_area_chromista <- sapply(chromista_permdisps, function(df){
  unlist(select(df, p_area))
})
mean_p_area_chromista <- mean(p_area_chromista)
prop_area_0.05 <- mean(p_area_chromista <= 0.05)
prop_area_0.01 <- mean(p_area_chromista <= 0.01)
prop_area_0.001 <- mean(p_area_chromista <= 0.001)

##### Save in pdf
pdf('../../outputs/INTER_chromista_MOTU_PCoA_p_Area_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_area_chromista, col = 'skyblue', xlim = c(0, 1),
     main = 'Area p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.1, y = max(hist(p_area_chromista, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_area_chromista, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_chromista, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_area_0.05, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_chromista, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_area_0.01, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_chromista, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_area_0.001, 3)), adj = 0)
dev.off()
##### NOTE: The effect of Area is not significant

### Plantae
plantae_permdisps <- kingdom_permdisps(plantae_beta$jac_dist$beta.jac)

### Check p-value distributions
#### Habitat
p_habitat_plantae <- sapply(plantae_permdisps, function(df){
  unlist(select(df, p_habitat))
})
mean_p_habitat_plantae <- mean(p_habitat_plantae)
prop_habitat_0.05 <- mean(p_habitat_plantae <= 0.05)
prop_habitat_0.01 <- mean(p_habitat_plantae <= 0.01)
prop_habitat_0.001 <- mean(p_habitat_plantae <= 0.001)

##### Save in pdf
pdf('../../outputs/INTER_plantae_MOTU_PCoA_p_Habitat_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_habitat_plantae, col = 'skyblue', xlim = c(0, 1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.1, y = max(hist(p_habitat_plantae, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_habitat_plantae, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_plantae, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_habitat_0.05, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_plantae, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_habitat_0.01, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_habitat_plantae, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_habitat_0.001, 3)), adj = 0)
dev.off()
##### NOTE: The effect of Habitat is not significant

#### Area
p_area_plantae <- sapply(plantae_permdisps, function(df){
  unlist(select(df, p_area))
})
mean_p_area_plantae <- mean(p_area_plantae)
prop_area_0.05 <- mean(p_area_plantae <= 0.05)
prop_area_0.01 <- mean(p_area_plantae <= 0.01)
prop_area_0.001 <- mean(p_area_plantae <= 0.001)

##### Save in pdf
pdf('../../outputs/INTER_plantae_MOTU_PCoA_p_Area_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_area_plantae, col = 'skyblue', xlim = c(0, 1),
     main = 'Area p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.1, y = max(hist(p_area_plantae, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_area_plantae, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_plantae, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_area_0.05, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_plantae, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_area_0.01, 3)), adj = 0)
text(x = 0.1, y = max(hist(p_area_plantae, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_area_0.001, 3)), adj = 0)
dev.off()
##### NOTE: The effect of Area is not significant
# _____________________________________________



# ---- 3. β-diversity on MOTUs by phyla ----
## 3-1. Select the phyla
MOTUs_by_phylum <- MOTU_table_final %>%
  filter(!is.na(phylum)) %>%
  group_by(phylum) %>%
  summarise(n_MOTUs = n(), .groups = 'drop')
MOTUs_by_phylum

### Restrict to phyla with >= 10 MOTUs
kept_phylum <- MOTUs_by_phylum %>% 
  filter(n_MOTUs >= 10) %>% 
  select(phylum) %>% 
  unlist()
kept_phylum


## 3-2. Compute MOTU dissimilarity matrix per phylum
phylum_MOTU_jac <- function(selected_phylum){
  
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
  
  ### Compute MOTU dissimilarity matrix
  jac_dist <- beta.pair(final_matrix_new, 'jaccard')
  
  ### PCoA ordination
  pcoa_jac <- cmdscale(jac_dist$beta.jac, k = dim(final_matrix_new)[1] - 1, eig = T)
  
  ### Check influence of richness on patterns
  #### Compare Jaccard and Raup-Crick distances
  raup_dist <- raupcrick(final_matrix_new)
  print(mantel(jac_dist$beta.jac, raup_dist))
  
  #### Compare ordinations
  pcoa_raup <- cmdscale(raup_dist, k = dim(final_matrix_new)[1] - 1, eig = T)
  print(protest(pcoa_jac, pcoa_raup))
  
  return(list(jac_dist = jac_dist,
              pcoa_jac = pcoa_jac))
}

### Annelida
annelida_beta <- phylum_MOTU_jac('Annelida')
#### NOTE: Matrices are correlated: 0.8175
#### NOTE: Ordinations are highly correlated: 0.8803
#### NOTE: The effect of richness is negligible on Jaccard dissimilarity patterns

### Arthropoda
arthropoda_beta <- phylum_MOTU_jac('Arthropoda')
#### NOTE: Matrices are correlated: 0.6029
#### NOTE: Ordinations are highly correlated: 0.8158
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns

### Chordata
chordata_beta <- phylum_MOTU_jac('Chordata')
#### NOTE: Matrices are correlated: 0.8551
#### NOTE: Ordinations are highly correlated: 0.9395
#### NOTE: The effect of richness is negligible on Jaccard dissimilarity patterns

### Cnidaria
cnidaria_beta <- phylum_MOTU_jac('Cnidaria')
#### NOTE: Matrices are correlated: 0.6108
#### NOTE: Ordinations are highly correlated: 0.8902
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns

### Mollusca
mollusca_beta <- phylum_MOTU_jac('Mollusca')
#### NOTE: Matrices are correlated: 0.6815
#### NOTE: Ordinations are highly correlated: 0.8228
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns

### Nemertea
nemertea_beta <- phylum_MOTU_jac('Nemertea')
#### NOTE: Matrices are correlated: 0.9427
#### NOTE: Ordinations are highly correlated: 0.945
#### NOTE: The effect of richness is negligible on Jaccard dissimilarity patterns

### Bacillariophyta
bacillariophyta_beta <- phylum_MOTU_jac('Bacillariophyta')
#### NOTE: Matrices are correlated: 0.7175
#### NOTE: Ordinations are highly correlated: 0.8684
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns

### Ochrophyta
ochrophyta_beta <- phylum_MOTU_jac('Ochrophyta')
#### NOTE: Matrices are correlated: 0.7849
#### NOTE: Ordinations are highly correlated: 0.8734
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns

### Oomycota
oomycota_beta <- phylum_MOTU_jac('Oomycota')
#### NOTE: Matrices are correlated: 0.8539
#### NOTE: Ordinations are highly correlated: 0.9267
#### NOTE: The effect of richness is negligible on Jaccard dissimilarity patterns

### Rhodophyta
rhodophyta_beta <- phylum_MOTU_jac('Rhodophyta')
#### NOTE: Matrices are correlated: 0.6187
#### NOTE: Ordinations are highly correlated: 0.8602
#### NOTE: The effect of richness is moderate on Jaccard dissimilarity patterns

## 3-3. Visualize PCoA per phylum
phylum_pcoa <- function(pcoa_jac){
  
  ### Extract PCoA data
  pcoa_df <- as.data.frame(pcoa_jac$points[, 1:2])
  colnames(pcoa_df) <- c('PCoA1', 'PCoA2')
  pcoa_df$sample <- rownames(pcoa_df)
  
  #### Add metadata
  pcoa_df <- pcoa_df %>%
    left_join(metadata_new %>%
                select(sample, Habitat, Area),
              by = c('sample' = 'sample'))
  
  #### Record the percentage of variance explained by axes 1 and 2
  eig_vals <- pcoa_jac$eig
  var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])
  perc1 <- round(100 * var_explained[1], 1)
  perc2 <- round(100 * var_explained[2], 1)
  
  ### Calculate habitat's centroid
  #### Extract axes names
  axes <- c('PCoA1', 'PCoA2')
  
  #### Extract coordinates
  seaport_coords <- subset(pcoa_df, Habitat == 'Seaport', select = axes)
  sand_coords <- subset(pcoa_df, Habitat == 'Sand', select = axes)
  rock_coords <- subset(pcoa_df, Habitat == 'Rock', select = axes)
  
  #### Calculate mean across columns (= centroid) for each habitat
  seaport_centroid <- colMeans(seaport_coords)
  sand_centroid <- colMeans(sand_coords)
  rock_centroid <- colMeans(rock_coords)
  
  #### Merge habitats
  centroids <- rbind(seaport_centroid,
                     sand_centroid,
                     rock_centroid) %>% 
    as.data.frame() %>% 
    mutate(Habitat = c('Seaport', 'Sand', 'Rock'))
  colnames(centroids)[1:2] <- c('X1_centroid', 'X2_centroid')
  
  #### Add the centroids to PCoA dataframe
  pcoa_df <- pcoa_df %>%
    left_join(centroids, by = 'Habitat') %>% 
    mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
           Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))
  
  ### Draw ellipsoids around habitats in the PCoA
  ellipse_df <- pcoa_df %>%
    group_split(Habitat) %>% # we divide the dataframe into a list of dataframe (one per habitat)
    map_dfr(~ {hull <- ellipsoidhull(as.matrix(.x[, c("PCoA1", "PCoA2")])) # compute the ellipsoid hull for each dataframe
    pts <- predict(hull) # record the predicted points of the hull
    tibble(X1 = pts[,1],
           X2 = pts[,2],
           Habitat = unique(.x$Habitat))})
  
  ### Plot PCoA
  INTER_phylum_MOTU_PCoA_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2,
                                                      fill = Habitat,
                                                      color = Habitat,
                                                      shape = Area)) +
    geom_point(size = 3, color = 'black') +
    geom_segment(aes(xend = X1_centroid,
                     yend = X2_centroid)) +
    # geom_polygon(data = ellipse_df,
    #              aes(x = X1, y = X2,
    #                  fill = Habitat, group = Habitat),
    #              alpha = 0.4, color = 'black', inherit.aes = FALSE) +
    scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
    scale_color_manual(values = wes_palette('Moonrise2')[1:3]) +
    scale_shape_manual(values = c(21, 22, 25)) +
    labs(x = paste0('PCoA 1 (', perc1, '%)'),
         y = paste0('PCoA 2 (', perc2, '%)'),
         fill = 'Habitat', shape = 'Area') +
    guides(fill = guide_legend(title.position = 'top',
                               title.hjust = 0.5, nrow = 1),
           shape = guide_legend(title.position = 'top',
                                title.hjust = 0.5, nrow = 1)) +
    theme_classic() +
    theme(panel.border = element_rect(fill = NA),
          axis.text = element_text(colour = 'black', size = 10),
          axis.title = element_text(size = 15, colour = 'black'),
          legend.background = element_rect(color = 'black', fill = 'white'),
          legend.position = 'top',
          legend.title.align = 0.5,
          legend.title = element_text(face = 'bold'))
  
  return(INTER_phylum_MOTU_PCoA_plot)
}

### Annelida
INTER_annelida_MOTU_PCoA_plot <- phylum_pcoa(annelida_beta$pcoa_jac)
INTER_annelida_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_annelida_MOTU_PCoA_plot,
       file = '../../outputs/INTER_annelida_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Arthropoda
INTER_arthropoda_MOTU_PCoA_plot <- phylum_pcoa(arthropoda_beta$pcoa_jac)
INTER_arthropoda_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_arthropoda_MOTU_PCoA_plot,
       file = '../../outputs/INTER_arthropoda_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Chordata
INTER_chordata_MOTU_PCoA_plot <- phylum_pcoa(chordata_beta$pcoa_jac)
INTER_chordata_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_chordata_MOTU_PCoA_plot,
       file = '../../outputs/INTER_chordata_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Cnidaria
INTER_cnidaria_MOTU_PCoA_plot <- phylum_pcoa(cnidaria_beta$pcoa_jac)
INTER_cnidaria_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_cnidaria_MOTU_PCoA_plot,
       file = '../../outputs/INTER_cnidaria_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Mollusca
INTER_mollusca_MOTU_PCoA_plot <- phylum_pcoa(mollusca_beta$pcoa_jac)
INTER_mollusca_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_mollusca_MOTU_PCoA_plot,
       file = '../../outputs/INTER_mollusca_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Nemertea
INTER_nemertea_MOTU_PCoA_plot <- phylum_pcoa(nemertea_beta$pcoa_jac)
INTER_nemertea_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_nemertea_MOTU_PCoA_plot,
       file = '../../outputs/INTER_nemertea_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Bacillariophyta
INTER_bacillariophyta_MOTU_PCoA_plot <- phylum_pcoa(bacillariophyta_beta$pcoa_jac)
INTER_bacillariophyta_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_bacillariophyta_MOTU_PCoA_plot,
       file = '../../outputs/INTER_bacillariophyta_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Ochrophyta
INTER_ochrophyta_MOTU_PCoA_plot <- phylum_pcoa(ochrophyta_beta$pcoa_jac)
INTER_ochrophyta_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_ochrophyta_MOTU_PCoA_plot,
       file = '../../outputs/INTER_ochrophyta_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Oomycota
INTER_oomycota_MOTU_PCoA_plot <- phylum_pcoa(oomycota_beta$pcoa_jac)
INTER_oomycota_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_oomycota_MOTU_PCoA_plot,
       file = '../../outputs/INTER_oomycota_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')

### Rhodophyta
INTER_rhodophyta_MOTU_PCoA_plot <- phylum_pcoa(rhodophyta_beta$pcoa_jac)
INTER_rhodophyta_MOTU_PCoA_plot

#### Save in pdf
ggsave(INTER_rhodophyta_MOTU_PCoA_plot,
       file = '../../outputs/INTER_rhodophyta_MOTU_PCoA_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')
# __________________________________________