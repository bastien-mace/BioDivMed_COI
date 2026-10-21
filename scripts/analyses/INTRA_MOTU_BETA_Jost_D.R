# ---- 0. Load packages ----
library(tidyverse)
library(cluster) # draw ellipsoids
library(wesanderson)
library(vegan) # beta-diversity analyses
# __________________________



# ---- 1. Import data ----
## 1-1. Presence/Absence matrix (Subset MOTUs)
### Restrict MOTU table to samples' columns
MOTU_table_subset <- read.csv('../data/MOTU_table_subset.csv')
pa_matrix <- MOTU_table_subset %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

### Transform into presence-absence_matrix
pa_matrix[pa_matrix > 0] <- 1
colSums(pa_matrix)
final_matrix <- t(pa_matrix) # transpose matrix to get samples in rows and MOTUs in columns


## 1-2. ASV table
ASV_table_subset <- read.csv('../data/ASV_table_final.csv') %>% 
  filter(MOTU %in% MOTU_table_subset$id)


## 1-3. Site metadata
metadata <- read.csv2('../data/site_metadata.csv') %>%
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



# ---- 2. Compute Jost's D ----
Jost_D <- function(MOTU_pa_table, ASV_table, sample_pattern = '^SPY') {
  ## Add an id column to MOTU table
  MOTU_pa_table <- t(MOTU_pa_table) %>%
    as.data.frame() %>% 
    rownames_to_column('id')
  
  ## Create an empty mean Jost's D matrix (Sample x Sample)
  samples <- grep(sample_pattern, colnames(MOTU_pa_table), value = TRUE)
  n_samples <- length(samples)
  mean_D_matrix <- matrix(NA, n_samples, n_samples)
  rownames(mean_D_matrix) <- samples
  colnames(mean_D_matrix) <- samples
  
  ## Create an empty matrix with the number of MOTUs compared
  nb_matrix <- matrix(NA, n_samples, n_samples)
  rownames(nb_matrix) <- samples
  colnames(nb_matrix) <- samples
  
  ## Create a list to store all MOTU distance matrices
  MOTU_D_matrices <- list()
  
  ## Loop over sample pairs to calculate Jost's D
  for(i in 1:(n_samples - 1)){
    for(j in (i + 1):n_samples){
      s1 <- samples[i]
      s2 <- samples[j]
      
      ### Filter MOTUs present in both samples
      shared_MOTUs <- MOTU_pa_table %>%
        filter(.data[[s1]] > 0 & .data[[s2]] > 0) %>%
        pull(id)
      
      if(length(shared_MOTUs) == 0) next # stop here if no MOTUs are shared
      
      ### Record Jost's D values for the sample pair
      D_vals <- c()
      
      ### Extract ASV belonging to these MOTUs (MOTU by MOTU)
      for(m in shared_MOTUs){
        ASV_table_shared <- ASV_table %>% 
          filter(MOTU == m) %>% 
          select(id, starts_with('SPY')) %>% 
          column_to_rownames('id')
        ASV_table_shared <- t(ASV_table_shared) 
        
        #### Pool the problematic sample
        pooled_row <- matrix(as.numeric(ceiling(colSums(ASV_table_shared[samples_to_pool, , drop = FALSE]))),
                             nrow = 1,
                             dimnames = list(pooled_name, colnames(ASV_table_shared)))
        final_ASV_table_shared <- ASV_table_shared[!rownames(ASV_table_shared) %in% samples_to_pool, , drop = FALSE]
        final_ASV_table_shared <- rbind(final_ASV_table_shared, pooled_row)
        
        if(nrow(final_ASV_table_shared) < 2) next # stop here if the MOTU is not polymorphic
        
        #### Restrict the matrix to samples s1 and s2
        mat <- final_ASV_table_shared[c(s1, s2), , drop = FALSE]
        mat <- t(mat)
        
        if(sum(mat[, 1] > 0) < 2) next # stop here if there are less than 2 ASVs in sample s1
        
        if(sum(mat[, 2] > 0) < 2) next # stop here if there are less than 2 ASVs in sample s2
        
        #### Calculate haplotype frequencies for MOTU m
        p1 <- mat[, 1] / sum(mat[, 1]) # frequency in sample s1
        p2 <- mat[, 2] / sum(mat[, 2]) # frequency in sample s2
        
        #### Calculate haplotype diversities (Hd) for MOTU m
        Hd1 <- 1 - sum(p1^2) # Hd in sample s1
        Hd2 <- 1 - sum(p2^2) # Hd in sample s2
        
        #### Calculate mean haplotype diversity (Hs)
        Hs <- mean(c(Hd1, Hd2))
        
        #### Calculate total haplotype diversity (Ht)
        ptot <- (p1 + p2) / 2
        Ht <- 1 - sum(ptot^2)
        
        #### Calculate Jost's D
        D_m <- (Ht - Hs) / (1 - Hs)
        
        #### Add the Jost's D value to the sample pair vector
        D_vals <- c(D_vals, D_m)
        
        #### Initialize MOTU matrix if first time seeing this MOTU
        if(!m %in% names(MOTU_D_matrices)){
          MOTU_D_matrices[[m]] <- matrix(NA, n_samples, n_samples)
          rownames(MOTU_D_matrices[[m]]) <- samples
          colnames(MOTU_D_matrices[[m]]) <- samples
          diag(MOTU_D_matrices[[m]]) <- 0
        }
        
        #### Store Jost's D for this MOTU and sample pair
        MOTU_D_matrices[[m]][s1, s2] <- D_m
        MOTU_D_matrices[[m]][s2, s1] <- D_m
      }
      
      if(length(D_vals) > 0){
        #### Calculate mean Jost's D
        mean_D_matrix[s1, s2] <- mean(D_vals)
        mean_D_matrix[s2, s1] <- mean(D_vals)
        
        ### Record number of MOTUs compared
        nb_matrix[s1, s2] <- length(D_vals)
        nb_matrix[s2, s1] <- length(D_vals)
      }
    }
  }
  
  ## Convert MOTU matrices to dist objects
  MOTU_D_matrices <- lapply(MOTU_D_matrices, as.dist)
  
  return(list(mean_D_matrix = as.dist(mean_D_matrix),
              nb_matrix = as.dist(nb_matrix),
              MOTU_D_matrices = MOTU_D_matrices))
}
D <- Jost_D(final_matrix_new, ASV_table_subset)

### Save the Jost's D distance matrix
# save(D, file = '../data/beta_GD.Rdata')
# _____________________________



# ---- 3. Visualize PCoA ordination ----
## 3-1. Perform PCoA
pcoa <- cmdscale(D$mean_D_matrix, eig = TRUE)

### Extract points coordinates in a dataframe
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c('PCoA1', 'PCoA2')
pcoa_df$sample <- rownames(pcoa_df)

### Add metadata
pcoa_df <- pcoa_df %>%
  left_join(metadata_new %>%
              select(sample, Habitat, Area),
            by = c('sample' = 'sample'))

### Record the percentage of variance explained by axes 1 and 2
eig_vals <- pcoa$eig
var_explained <- eig_vals / sum(eig_vals[eig_vals > 0])
perc1 <- round(100 * var_explained[1], 1)
perc2 <- round(100 * var_explained[2], 1)


## 3-2. Calculate habitat's centroid
### Extract axes names
axes <- c('PCoA1', 'PCoA2')

### Extract coordinates
seaport_coords <- subset(pcoa_df, Habitat == 'Seaport', select = axes)
sand_coords <- subset(pcoa_df, Habitat == 'Sand', select = axes)
rock_coords <- subset(pcoa_df, Habitat == 'Rock', select = axes)

### Calculate mean across columns (= centroid) for each habitat
seaport_centroid <- colMeans(seaport_coords)
sand_centroid <- colMeans(sand_coords)
rock_centroid <- colMeans(rock_coords)

### Merge habitats
centroids <- rbind(seaport_centroid,
                   sand_centroid,
                   rock_centroid) %>% 
  as.data.frame() %>% 
  mutate(Habitat = c('Seaport', 'Sand', 'Rock'))
colnames(centroids)[1:2] <- c('X1_centroid', 'X2_centroid')

### Add the centroids to PCoA dataframe
pcoa_df <- pcoa_df %>%
  left_join(centroids, by = 'Habitat') %>% 
  mutate(Area = factor(Area, levels = c('West', 'Central', 'East')),
         Habitat = factor(Habitat, levels = c('Seaport', 'Rock', 'Sand')))


## 3-3. Draw ellipsoids around habitats in the PCoA
ellipse_df <- pcoa_df %>%
  group_split(Habitat) %>% # we divide the dataframe into a list of dataframe (one per habitat)
  map_dfr(~ {hull <- ellipsoidhull(as.matrix(.x[, c("PCoA1", "PCoA2")])) # compute the ellipsoid hull for each dataframe
  pts <- predict(hull) # record the predicted points of the hull
  tibble(X1 = pts[,1],
         X2 = pts[,2],
         Habitat = unique(.x$Habitat))})


## 3-4. Plot PCoA
INTRA_MOTU_PCoA_Jost_plot <- ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2,
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
  ylim( -0.1, 0.15) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA),
        axis.text = element_text(colour = 'black', size = 10),
        axis.title = element_text(size = 15, colour = 'black'),
        legend.background = element_rect(color = 'black', fill = 'white'),
        legend.position = 'top',
        legend.title.align = 0.5,
        legend.title = element_text(face = 'bold'))
INTRA_MOTU_PCoA_Jost_plot

### Save in pdf
ggsave(INTRA_MOTU_PCoA_Jost_plot,
       file = '../outputs/INTRA_MOTU_PCoA_Jost_plot.pdf',
       height = 12,
       width = 12,
       unit = 'cm')
# ______________________________________



# ---- 4. PERMANOVA ----
## 4-1. Bootstrap PERMANOVAs to account for site replication
set.seed(2026)
nboot <- 999

### Create an empty list
permanovas <- vector('list', nboot)

### Perform 999 bootstrap iterations
for (i in seq_len(nboot)){
  
  #### Randomly select one sample per site
  meta_sub <- metadata_new %>%
    group_by(pool) %>%
    slice_sample(n = 1) %>%
    ungroup()
  samples <- meta_sub %>% 
    pull() %>%
    as.vector()
  
  #### Subset the distance matrix to these samples
  D_sub <- as.dist(as.matrix(D$mean_D_matrix)[samples, samples, drop = FALSE])
  
  #### Compute PERMANOVA
  fit <- adonis2(D_sub ~ Habitat + Area,
                 data = meta_sub,
                 by = 'margin')
  
  #### Record the results
  permanovas[[i]] <- data.frame(term = rownames(fit),
                                R2 = fit$R2,
                                p = fit$`Pr(>F)`)
}


## 4-2. Check R² and p-value distributions
### 4-2-1. Habitat
#### 4-2-1-1. R²
##### Extract values
r2_habitat <- sapply(permanovas, function(df){
  filter(df, term == 'Habitat')$R2
})
mean_r2_habitat <- mean(r2_habitat)

##### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_R2_Habitat_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_habitat, col = 'skyblue',
     main = 'Habitat R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_habitat, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_habitat, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_habitat, 3)), adj = 0)
dev.off()

#### 4-2-1-2. p-value
##### Extract values
p_habitat <- sapply(permanovas, function(df){
  filter(df, term == 'Habitat')$p
})
mean_p_habitat <- mean(p_habitat)
prop_habitat_0.05 <- mean(p_habitat <= 0.05)
prop_habitat_0.01 <- mean(p_habitat <= 0.01)
prop_habitat_0.001 <- mean(p_habitat <= 0.001)

##### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_p_Habitat_PERMANOVA.pdf',
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

##### NOTE: 56.2% of the variance is explained by the Habitat variable, and the effect is significant ***

### 4-2-2. Area
#### 4-2-2-1. R²
##### Extract values
r2_area <- sapply(permanovas, function(df){
  filter(df, term == 'Area')$R2
})
mean_r2_area <- mean(r2_area)

##### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_R2_Area_PERMANOVA.pdf',
    height = 4,
    width = 6)
hist(r2_area, col = 'skyblue',
     main = 'Area R² bootstrap distribution', xlab = 'Bootstrapped R²')
text(x = max(hist(r2_area, plot = FALSE)$breaks) * 0.9, y = max(hist(r2_area, plot = FALSE)$counts),
     labels = paste0('Mean R²: ', round(mean_r2_area, 3)), adj = 0)
dev.off()

#### 4-2-2-2. p-value
##### Extract values
p_area <- sapply(permanovas, function(df){
  filter(df, term == 'Area')$p
})
mean_p_area <- mean(p_area)
prop_area_0.05 <- mean(p_area <= 0.05)
prop_area_0.01 <- mean(p_area <= 0.01)
prop_area_0.001 <- mean(p_area <= 0.001)

##### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_p_Area_PERMANOVA.pdf',
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

##### NOTE: 32.4% of the variance is explained by the Area variable, and the effect is significant ***
# ______________________



# ---- 5. PERMDISP ----
## 5-1. Bootstrap PERMDISPs to account for site replication
set.seed(2026)
nboot <- 999

### Create an empty list
permdisps <- vector('list', nboot)

### Perform 999 bootstrap iterations
for (i in seq_len(nboot)){
  
  #### Randomly select one sample per site
  meta_sub <- metadata_new %>%
    group_by(pool) %>%
    slice_sample(n = 1) %>%
    ungroup()
  samples <- meta_sub %>% 
    pull() %>%
    as.vector()
  
  #### Subset the distance matrix to these samples
  D_sub <- as.dist(as.matrix(D$mean_D_matrix)[samples, samples, drop = FALSE])
  
  #### Compute PERMDISP
  habitat_fit <- anova(betadisper(D_sub, meta_sub$Habitat))
  area_fit <- anova(betadisper(D_sub, meta_sub$Area))
  
  #### Record the results
  permdisps[[i]] <- data.frame(p_habitat = habitat_fit$`Pr(>F)`[1],
                               p_area = area_fit$`Pr(>F)`[1])
}


## 5-2. Check p-value distributions
### 5-2-1. Habitat
p_habitat <- sapply(permdisps, function(df){
  unlist(select(df, p_habitat))
})
mean_p_habitat <- mean(p_habitat)
prop_habitat_0.05 <- mean(p_habitat <= 0.05)
prop_habitat_0.01 <- mean(p_habitat <= 0.01)
prop_habitat_0.001 <- mean(p_habitat <= 0.001)

#### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_p_Habitat_PERMDISP.pdf',
    height = 4,
    width = 6)
hist(p_habitat, col = 'skyblue', xlim = c(0, 1),
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
#### NOTE: The effect of Habitat is significant *

### 5-2-2. Area
p_area <- sapply(permdisps, function(df){
  unlist(select(df, p_area))
})
mean_p_area <- mean(p_area)
prop_area_0.05 <- mean(p_area <= 0.05)
prop_area_0.01 <- mean(p_area <= 0.01)
prop_area_0.001 <- mean(p_area <= 0.001)

#### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_p_Area_PERMDISP.pdf',
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
#### NOTE: The effect of Area is not significant


## 5-3. Identify if the significant pattern is explained by one specific habitat
### 5-3-1. Compute pairwise post-hoc tests for Habitat
set.seed(2026)
nboot <- 999

### Create an empty list
pairwise_permdisps <- vector('list', nboot)

### Perform 999 bootstrap iterations
for (i in seq_len(nboot)){
  
  #### Randomly select one sample per site
  meta_sub <- metadata_new %>%
    group_by(pool) %>%
    slice_sample(n = 1) %>%
    ungroup()
  samples <- meta_sub %>% 
    pull() %>%
    as.vector()
  
  #### Subset the distance matrix to these samples
  D_sub <- as.dist(as.matrix(D$mean_D_matrix)[samples, samples, drop = FALSE])
  
  #### Compute PERMDISP with post-hoc Tukey HSD test
  habitat_pairwise <- TukeyHSD(betadisper(D_sub, meta_sub$Habitat))
  
  #### Record the results
  pairwise_permdisps[[i]] <- data.frame(p_Sand_Rock = habitat_pairwise$group['Sand-Rock', 'p adj'],
                                        p_Seaport_Rock = habitat_pairwise$group['Seaport-Rock', 'p adj'],
                                        p_Seaport_Sand = habitat_pairwise$group['Seaport-Sand', 'p adj'])
}

### 5-3-2. Check p-value distributions
#### 5-3-2-1. Sand vs Rock
p_Sand_Rock <- sapply(pairwise_permdisps, function(df){
  unlist(select(df, p_Sand_Rock))
})
mean_p_Sand_Rock <- mean(p_Sand_Rock)
prop_Sand_Rock_0.05 <- mean(p_Sand_Rock <= 0.05)
prop_Sand_Rock_0.01 <- mean(p_Sand_Rock <= 0.01)
prop_Sand_Rock_0.001 <- mean(p_Sand_Rock <= 0.001)

#### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_p_Sand_Rock_Tukey.pdf',
    height = 4,
    width = 6)
hist(p_Sand_Rock, col = 'skyblue', xlim = c(0, 1),
     main = 'Sand-Rock p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_Sand_Rock, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_Sand_Rock, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Sand_Rock, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_Sand_Rock_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Sand_Rock, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_Sand_Rock_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Sand_Rock, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_Sand_Rock_0.001, 3)), adj = 0)
dev.off()
#### NOTE: The difference between Sand and Rock is not significant

#### 5-3-2-2. Seaport vs Rock
p_Seaport_Rock <- sapply(pairwise_permdisps, function(df){
  unlist(select(df, p_Seaport_Rock))
})
mean_p_Seaport_Rock <- mean(p_Seaport_Rock)
prop_Seaport_Rock_0.05 <- mean(p_Seaport_Rock <= 0.05)
prop_Seaport_Rock_0.01 <- mean(p_Seaport_Rock <= 0.01)
prop_Seaport_Rock_0.001 <- mean(p_Seaport_Rock <= 0.001)

#### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_p_Seaport_Rock_Tukey.pdf',
    height = 4,
    width = 6)
hist(p_Seaport_Rock, col = 'skyblue', xlim = c(0, 1),
     main = 'Seaport-Rock p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_Seaport_Rock, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_Seaport_Rock, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Seaport_Rock, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_Seaport_Rock_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Seaport_Rock, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_Seaport_Rock_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Seaport_Rock, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_Seaport_Rock_0.001, 3)), adj = 0)
dev.off()
#### NOTE: The difference between Seaport and Rock is not significant

#### 5-3-2-3. Seaport vs Sand
p_Seaport_Sand <- sapply(pairwise_permdisps, function(df){
  unlist(select(df, p_Seaport_Sand))
})
mean_p_Seaport_Sand <- mean(p_Seaport_Sand)
prop_Seaport_Sand_0.05 <- mean(p_Seaport_Sand <= 0.05)
prop_Seaport_Sand_0.01 <- mean(p_Seaport_Sand <= 0.01)
prop_Seaport_Sand_0.001 <- mean(p_Seaport_Sand <= 0.001)

#### Save in pdf
pdf('../outputs/INTRA_MOTU_PCoA_Jost_p_Seaport_Sand_Tukey.pdf',
    height = 4,
    width = 6)
hist(p_Seaport_Sand, col = 'skyblue', xlim = c(0, 1),
     main = 'Seaport-Sand p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(p_Seaport_Sand, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p_Seaport_Sand, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Seaport_Sand, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_Seaport_Sand_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Seaport_Sand, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_Seaport_Sand_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(p_Seaport_Sand, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_Seaport_Sand_0.001, 3)), adj = 0)
dev.off()
#### NOTE: The difference between Seaport and Sand is significant *
# _____________________