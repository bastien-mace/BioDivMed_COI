# ---- 0. Load packages ----
library(tidyverse)
library(vegan) # beta-diversity analyses
library(wesanderson)
# __________________________



# ---- 1. Import data ----
## 1-1. Site metadata
metadata <- read.csv2('../../data/site_metadata.csv') %>%
  pivot_longer(cols = c(spygen_code_1, spygen_code_2),
               values_to = 'sample') %>% 
  arrange(sample)

### Pool the problematic sample
samples_to_pool <- c('SPY232133', 'SPY232146')
pooled_name <- 'SPY232133_SPY232146'

### Modify metadata
metadata_new <- metadata[!metadata$sample %in% samples_to_pool, ]
pooled_meta <- metadata[metadata$sample == samples_to_pool[1], ]
pooled_meta$sample <- pooled_name
metadata_new <- rbind(metadata_new, pooled_meta)


## 1-2. Load distance matrices
### Jaccard matrix (β-SD)
load('../../data/beta_SD.Rdata') # built from INTER_MOTU_BETA.R
beta_SD <- jac_dist$beta.jtu

### Jost's D matrix (β-GD)
load('../../data/beta_GD.Rdata') # built from INTRA_MOTU_BETA_Jost_D.R
beta_GD <- D$mean_D_matrix
# ________________________



# ---- 2. β-SGDC ----
## 2-1. Bootstrap Mantel tests to account for site replication
set.seed(2026)
nboot <- 999

### Create an empty list
mantels <- data.frame(r = as.numeric(),
                      p = as.numeric())

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
  
  #### Subset the distance matrices to these samples
  beta_SD_sub<- as.dist(as.matrix(beta_SD)[samples, samples, drop = FALSE])
  beta_GD_sub <- as.dist(as.matrix(beta_GD)[samples, samples, drop = FALSE])
  
  #### Compute Mantel tests
  mt <- mantel(beta_SD_sub, beta_GD_sub)
  mt_df <- data.frame(r = mt$statistic,
                      p = mt$signif)
  
  #### Record the results
  mantels <- rbind(mantels,
                   mt_df)
}

## 2-2. Check Mantel r and p-value distributions
#### 2-2-1. Mantel r
##### Extract values
mean_r <- mean(mantels$r)

#### Save in pdf
pdf('../../outputs/BETA_SGDC_r.pdf',
    height = 4,
    width = 6)
hist(mantels$r, col = 'skyblue',
     main = 'Mantel r bootstrap distribution', xlab = 'Bootstrapped Mantel r')
text(x = max(hist(mantels$r, plot = FALSE)$breaks) * 0.9, y = max(hist(mantels$r, plot = FALSE)$counts),
     labels = paste0('Mean Mantel r: ', round(mean_r, 3)), adj = 0)
dev.off()

### 2-2-2. p-value
#### Extract values
mean_p <- mean(mantels$p)
prop_0.05 <- mean(mantels$p <= 0.05)
prop_0.01 <- mean(mantels$p <= 0.01)
prop_0.001 <- mean(mantels$p <= 0.001)

#### Save in pdf
pdf('../../outputs/BETA_SGDC_p.pdf',
    height = 4,
    width = 6)
hist(mantels$p, col = 'skyblue', xlim = c(0, 0.1),
     main = 'Habitat p-value bootstrap distribution', xlab = 'Bootstrapped p-value')
abline(v = 0.05, col = 'red', lty = 2)
text(x = 0.055, y = max(hist(mantels$p, plot = FALSE)$counts) * 1,
     labels = paste0('Mean p-value: ', round(mean_p, 3)), adj = 0)
text(x = 0.055, y = max(hist(mantels$p, plot = FALSE)$counts) * 0.9,
     labels = paste0('p <= 0.05: ', round(prop_0.05, 3)), adj = 0)
text(x = 0.055, y = max(hist(mantels$p, plot = FALSE)$counts) * 0.8,
     labels = paste0('p <= 0.01: ', round(prop_0.01, 3)), adj = 0)
text(x = 0.055, y = max(hist(mantels$p, plot = FALSE)$counts) * 0.7,
     labels = paste0('p <= 0.001: ', round(prop_0.001, 3)), adj = 0)
dev.off()

#### NOTE: The correlation is strong (Mantel r = 0.62) and significant ***
# ___________________________________



# ---- 3. Plot the correlation ----
## 3-1. Remove comparisons of samples from the same pool
same_pool <- outer(metadata_new$pool, metadata_new$pool, FUN = '==')
beta_SD_mat <- as.matrix(beta_SD)
beta_GD_mat <- as.matrix(beta_GD)
beta_SD_mat[same_pool] <- NA
beta_GD_mat[same_pool] <- NA


## 3-2. Precise the habitats of the sample compared
### Create a comparison matrix
hab1 <- outer(metadata_new$Habitat, metadata_new$Habitat, FUN = function(x, y) x)
hab2 <- outer(metadata_new$Habitat, metadata_new$Habitat, FUN = function(x, y) y)
hab_comparison <- matrix(NA, nrow = nrow(hab1), ncol = ncol(hab1))

### Record Seaport vs Seaport comparisons
hab_comparison[hab1 == 'Seaport' & hab2 == 'Seaport'] <- 'Seaport vs Seaport'

### Record Seaport vs Outside comparisons
hab_comparison[(hab1 == 'Seaport' & hab2 %in% c('Sand', 'Rock')) |
                 (hab2 == 'Seaport' & hab1 %in% c('Sand', 'Rock'))] <- 'Seaport vs Outside'

### Record Outside vs Outside comparisons
hab_comparison[hab1 %in% c('Sand', 'Rock') & hab2 %in% c('Sand', 'Rock')] <- 'Outside vs Outside'


## 3-3. Create a dataframe recording the values
beta_values <- data.frame(SD = as.numeric(beta_SD_mat[lower.tri(beta_SD_mat)]),
                          GD = as.numeric(beta_GD_mat[lower.tri(beta_GD_mat)]),
                          Comparison = hab_comparison[lower.tri(hab_comparison)]) %>%
  drop_na() %>%
  mutate(Comparison = factor(Comparison, levels = c('Seaport vs Seaport', 'Seaport vs Outside', 'Outside vs Outside')))


## 3-3. Plot the filtered data
BETA_SGDC_plot <- ggplot(beta_values,
                         aes(x = SD, y = GD,
                             fill = Comparison)) +
  geom_point(size = 2, shape = 21, color = 'black') +
  scale_fill_manual(values = wes_palette('IsleofDogs1')) +
  guides(fill = guide_legend(title.position = 'top',
                             title.hjust = 0.5, nrow = 1)) +
  xlab('Jaccard dissimilarity (β-SD)') +
  ylab("Jost's D differentiation (β-GD)") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.y = element_text(colour = 'black', size = 15),
        axis.title.x = element_text(colour = 'black', size = 15),
        legend.background = element_rect(color = 'black', fill = 'white'),
        legend.position = 'bottom',
        legend.title = element_text(face = 'bold'))
BETA_SGDC_plot

### Save in pdf
ggsave(BETA_SGDC_plot,
       file = '../../outputs/BETA_SGDC_plot.pdf',
       height = 10,
       width = 10,
       unit = 'cm')
#  ________________________________