# ---- 0. Load packages ----
library(tidyverse)
library(eulerr) # display Euler diagrams
library(wesanderson)
# __________________________



# ---- 1. Import data ----
## 1-1. Presence/Absence MOTU matrix
### Restrict metabarcoding data to samples
MOTU_table_final <- read.csv('../../data/MOTU_table_final.csv') # load data
pa_MOTU_matrix <- MOTU_table_final %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

### Transform into presence-absence_matrix
pa_MOTU_matrix[pa_MOTU_matrix > 0] <- 1
final_MOTU_matrix <- t(pa_MOTU_matrix) # transpose matrix to get samples in rows and MOTUs in columns


## 1-2. Presence/Absence ASV matrix
### Restrict metabarcoding data to samples
ASV_table_final <- read.csv('../../data/ASV_table_final.csv') # load data
pa_ASV_matrix <- ASV_table_final %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

### Transform into presence-absence_matrix
pa_ASV_matrix[pa_ASV_matrix > 0] <- 1
final_ASV_matrix <- t(pa_ASV_matrix) # transpose matrix to get samples in rows and MOTUs in columns


## 1-3. Site metadata
metadata <- read.csv2('../../data/site_metadata.csv') %>%
  pivot_longer(cols = c(spygen_code_1, spygen_code_2),
               values_to = 'sample') %>% 
  mutate(sample = factor(sample, levels = rownames(final_MOTU_matrix))) %>%
  arrange(sample)


## 1-4. Pool samples
### Transform Presence/Absence matrices
#### Final MOTU matrix
final_MOTU_matrix_new <- t(sapply(split(seq_len(nrow(final_MOTU_matrix)), metadata$pool),
                             function(idx){
                               as.integer(colSums(final_MOTU_matrix[idx, , drop = FALSE]) > 0)
                               }))
colnames(final_MOTU_matrix_new) <- colnames(final_MOTU_matrix)

#### Final ASV matrix
final_ASV_matrix_new <- t(sapply(split(seq_len(nrow(final_ASV_matrix)), metadata$pool),
                                 function(idx){
                                   as.integer(colSums(final_ASV_matrix[idx, , drop = FALSE]) > 0)
                                   }))
colnames(final_ASV_matrix_new) <- colnames(final_ASV_matrix)

### Transform site metadata
metadata_new <- metadata %>% 
  select(-c(sample, name)) %>% 
  unique()
# ________________________



# ---- 2. MOTU Euler diagrams ----
## 2-1. Display Euler diagram by Area
### 2-1-1. Record MOTUs detected in each area
#### West
west_MOTUs <- colnames(final_MOTU_matrix_new)[colSums(final_MOTU_matrix_new[metadata_new$Area == 'West', ]) > 0]
#### Central
central_MOTUs <- colnames(final_MOTU_matrix_new)[colSums(final_MOTU_matrix_new[metadata_new$Area == 'Central', ]) > 0]
#### East
east_MOTUs <- colnames(final_MOTU_matrix_new)[colSums(final_MOTU_matrix_new[metadata_new$Area == 'East', ]) > 0]

### 2-1-2. List the three diagram areas
MOTU_Area_list <- list(West = west_MOTUs, 
                       East = east_MOTUs,
                       Central = central_MOTUs)

### 2-1-3. Generate the diagram
Euler_MOTU_Area <- euler(MOTU_Area_list)
Euler_MOTU_Area_diagram <- plot(Euler_MOTU_Area,
                                fills = wes_palette('AsteroidCity1')[c(1, 4, 3)],
                                alpha = 0.75,
                                quantities = TRUE,
                                labels = list(cex = 1, font = 2, pos = 3),
                                edges = list(lwd = 2))
Euler_MOTU_Area_diagram

### Save in pdf
ggsave(Euler_MOTU_Area_diagram,
       file = '../../outputs/Euler_MOTU_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 2-2. Display Euler diagram by Habitat
### 2-2-1. Record MOTUs detected in each habitat
#### Seaport
seaport_MOTUs <- colnames(final_MOTU_matrix_new)[colSums(final_MOTU_matrix_new[metadata_new$Habitat == 'Seaport', ]) > 0]
#### Rock
rock_MOTUs <- colnames(final_MOTU_matrix_new)[colSums(final_MOTU_matrix_new[metadata_new$Habitat == 'Rock', ]) > 0]
#### Sand
sand_MOTUs <- colnames(final_MOTU_matrix_new)[colSums(final_MOTU_matrix_new[metadata_new$Habitat == 'Sand', ]) > 0]

### 2-2-2. List the three diagram areas
MOTU_Habitat_list <- list(Rock = rock_MOTUs,
                          Sand = sand_MOTUs,
                          Seaport = seaport_MOTUs)

### 2-2-3. Generate the diagram
Euler_MOTU_Habitat <- euler(MOTU_Habitat_list)
Euler_MOTU_Habitat_diagram <- plot(Euler_MOTU_Habitat,
                                   fills = wes_palette('Moonrise2')[c(2, 3, 1)],
                                   alpha = 0.75,
                                   quantities = TRUE,
                                   labels = list(cex = 1, font = 2, pos = 3),
                                   edges = list(lwd = 2))
Euler_MOTU_Habitat_diagram

### Save in pdf
ggsave(Euler_MOTU_Habitat_diagram,
       file = '../../outputs/Euler_MOTU_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')
# ________________________________



# ---- 3. ASV Euler diagrams ----
## 3-1. Display Euler diagram by Area
### 3-1-1. Record ASVs detected in each area
#### West
west_ASVs <- colnames(final_ASV_matrix_new)[colSums(final_ASV_matrix_new[metadata_new$Area == 'West', ]) > 0]
#### Central
central_ASVs <- colnames(final_ASV_matrix_new)[colSums(final_ASV_matrix_new[metadata_new$Area == 'Central', ]) > 0]
#### East
east_ASVs <- colnames(final_ASV_matrix_new)[colSums(final_ASV_matrix_new[metadata_new$Area == 'East', ]) > 0]

### 3-1-2. List the three diagram areas
ASV_Area_list <- list(West = west_ASVs, 
                      East = east_ASVs,
                      Central = central_ASVs)

### 3-1-3. Generate the diagram
Euler_ASV_Area <- euler(ASV_Area_list)
Euler_ASV_Area_diagram <- plot(Euler_ASV_Area,
                               fills = wes_palette('AsteroidCity1')[c(1, 4, 3)],
                               alpha = 0.75,
                               quantities = TRUE,
                               labels = list(cex = 1, font = 2, pos = 3),
                               edges = list(lwd = 2))
Euler_ASV_Area_diagram

### Save in pdf
ggsave(Euler_ASV_Area_diagram,
       file = '../../outputs/Euler_ASV_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 3-2. Display Euler diagram by Habitat
### 3-2-1. Record ASVs detected in each habitat
#### Seaport
seaport_ASVs <- colnames(final_ASV_matrix_new)[colSums(final_ASV_matrix_new[metadata_new$Habitat == 'Seaport', ]) > 0]
#### Rock
rock_ASVs <- colnames(final_ASV_matrix_new)[colSums(final_ASV_matrix_new[metadata_new$Habitat == 'Rock', ]) > 0]
#### Sand
sand_ASVs <- colnames(final_ASV_matrix_new)[colSums(final_ASV_matrix_new[metadata_new$Habitat == 'Sand', ]) > 0]

### 3-2-2. List the three diagram areas
ASV_Habitat_list <- list(Rock = rock_ASVs,
                         Sand = sand_ASVs,
                         Seaport = seaport_ASVs)

### 3-2-3. Generate the diagram
Euler_ASV_Habitat <- euler(ASV_Habitat_list)
Euler_ASV_Habitat_diagram <- plot(Euler_ASV_Habitat,
                                  fills = wes_palette('Moonrise2')[c(2, 3, 1)],
                                  alpha = 0.75,
                                  quantities = TRUE,
                                  labels = list(cex = 1, font = 2, pos = 3),
                                  edges = list(lwd = 2))
Euler_ASV_Habitat_diagram

### Save in pdf
ggsave(Euler_ASV_Habitat_diagram,
       file = '../../outputs/Euler_ASV_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')
# ________________________________



# ---- 4. ASVs' subset MOTUs Euler diagrams ----
## 4-1. Subset MOTUs
### 4-1-1. Select polymorphic MOTUs
polymorphic_MOTUs <- ASV_table_final %>%
  count(MOTU, name = 'n_ASV') %>% 
  filter(n_ASV >= 2) %>% # only keep MOTUs with at least 2 ASVs
  pull(MOTU)

### 4-1-2. Select MOTUs shared between all areas and all habitats
Area_shared_MOTUs <- intersect(central_MOTUs, intersect(west_MOTUs, east_MOTUs))
Habitat_shared_MOTUs <- intersect(seaport_MOTUs, intersect(rock_MOTUs, sand_MOTUs))
shared_MOTUs <- intersect(Area_shared_MOTUs, Habitat_shared_MOTUs)

### 4-1-3. Restrict to the intersection of both subsets
MOTUs_to_keep <- intersect(poly, shared_MOTUs)

### 4-1-4. Record the subset MOTU table
MOTU_table_subset <- MOTU_table_final %>% 
  filter(id %in% MOTUs_to_keep)
# write.csv(MOTU_table_subset, '../data/MOTU_table_subset.csv', row.names = F)


## 4-2. Examine taxonomic composition of subset MOTUs
MOTU_table_subset_kingdoms <- MOTU_table_subset %>%
  count(kingdom) %>%
  mutate(percent = n / sum(n) * 100,
         kingdom = ifelse(is.na(kingdom), 'Unassigned', kingdom),
         kingdom = factor(kingdom, levels = c('Metazoa', 'Chromista', 'Plantae', 'Fungi', 'Unassigned')))
  
MOTU_table_subset_composition_plot <- ggplot(MOTU_table_subset_kingdoms,
                                             aes(x = '', y = percent, fill = kingdom)) +
  geom_col(width = 1, color = 'black') +
  scale_fill_manual(values = c('#b3cde3ff', '#ccebc5ff', '#fbb4aeff', 'black', 'grey70')) +
  coord_polar(theta = 'y') +
  geom_text(data = subset(MOTU_table_subset_kingdoms, kingdom == 'Fungi'),
            aes(label = 'Fungi')) +
  geom_text(data = subset(MOTU_table_subset_kingdoms, kingdom == 'Unassigned'),
            aes(label = 'Unassigned\nEukaryota')) +
  theme_void() +
  theme(legend.position = 'none')
MOTU_table_subset_composition_plot

### Save in pdf
ggsave(MOTU_table_subset_composition_plot,
       file = '../../outputs/MOTU_table_subset_composition_plot.pdf',
       height = 5,
       width = 5,
       unit = 'cm')


## 4-3. Presence/Absence ASVs' subset MOTUs matrix
### Restrict metabarcoding data to samples
ASV_table_final <- read.csv('../../data/ASV_table_final.csv') # load data
pa_ASV_subset_matrix <- ASV_table_final %>%
  filter(MOTU %in% MOTUs_to_keep) %>% 
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

### Transform into presence-absence_matrix
pa_ASV_subset_matrix[pa_ASV_subset_matrix > 0] <- 1
final_ASV_subset_matrix <- t(pa_ASV_subset_matrix) # transpose matrix to get samples in rows and MOTUs in columns

#### Pool samples
final_ASV_subset_matrix_new <- t(sapply(split(seq_len(nrow(final_ASV_subset_matrix)), metadata$pool),
                                        function(idx){
                                          as.integer(colSums(final_ASV_subset_matrix[idx, , drop = FALSE]) > 0)
                                          }))
colnames(final_ASV_subset_matrix_new) <- colnames(final_ASV_subset_matrix)


## 4-4. Display Euler diagram by Area
### 4-4-1. Record ASVs detected in each area
#### West
west_ASVs_subset <- colnames(final_ASV_subset_matrix_new)[colSums(final_ASV_subset_matrix_new[metadata_new$Area == 'West', ]) > 0]
#### Central
central_ASVs_subset <- colnames(final_ASV_subset_matrix_new)[colSums(final_ASV_subset_matrix_new[metadata_new$Area == 'Central', ]) > 0]
#### East
east_ASVs_subset <- colnames(final_ASV_subset_matrix_new)[colSums(final_ASV_subset_matrix_new[metadata_new$Area == 'East', ]) > 0]

### 4-4-2. List the three diagram areas
ASV_subset_Area_list <- list(West = west_ASVs_subset,
                             East = east_ASVs_subset,
                             Central = central_ASVs_subset)

### 4-4-3. Generate the diagram
Euler_ASV_subset_Area <- euler(ASV_subset_Area_list)
Euler_ASV_subset_Area_diagram <- plot(Euler_ASV_subset_Area,
                                      fills = wes_palette('AsteroidCity1')[c(1, 4, 3)],
                                      alpha = 0.75,
                                      quantities = TRUE,
                                      labels = list(cex = 1, font = 2, pos = 3),
                                      edges = list(lwd = 2))
Euler_ASV_subset_Area_diagram

### Save in pdf
ggsave(Euler_ASV_subset_Area_diagram,
       file = '../../outputs/Euler_ASV_subset_Area_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')


## 4-5. Display Euler diagram by Habitat
### 4-5-1. Record ASVs detected in each habitat
#### Seaport
seaport_ASVs_subset <- colnames(final_ASV_subset_matrix_new)[colSums(final_ASV_subset_matrix_new[metadata_new$Habitat == 'Seaport', ]) > 0]
#### Rock
rock_ASVs_subset <- colnames(final_ASV_subset_matrix_new)[colSums(final_ASV_subset_matrix_new[metadata_new$Habitat == 'Rock', ]) > 0]
#### Sand
sand_ASVs_subset <- colnames(final_ASV_subset_matrix_new)[colSums(final_ASV_subset_matrix_new[metadata_new$Habitat == 'Sand', ]) > 0]

### 4-5-2. List the three diagram areas
ASV_subset_Habitat_list <- list(Rock = rock_ASVs_subset,
                                Sand = sand_ASVs_subset,
                                Seaport = seaport_ASVs_subset)

### 4-5-3. Generate the diagram
Euler_ASV_subset_Habitat <- euler(ASV_subset_Habitat_list)
Euler_ASV_subset_Habitat_diagram <- plot(Euler_ASV_subset_Habitat,
                                         fills = wes_palette('Moonrise2')[c(2, 3, 1)],
                                         alpha = 0.75,
                                         quantities = TRUE,
                                         labels = list(cex = 1, font = 2, pos = 3),
                                         edges = list(lwd = 2))
Euler_ASV_subset_Habitat_diagram

### Save in pdf
ggsave(Euler_ASV_subset_Habitat_diagram,
       file = '../../outputs/Euler_ASV_subset_Habitat_plot.pdf',
       height = 8,
       width = 8,
       unit = 'cm')
# ________________________________