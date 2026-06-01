# ---- 0. Load packages ----
library(tidyverse)
library(rnaturalearth)
library(sf)
library(scatterpie)
library(ggspatial)
library(wesanderson)
library(grid)
# __________________________



# ---- 1. Prepare data ----
## 1-1. Mediterranean frame
### 1-1-1. Extract the study area (Occitanie)
world <- ne_countries(scale = 'large')
box <- st_bbox(c(xmin = 2.9, xmax = 4.5,
                 ymin = 42.5, ymax = 43.7),
               crs = st_crs(4326))
box_sf <- st_as_sfc(box)

### 1-1-2. Extract the areas' coordinates
sampled_sites <- read.csv2('../../data/site_metadata.csv') %>% 
  filter(Habitat == 'Seaport') %>% 
  select(Name, Longitude, Latitude) %>% 
  mutate(Seaport = as.numeric(1),
         Rock = as.numeric(1),
         Sand = as.numeric(1))
sampled_sites_sf <- st_as_sf(sampled_sites,
                             coords = c('Longitude', 'Latitude'),
                             crs = 4326)
coords <- st_coordinates(sampled_sites_sf)
sampled_sites_sf$lon <- coords[,1]
sampled_sites_sf$lat <- coords[,2]

### 1-1-3. Prepare data for scatterpies
pie_west <- filter(sampled_sites_sf, Name == 'Saint-Cyprien')
pie_central <- filter(sampled_sites_sf, Name == 'Agde')
pie_east <- filter(sampled_sites_sf, Name == 'Le Grau-du-Roi')


## 1-2. NIS presence/absence data
### 1-2-1. Presence/Absence matrix
#### Restrict MOTU table to samples' columns
MOTU_table_NIS <- read.csv('../../data/MOTU_table_final.csv') %>%
  filter(NIS == 'NIS' | ltg_name == 'Schizoporella pungens')
pa_matrix <- MOTU_table_NIS %>%
  select(id, starts_with('SPY')) %>% 
  column_to_rownames('id')

#### Transform into presence-absence_matrix
pa_matrix[pa_matrix > 0] <- 1
colSums(pa_matrix)
final_matrix <- t(pa_matrix) # transpose matrix to get samples in rows and MOTUs in columns

### 1-2-2. Site metadata
metadata <- read.csv2('../../data/site_metadata.csv') %>%
  pivot_longer(cols = c(spygen_code_1, spygen_code_2),
               values_to = 'sample') %>% 
  mutate(sample = factor(sample, levels = rownames(final_matrix))) %>%
  arrange(sample)

### 1-2-3. Pool the problematic sample
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

### 1-2-4 Merge MOTU table to metadata
pa_long <- final_matrix_new %>%
  as.data.frame() %>%
  rownames_to_column('sample') %>%
  pivot_longer(cols = -sample,
               names_to = 'id',
               values_to = 'PA') %>%
  left_join(metadata_new) %>% 
  left_join(MOTU_table_NIS %>% 
              select(id, ltg_name))
# _________________________



# ---- 2. Map all NIS MOTUs ----
unique(MOTU_table_NIS$ltg_name)

## 2-1. Amphibalanus amphitrite
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Amphibalanus amphitrite'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 1
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Amphibalanus amphitrite'), c(Area, Habitat)))
pie_central$Seaport <- 1
pie_central$Rock <- 1
pie_central$Sand <- 0

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Amphibalanus amphitrite'), c(Area, Habitat)))
pie_east$Seaport <- 1
pie_east$Rock <- 1
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Amphibalanus_amphitrite_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Amphibalanus amphitrite') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Amphibalanus_amphitrite_map

### Save in pdf
ggsave(Amphibalanus_amphitrite_map,
       file = '../../outputs/Amphibalanus_amphitrite_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-2. Coryne eximia
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Coryne eximia'), c(Area, Habitat)))
pie_west$Seaport <- 0
pie_west$Rock <- 0
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Coryne eximia'), c(Area, Habitat)))
pie_central$Seaport <- 0
pie_central$Rock <- 0
pie_central$Sand <- 1

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Coryne eximia'), c(Area, Habitat)))
pie_east$Seaport <- 0
pie_east$Rock <- 0
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Coryne_eximia_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[3]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Coryne eximia') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Coryne_eximia_map

### Save in pdf
ggsave(Coryne_eximia_map,
       file = '../../outputs/Coryne_eximia_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-3. Ruditapes philippinarum
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Ruditapes philippinarum'), c(Area, Habitat)))
pie_west$Seaport <- 0
pie_west$Rock <- 0
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Ruditapes philippinarum'), c(Area, Habitat)))
pie_central$Seaport <- 0
pie_central$Rock <- 1
pie_central$Sand <- 1

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Ruditapes philippinarum'), c(Area, Habitat)))
pie_east$Seaport <- 0
pie_east$Rock <- 1
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Ruditapes_philippinarum_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[2:3]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Ruditapes philippinarum') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Ruditapes_philippinarum_map

### Save in pdf
ggsave(Ruditapes_philippinarum_map,
       file = '../../outputs/Ruditapes_philippinarum_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-4. Caprella scaura
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Caprella scaura'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 0
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Caprella scaura'), c(Area, Habitat)))
pie_central$Seaport <- 0
pie_central$Rock <- 0
pie_central$Sand <- 0

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Caprella scaura'), c(Area, Habitat)))
pie_east$Seaport <- 0
pie_east$Rock <- 0
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Caprella_scaura_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Caprella scaura') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Caprella_scaura_map

### Save in pdf
ggsave(Caprella_scaura_map,
       file = '../../outputs/Caprella_scaura_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-5. Hydroides elegans
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Hydroides elegans'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 0
pie_west$Sand <- 1

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Hydroides elegans'), c(Area, Habitat)))
pie_central$Seaport <- 1
pie_central$Rock <- 0
pie_central$Sand <- 0

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Hydroides elegans'), c(Area, Habitat)))
pie_east$Seaport <- 1
pie_east$Rock <- 1
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Hydroides_elegans_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Hydroides elegans') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Hydroides_elegans_map

### Save in pdf
ggsave(Hydroides_elegans_map,
       file = '../../outputs/Hydroides_elegans_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-6. Oithona davisae
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Oithona davisae'), c(Area, Habitat)))
pie_west$Seaport <- 0
pie_west$Rock <- 1
pie_west$Sand <- 1

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Oithona davisae'), c(Area, Habitat)))
pie_central$Seaport <- 1
pie_central$Rock <- 1
pie_central$Sand <- 1

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Oithona davisae'), c(Area, Habitat)))
pie_east$Seaport <- 1
pie_east$Rock <- 1
pie_east$Sand <- 1

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Oithona_davisae_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Oithona davisae') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Oithona_davisae_map

### Save in pdf
ggsave(Oithona_davisae_map,
       file = '../../outputs/Oithona_davisae_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-7. Pseudodiaptomus marinus
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Pseudodiaptomus marinus'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 1
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Pseudodiaptomus marinus'), c(Area, Habitat)))
pie_central$Seaport <- 1
pie_central$Rock <- 1
pie_central$Sand <- 1

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Pseudodiaptomus marinus'), c(Area, Habitat)))
pie_east$Seaport <- 0
pie_east$Rock <- 1
pie_east$Sand <- 1

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Pseudodiaptomus_marinus_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Pseudodiaptomus marinus') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Pseudodiaptomus_marinus_map

### Save in pdf
ggsave(Pseudodiaptomus_marinus_map,
       file = '../../outputs/Pseudodiaptomus_marinus_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-8. Polycerella emertoni
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Polycerella emertoni'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 1
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Polycerella emertoni'), c(Area, Habitat)))
pie_central$Seaport <- 0
pie_central$Rock <- 1
pie_central$Sand <- 0

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Polycerella emertoni'), c(Area, Habitat)))
pie_east$Seaport <- 1
pie_east$Rock <- 0
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Polycerella_emertoni_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:2]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Polycerella emertoni') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Polycerella_emertoni_map

### Save in pdf
ggsave(Polycerella_emertoni_map,
       file = '../../outputs/Polycerella_emertoni_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-9. Paracartia grani
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Paracartia grani'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 1
pie_west$Sand <- 1

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Paracartia grani'), c(Area, Habitat)))
pie_central$Seaport <- 1
pie_central$Rock <- 1
pie_central$Sand <- 1

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Paracartia grani'), c(Area, Habitat)))
pie_east$Seaport <- 1
pie_east$Rock <- 1
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Paracartia_grani_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Paracartia grani') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Paracartia_grani_map

### Save in pdf
ggsave(Paracartia_grani_map,
       file = '../../outputs/Paracartia_grani_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-10. Amphibalanus eburneus
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Amphibalanus eburneus'), c(Area, Habitat)))
pie_west$Seaport <- 0
pie_west$Rock <- 0
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Amphibalanus eburneus'), c(Area, Habitat)))
pie_central$Seaport <- 0
pie_central$Rock <- 1
pie_central$Sand <- 0

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Amphibalanus eburneus'), c(Area, Habitat)))
pie_east$Seaport <- 0
pie_east$Rock <- 1
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Amphibalanus_eburneus_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[2]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Amphibalanus eburneus') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Amphibalanus_eburneus_map

### Save in pdf
ggsave(Amphibalanus_eburneus_map,
       file = '../../outputs/Amphibalanus_eburneus_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-11. Amathia verticillata
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Amathia verticillata'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 0
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Amathia verticillata'), c(Area, Habitat)))
pie_central$Seaport <- 1
pie_central$Rock <- 0
pie_central$Sand <- 0

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Amathia verticillata'), c(Area, Habitat)))
pie_east$Seaport <- 1
pie_east$Rock <- 0
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Amathia_verticillata_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Amathia verticillata') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Amathia_verticillata_map

### Save in pdf
ggsave(Amathia_verticillata_map,
       file = '../../outputs/Amathia_verticillata_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-12. Watersipora subtorquata
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Watersipora subtorquata'), c(Area, Habitat)))
pie_west$Seaport <- 0
pie_west$Rock <- 0
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Watersipora subtorquata'), c(Area, Habitat)))
pie_central$Seaport <- 0
pie_central$Rock <- 0
pie_central$Sand <- 0

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Watersipora subtorquata'), c(Area, Habitat)))
pie_east$Seaport <- 0
pie_east$Rock <- 1
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Watersipora_subtorquata_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[2]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Watersipora subtorquata/subatra') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Watersipora_subtorquata_map

### Save in pdf
ggsave(Watersipora_subtorquata_map,
       file = '../../outputs/Watersipora_subtorquata_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-13. Pseudopolydora paucibranchiata
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Pseudopolydora paucibranchiata'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 0
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Pseudopolydora paucibranchiata'), c(Area, Habitat)))
pie_central$Seaport <- 0
pie_central$Rock <- 1
pie_central$Sand <- 1

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Pseudopolydora paucibranchiata'), c(Area, Habitat)))
pie_east$Seaport <- 1
pie_east$Rock <- 1
pie_east$Sand <- 1

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Pseudopolydora_paucibranchiata_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:3]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Pseudopolydora paucibranchiata') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Pseudopolydora_paucibranchiata_map

### Save in pdf
ggsave(Pseudopolydora_paucibranchiata_map,
       file = '../../outputs/Pseudopolydora_paucibranchiata_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')


## 2-14. Schizoporella pungens
### West Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'West'), ltg_name == 'Schizoporella pungens'), c(Area, Habitat)))
pie_west$Seaport <- 1
pie_west$Rock <- 1
pie_west$Sand <- 0

### Central Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'Central'), ltg_name == 'Schizoporella pungens'), c(Area, Habitat)))
pie_central$Seaport <- 1
pie_central$Rock <- 0
pie_central$Sand <- 0

### East Area
unique(select(filter(filter(pa_long, PA == 1 & Area == 'East'), ltg_name == 'Schizoporella pungens'), c(Area, Habitat)))
pie_east$Seaport <- 1
pie_east$Rock <- 0
pie_east$Sand <- 0

### Merge pies
pie_data <- bind_rows(pie_west, pie_central, pie_east) %>%
  st_drop_geometry()

### PLot the map
Schizoporella_pungens_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_point(data = sampled_sites_sf,
             aes(x = lon, y = lat),
             color = 'black', size = 3) +
  geom_scatterpie(data = pie_data,
                  aes(x = lon, y = lat),
                  cols = c('Seaport', 'Rock', 'Sand'),
                  color = 'black',
                  pie_scale = 2) +
  geom_sf_text(data = sampled_sites_sf,
               aes(label = Name),
               size = 4,
               nudge_y = 0.1) +
  coord_sf(xlim = c(2.8, 4.55),
           ylim = c(42.5, 43.7),
           expand = FALSE) +
  scale_fill_manual(values = wes_palette('Moonrise2')[1:2]) +
  scale_x_continuous(breaks = c(3, 3.5, 4, 4.5)) +
  scale_y_continuous(breaks = c(42.5, 43, 43.5)) +
  xlab('Longitude') + ylab('Latitude') +
  ggtitle('Schizoporella pungens') +
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        plot.title = element_text(hjust = 0.5, face = 'italic', size = 15),
        legend.position = 'none')
Schizoporella_pungens_map

### Save in pdf
ggsave(Schizoporella_pungens_map,
       file = '../../outputs/Schizoporella_pungens_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')
# ______________________________