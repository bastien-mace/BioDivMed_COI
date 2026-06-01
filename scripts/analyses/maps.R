# ---- 0. Load packages ----
library(tidyverse)
library(rnaturalearth)
library(sf)
library(scatterpie)
library(ggspatial)
library(patchwork)
library(maptiles)
library(terra)
library(tidyterra)
library(wesanderson)
library(grid)
# __________________________



# ---- 1. Mediterranean map ----
## 1-1. Extract the study area (Occitanie)
world <- ne_countries(scale = 'large')
box <- st_bbox(c(xmin = 2.9, xmax = 4.5,
                 ymin = 42.5, ymax = 43.7),
               crs = st_crs(4326))
box_sf <- st_as_sfc(box)


## 1-2. Extract the areas' coordinates
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


## 1-3. Create the map
main_map <- ggplot() +
  geom_sf(data = world, fill = '#666666ff', color = 'black') +
  geom_scatterpie(data = st_drop_geometry(sampled_sites_sf),
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
  theme_classic() +
  theme(axis.title = element_text(size = 15),
        legend.position = 'none')
main_map


## 1-4. Create an inset representing France
france <- world %>%
  filter(admin == 'France')
inset_map <- ggplot() +
  geom_sf(data = france, fill = 'grey80', color = 'black') +
  geom_sf(data = box_sf, fill = NA, color = 'black', linewidth = 1) +
  coord_sf(xlim = c(-5, 10), ylim = c(41, 52), # restrict to European France
           expand = FALSE) +
  theme_void()
inset_map


## 1-5. Plot the final map
final_map <- main_map +
  inset_element(inset_map,
                left = 0, right = 0.35,
                bottom = 0.7, top = 1)
final_map

### Save in pdf
ggsave(final_map,
       file = '../../outputs/Mediterranean_map.pdf',
       height = 15,
       width = 15,
       unit = 'cm')
# ______________________________



# ---- 2. Local maps ----
## 2-1. Function to create OpenStreetMap views inside each area
make_local_map <- function(area){
  
  ### 2-1-1. Import data
  sampled_sites <- read.csv2('../../data/site_metadata.csv') %>% 
    filter(Area == area) %>% 
    select(Habitat, Longitude, Latitude) %>%
    mutate(Habitat = factor(Habitat,
                            levels = c('Seaport', 'Rock', 'Sand')))
  sampled_sites_sf <- st_as_sf(sampled_sites,
                               coords = c('Longitude', 'Latitude'),
                               crs = 4326)
  points <- st_transform(sampled_sites_sf, 3857) # transform in metric projection
  
  ### 2-1-2. Prepare a circle around the data points
  center <- st_centroid(st_union(points)) # calculate the area centroid
  dists <- st_distance(points, center) # calculate the distance to the center for each point
  radius <- max(dists) # retain the maximal distance
  radius <- as.numeric(radius) * 1.2 # create the circle radius
  circle <- st_buffer(center, dist = radius)
  
  ### 2-1-3. Create OpenStreetMap circle
  bbox <- st_bbox(circle) %>%
    st_as_sfc() %>%
    st_transform(4326)
  tiles <- get_tiles(bbox,
                     provider = 'OpenStreetMap',
                     crop = TRUE,
                     zoom = 16)
  tiles_proj <- project(tiles, 'EPSG:3857')
  tiles_crop <- crop(tiles_proj, vect(circle))
  tiles_mask <- mask(tiles_crop, vect(circle))
  points_circle <- st_intersection(points, circle)
  
  ### 2-1-3. Plot the area
  ggplot() +
    geom_spatraster_rgb(data = tiles_mask) +
    geom_sf(data = points_circle,
            aes(fill = Habitat),
            shape = 21, size = 4, color = 'black') +
    geom_sf(data = circle,
            fill = NA, color = 'black', linewidth = 1) +
    scale_fill_manual(name = 'Habitat',
                      values = wes_palette('Moonrise2')[1:3]) +
    annotation_scale(pad_x = unit(1, 'cm'),
                     pad_y = unit(1, 'cm'),
                     width_hint = 0.3) +
    annotation_north_arrow(pad_x = unit(1, 'cm'),
                           pad_y = unit(2, 'cm'),
                           which_north = 'true',
                           height = unit(1, 'cm'),
                           width  = unit(0.75, 'cm')) +
    coord_sf(crs = 3857) +
    guides(fill = guide_legend(title.position = 'top',
                               title.hjust = 0.5, nrow = 1)) +
    theme_void() +
    theme(legend.background = element_rect(color = 'black', fill = 'white'),
          legend.position = 'top',
          legend.title.align = 0.5,
          legend.title = element_text(face = 'bold'),
          legend.margin = margin(5, 10, 5, 10))
}


## 2-2. Plot each area
### 2-2-1. St-Cyprien (West)
west_map <- make_local_map(area = 'West')
west_map

#### Save in pdf
ggsave(west_map,
       file = '../../outputs/west_map.pdf',
       height = 12.5,
       width = 12.5,
       unit = 'cm')

### 2-2-2. Agde (Central)
central_map <- make_local_map(area = 'Central')
central_map

#### Save in pdf
ggsave(central_map,
       file = '../../outputs/central_map.pdf',
       height = 12.5,
       width = 12.5,
       unit = 'cm')

### 2-2-3. Le Grau-du-Roi
east_map <- make_local_map(area = 'East')
east_map

#### Save in pdf
ggsave(east_map,
       file = '../../outputs/east_map.pdf',
       height = 12.5,
       width = 12.5,
       unit = 'cm')
# _______________________