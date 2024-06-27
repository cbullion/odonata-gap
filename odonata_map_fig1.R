

#create map of inaturalist records of odonated in northeastern united states until end of 2019

library(sf)
library(sp)
library(ggplot2)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(ggspatial)

#Load the CSV file
data <- read.csv("observations-453205-1.csv")


#  Filter the dataset for relevant columns and clean the data
# The relevant columns are 'latitude', 'longitude', 'species_guess', and 'observed_on'
data <- data[!is.na(data$latitude) & !is.na(data$longitude), c("latitude", "longitude", "species_guess", "observed_on")]

# Convert observed_on to Date type
data$observed_on <- as.Date(data$observed_on)

#Create a spatial grid based on the extent of the data
xmin <- min(data$longitude)
xmax <- max(data$longitude)
ymin <- min(data$latitude)
ymax <- max(data$latitude)

bbox <- st_bbox(c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), crs = st_crs(4326))

# Create an empty grid with 0.15 degree cells
cellsize <- 0.15
grid <- st_make_grid(st_as_sfc(bbox), cellsize = cellsize)

# Convert grid to sf object
grid_sf <- st_sf(geometry = grid)

#Convert data to sf object
data_sf <- st_as_sf(data, coords = c("longitude", "latitude"), crs = 4326)

#Perform a spatial join to count the number of observations in each grid cell
intersections <- st_intersects(data_sf, grid_sf)

# Flatten the list of intersections and count occurrences of each grid cell
grid_ids <- unlist(intersections)
grid_counts <- table(grid_ids)

# Convert counts to data frame
counts_df <- as.data.frame(grid_counts)
names(counts_df) <- c("grid_id", "count")

# Add counts to the grid_sf object
grid_sf$count <- 0
grid_sf$count[as.integer(as.character(counts_df$grid_id))] <- counts_df$count

#Load political map, state boundaries, and lakes (including Lake Erie)
world <- ne_countries(scale = "medium", returnclass = "sf")
states <- ne_states(country = "United States of America", returnclass = "sf")
lakes <- ne_download(scale = "medium", type = "lakes", category = "physical", returnclass = "sf")

# Fix any invalid geometries in the lakes data
lakes <- st_make_valid(lakes)

#Crop the political map to the extent of the heatmap data
states <- st_crop(states, bbox)
world <- st_crop(world, bbox)
lakes <- st_crop(lakes, bbox)

# Define city and sampling site locations
cities <- data.frame(
  name = c("Cleveland OH", "", "Washington DC", "Norton VA"),#cut out new york, it's out of frame
  longitude = c(-81.6944, -74.0060, -77.0369, -82.6294),
  latitude = c(41.4993, 40.7128, 38.9072, 36.9337)
)

# Convert cities to sf object
cities_sf <- st_as_sf(cities, coords = c("longitude", "latitude"), crs = 4326)

sites <- data.frame(
  name = c("A", "B", "C", "D", "E"),
  longitude = c(-81.6, -81.5, -82.4, -83.0, -82.6),
  latitude = c(41.4, 40.1, 38.2, 37.3, 37.0)
)

# Convert sites to sf object
sites_sf <- st_as_sf(sites, coords = c("longitude", "latitude"), crs = 4326)



# Define the custom color gradient
custom_colors <- c("transparent", "#fee08b", "#E4A465", "#E49365", "#E48265", "#E46565")
custom_values <- c(0, 0.01, 0.1, 0.3, 0.6, 1)  # Exponential-like distribution

#Plot the heatmap with background map, state boundaries, city labels, and Lake Erie
ggplot() +
  geom_sf(data = world, fill = "cornsilk", color = "gray80") +
  geom_sf(data = lakes, fill = "lightblue", color = NA) +
  geom_sf(data = grid_sf, aes(fill = squish(count, range = c(0, 150))), color = NA) +
  geom_sf(data = states, fill = NA, color = "black") + # Add state boundaries on top
  geom_sf(data = cities_sf, color = "black", fill="lightskyblue3",size = 5, pch=21) +
  geom_label(data = cities, aes(x = longitude, y = latitude, label = name), 
             color = "black", size = 4, fill = "cornsilk", label.size = NA, label.padding = unit(0.15, "lines"), 
             nudge_y = -0.3, nudge_x = -1.3) + # Add white background to labels
  geom_sf(data = sites_sf, color = "black", size = 4, pch=18) +
  geom_text(data = sites, aes(x = longitude, y = latitude, label = name), 
            color = "black", size = 5, nudge_y = 0.2, nudge_x = 0.3, fontface="bold") +
  scale_fill_gradientn(
    name = "Observations",
    colors = custom_colors,
    values = custom_values,
    limits = c(0, 150), # Scale from 0 to 150 observations
    oob = squish,       # Apply squish to out-of-bounds values
    na.value = "transparent"  # Handle NA values
  ) +
  coord_sf(xlim = c(xmin+0.1, xmax-0.1), ylim = c(ymin+0.1, ymax-0.1), expand = FALSE) + # Zoom to data extent
  labs(
    title = "",
    x = "Longitude",
    y = "Latitude"
  ) +
  annotation_north_arrow(
    location = "br", 
    which_north = "true", 
    style = north_arrow_orienteering(),
    pad_x = unit(0.25, "in"), 
    pad_y = unit(0.25, "in"),
    height = unit(1, "cm"),
    width = unit(1, "cm")
  ) +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue", color = NA),
    panel.grid.major = element_line(color = "white", size = 0.5)
  )