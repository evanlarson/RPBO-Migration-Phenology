options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
  library(dplyr)
})

fig_dir <- "figures_jul_oct"
if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Approximate weather station coordinates (town centers) in BC.
stations <- data.frame(
  station = c(
    "Victoria", "Vancouver", "Courtenay",
    "Port Hardy", "Bella Coola", "Prince Rupert",
    "Kamloops", "Princeton", "Pemberton",
    "Prince George", "Williams Lake", "Smithers"
  ),
  region = c(
    "Southern Coastal BC", "Southern Coastal BC", "Southern Coastal BC",
    "Central Coast", "Central Coast", "Central Coast",
    "South Interior", "South Interior", "South Interior",
    "Central Interior / Sub-Boreal", "Central Interior / Sub-Boreal", "Central Interior / Sub-Boreal"
  ),
  lon = c(
    -123.3656, -123.1207, -124.9936,
    -127.4977, -126.7532, -130.3208,
    -120.3273, -120.5108, -122.8050,
    -122.7497, -122.1417, -127.1685
  ),
  lat = c(
    48.4284, 49.2827, 49.6841,
    50.7212, 52.3721, 54.3150,
    50.6745, 49.4580, 50.3200,
    53.9171, 52.1417, 54.7804
  ),
  label_dx = c(
    1.25, 1.00, 1.10,
    1.00, 1.05, -1.10,
    0.85, 1.15, 1.10,
    1.10, 1.10, -1.05
  ),
  label_dy = c(
    -0.10, 0.35, 0.10,
    0.55, 0.25, 0.35,
    -0.20, -0.35, 0.35,
    0.20, -0.15, 0.20
  ),
  stringsAsFactors = FALSE
)

cache_dir <- "data_cache"
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
}

load_bc_boundary <- function(cache_file) {
  if (file.exists(cache_file)) {
    bc <- readRDS(cache_file)
  } else {
    states <- rnaturalearth::ne_download(
      scale = 50,
      type = "states",
      category = "cultural",
      returnclass = "sf"
    )
    bc <- states %>%
      filter(admin == "Canada", name == "British Columbia") %>%
      sf::st_transform(4326)
    saveRDS(bc, cache_file)
  }
  bc
}

bc_boundary <- load_bc_boundary(file.path(cache_dir, "bc_boundary_ne50.rds"))
base_countries <- rnaturalearth::ne_countries(scale = 50, returnclass = "sf") %>%
  filter(admin %in% c("Canada", "United States of America"))

region_colors <- c(
  "Southern Coastal BC" = "#2a9d8f",
  "Central Coast" = "#457b9d",
  "South Interior" = "#f4a261",
  "Central Interior / Sub-Boreal" = "#e76f51"
)
region_breaks <- names(region_colors)

p <- ggplot() +
  geom_sf(
    data = base_countries,
    fill = "#f7fafc",
    color = "#b8c2cc",
    linewidth = 0.30
  ) +
  geom_sf(
    data = bc_boundary,
    fill = "#e8eef5",
    color = NA
  ) +
  geom_segment(
    data = stations,
    aes(
      x = lon, y = lat,
      xend = lon + label_dx,
      yend = lat + label_dy,
      color = region
    ),
    linewidth = 0.45,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  geom_point(
    data = stations,
    aes(x = lon, y = lat, fill = region),
    shape = 21,
    size = 3.0,
    stroke = 0.5,
    color = "black"
  ) +
  geom_label(
    data = stations,
    aes(
      x = lon + label_dx,
      y = lat + label_dy,
      label = station,
      color = region
    ),
    size = 3.2,
    linewidth = 0.2,
    fill = "white",
    alpha = 0.97,
    label.r = grid::unit(0.08, "lines"),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = region_colors, breaks = region_breaks, labels = region_breaks, name = "Region") +
  scale_color_manual(values = region_colors, breaks = region_breaks, labels = region_breaks, guide = "none") +
  coord_sf(xlim = c(-132.6, -117.4), ylim = c(47.7, 56.2), expand = FALSE) +
  labs(
    title = "British Columbia Weather Stations",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#dbe2ea", linewidth = 0.35),
    legend.position = "bottom",
    axis.title = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14)
  )

out_file <- file.path(fig_dir, "figure_bc_weather_stations_map.png")
ggsave(out_file, p, width = 11, height = 8.5, dpi = 300)

cat("Wrote:", out_file, "\n")
