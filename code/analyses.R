library(tidyverse)

# Import the data ----
asv_table <- read_csv(
  file = file.path("data", "raw", "asv_18S.csv"),
  show_col_types = FALSE
)

metadata <- read_csv(
  file = file.path("data", "raw", "metadata_18S.csv"),
  show_col_types = FALSE
)

# Combine metadata with asv table ----
dna_data <- metadata |>
  full_join(
    asv_table |>
      pivot_longer(
        cols = where(is.numeric), # Pivot  numeric columns (abundance values)
        names_to = "source_material_ID",
        values_to = "Abundance"
      ),
    by = join_by(source_material_ID)
  )

# Get the station coordinates ----
if (!file.exists(file.path("data", "processed", "station_coordinate.csv"))) {
  station_coordinates <- metadata |>
    filter(organism != "WP2") |>
    select(collection_date, lat_lon) |>
    unique() |>
    mutate(month = month(collection_date), year = year(collection_date)) |>
    left_join(
      tibble(month_abb = month.abb, month = 1:12),
      by = join_by(month)
    ) |>
    separate(lat_lon, into = c("lat", "lon"), sep = " ")
  station_coordinates |>
    write_csv(file = file.path("data", "processed", "station_coordinate.csv"))
} else {
  station_coordinates <- read_csv(
    file.path("data", "processed", "station_coordinate.csv"),
    show_col_types = FALSE
  )
}
# Detach the data that is no longer needed
rm(asv_table, metadata)

# Plot map ----
# Get the ices shapefile dataset for the Baltic Sea
## Dowload the shapefile data if not downloaded yet
if (
  !dir.exists(file.path("data", "imported", "ICES_areas")) |
    !dir.exists(file.path("data", "imported", "ICES_rectangles")) |
    !dir.exists(file.path("data", "imported", "HELCOM_subbasins"))
) {
  # 1. ICES Areas
  # The URL for getting the data for ICES areas
  areas_url <- "https://gis.ices.dk/shapefiles/ICES_areas.zip"
  # and the path to save it
  areas_zip_path <- file.path("data", "imported", "ICES_areas.zip")

  # Download
  download.file(areas_url, areas_zip_path, mode = "wb")

  # Unzip and delete the zip file
  areas_unziped_dir <- file.path("data", "imported", "ICES_areas")
  unzip(areas_zip_path, exdir = areas_unziped_dir)
  unlink(areas_zip_path)

  rm(areas_url, areas_zip_path, areas_unziped_dir)

  # 2. ICES rectangles
  # The URL for getting the data for ICES areas
  rect_url <- "https://gis.ices.dk/shapefiles/ICES_rectangles.zip"
  # and the path to save it
  rect_zip_path <- file.path("data", "imported", "ICES_rectangles.zip")

  # Download
  download.file(rect_url, rect_zip_path, mode = "wb")

  # Unzip and delete the zip file
  rect_unziped_dir <- file.path("data", "imported", "ICES_rectangles")
  unzip(rect_zip_path, exdir = rect_unziped_dir)
  unlink(rect_zip_path)

  rm(rect_url, rect_zip_path, rect_unziped_dir)

  # 3. HELCOM sub basins
  # The URL for getting the data for HELCOM sub basins
  HELCOM_url <- "https://gis.ices.dk/shapefiles/HELCOM_subbasins.zip"
  # and the path to save it
  HELCOM_zip_path <- file.path("data", "imported", "HELCOM_subbasins.zip")

  # Download
  download.file(HELCOM_url, HELCOM_zip_path, mode = "wb")

  # Unzip and delete the zip file
  HELCOM_unziped_dir <- file.path("data", "imported", "HELCOM_subbasins")
  unzip(HELCOM_zip_path, exdir = HELCOM_unziped_dir)
  unlink(HELCOM_zip_path)

  rm(HELCOM_url, HELCOM_zip_path, HELCOM_unziped_dir)
}

# Process shapefile data
suppressPackageStartupMessages(library(sf))
baltic_sea_shp <-
  read_sf(list.files(
    file.path("data", "imported", "ICES_areas"),
    pattern = "\\.shp$",
    full.names = TRUE
  )) |>
  filter(SubDivisio %in% 24:32) |>
  group_by(Major_FA) |>
  summarise(geometry = st_union(geometry)) |>
  ungroup()
rectangle_shp <- # Shapefile file with ices rectangle
  read_sf(list.files(
    file.path("data", "imported", "ICES_rectangles"),
    pattern = "\\.shp$",
    full.names = TRUE
  )) |>
  filter(Ecoregion %in% c("Baltic Sea"))

# Combine rectangles with Baltic Sea contour
if (st_crs(baltic_sea_shp) != st_crs(rectangle_shp)) {
  rectangle_shp <- st_transform(rectangle_shp, st_crs(baltic_sea_shp))
}
suppressWarnings(
  fish_shp <- # Combined dataset
    st_intersection(baltic_sea_shp, rectangle_shp)
)

# Transform the station coordinates in shapefile
survey_shp <-
  station_coordinates |>
  mutate(year_month = paste(month_abb, year, sep = " ")) |>
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

# Plot
map <-
  ggplot() +
  geom_sf(
    data = fish_shp,
    fill = NA,
    col = "black"
  ) +
  geom_sf(
    data = survey_shp,
    aes(fill = year_month, shape = year_month),
    color = "black"
  ) +
  # Manual fill scale with labels and colors
  scale_fill_manual(
    values = c("#222E50", "#C9D8AB", "#699051", "#994636"),
    name = ""
  ) +
  scale_shape_manual(values = c(22, 23, 23, 25), name = "") +
  # Custom legend symbols
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text = element_text(color = "black", size = 10),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", fill = "white"),
    plot.background = element_rect(fill = "white", color = NA)
  )
# Save
ggsave(filename = file.path("output", "map.pdf"), plot = map)
ggsave(filename = file.path("output", "map.png"), plot = map)

# Detach packages and files not needed anymore
rm(
  baltic_sea_shp,
  fish_shp,
  map,
  rectangle_shp,
  station_coordinates,
  survey_shp
)
detach("package:sf", unload = TRUE)

# Analyse the DNA data ----
# Transform read count to relative abundance
fish_data <-
  dna_data |>
  filter(organism != "WP2", !is.na(library_ID), Abundance > 0) |>
  group_by(library_ID) |>
  mutate(rra = Abundance / sum(Abundance)) |>
  ungroup() |>
  mutate(month = month(collection_date), year = year(collection_date)) |>
  left_join(
    tibble(month_abb = month.abb, month = 1:12),
    by = join_by(month)
  ) |>
  mutate(survey = paste(month_abb, year))
## Metazoans vs non-metazoans ----
summary_metazoa <- fish_data |>
  mutate(diet = case_when(Subdivision == "Metazoa" ~ TRUE, .default = FALSE)) |>
  group_by(library_ID, diet, organism, survey) |>
  summarise(tot = sum(rra), .groups = "drop")

metazoa_prop <-
  summary_metazoa |>
  filter(diet == F) |>
  group_by(organism, survey) |>
  summarise(avg = mean(tot * 100), sd = sd(tot * 100), .groups = "drop") |>
  mutate(min = pmax(0, avg - sd), max = avg + sd) |>
  ggplot(aes(
    x = organism,
    y = avg,
    min = min,
    max = max,
    fill = survey,
    shape = survey
  )) +
  geom_jitter(
    data = summary_metazoa |> filter(diet == F),
    inherit.aes = FALSE,
    aes(x = organism, y = tot * 100, col = survey),
    position = position_jitterdodge(
      jitter.width = 0.1,
      dodge.width = 0.7,
      jitter.height = 0,
      seed = 999
    ),
    alpha = .7
  ) +
  geom_errorbar(
    col = "black",
    position = position_dodge(width = .7),
    width = .2
  ) +
  geom_point(position = position_dodge(width = .7), size = 3) +
  scale_fill_manual(
    values = c("#222E50", "#C9D8AB", "#699051", "#994636"),
    name = ""
  ) +
  scale_color_manual(
    values = c("#222E50", "#C9D8AB", "#699051", "#994636"),
    name = ""
  ) +
  scale_shape_manual(values = c(22, 23, 23, 25), name = "") +
  labs(x = NULL, y = "Proportion of non-metazoans reads (%; mean ± sd)") +
  theme_bw() +
  theme(panel.grid = element_blank())
ggsave(
  plot = metazoa_prop,
  filename = file.path("output", "metazoa_prop.pdf"),
  width = 6,
  height = 4
)

## Focus on the non-metazoans
nonmetazoan_data <- fish_data |>
  mutate(diet = case_when(Subdivision == "Metazoa" ~ TRUE, .default = FALSE)) |>
  filter(diet == FALSE)
barplot_non_metazoan_rra <-
  nonmetazoan_data |>
  group_by(organism, survey, Division, Subdivision, library_ID) |>
  summarise(rra = sum(rra), .groups = "drop_last") |>
  summarise(rra = mean(rra), .groups = "drop") |>
  group_by(organism, survey) |>
  mutate(rra = rra / sum(rra)) |>
  ungroup() |>
  ggplot(aes(x = organism, y = rra, fill = Division)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ survey) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = NA)
  ) +
  labs(x = NULL, y = "RRA")
ggsave(
  plot = barplot_non_metazoan_rra,
  filename = file.path("output", "barplot_non_metazoan_rra.pdf"),
  width = 7,
  height = 5
)
