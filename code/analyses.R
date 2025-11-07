#!/usr/bin/env Rscript

# Set up the working directory
if (!dir.exists(file.path("data", "imported")))
  dir.create(file.path("data", "imported"))
if (!dir.exists(file.path("data", "processed")))
  dir.create(file.path("data", "processed"))
if (!dir.exists(file.path("output", "figure")))
  dir.create(file.path("output", "figure"), recursive = TRUE)
if (!dir.exists(file.path("output", "table")))
  dir.create(file.path("output", "table"))

# Function that check if the package exist, or it will download it, and then load the library
load_or_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
load_or_install("tidyverse")

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
## Download the shapefile data if not downloaded yet
if (
  !dir.exists(file.path("data", "imported", "ICES_areas")) |
    !dir.exists(file.path("data", "imported", "ICES_rectangles")) |
    !dir.exists(file.path("data", "imported", "HELCOM_subbasins"))
) {
  message("Downloading shapefiles for the map")
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
load_or_install("sf")
message("Plotting the map")
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
ggsave(
  filename = file.path("output", "figure", "map.pdf"),
  plot = map,
  width = 7,
  height = 7
)
ggsave(
  filename = file.path("output", "figure", "map.png"),
  plot = map,
  width = 7,
  height = 7
)

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
message("Analysing DNA data")
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
  mutate(survey = paste(month_abb, year)) |>
  # remove fish reads
  filter(Family != "Teleostei") |>
  # filter out the samples with less than 10'000 reads
  group_by(library_ID) |>
  filter(sum(Abundance) > 10000) |>
  ungroup()

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
  filename = file.path("output", "figure", "metazoa_prop.pdf"),
  width = 6,
  height = 4
)

## Focus on the non-metazoans ----
nonmetazoan_data <- fish_data |>
  mutate(diet = case_when(Subdivision == "Metazoa" ~ TRUE, .default = FALSE)) |>
  filter(diet == FALSE)

summary_non_metazoan <-
  nonmetazoan_data |>
  group_by(organism, survey, Division, library_ID) |>
  summarise(rra = sum(rra), .groups = "drop_last") |>
  summarise(rra = mean(rra), .groups = "drop") |>
  group_by(organism, survey) |>
  mutate(rra = rra / sum(rra)) |>
  ungroup()

barplot_non_metazoan_rra <-
  summary_non_metazoan |>
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
  filename = file.path("output", "figure", "barplot_non_metazoan_rra.pdf"),
  width = 7,
  height = 5
)

summary_non_metazoan |>
  mutate(rra = round(rra, 3)) |>
  pivot_wider(names_from = Division, values_from = rra, values_fill = 0) |>
  write_csv(file = file.path("output", "table", "non_metazoan_division.csv"))

## Zoom in alveolata ----
alveolata_data <-
  nonmetazoan_data |>
  # Fill taxonomy so there is no NA
  mutate(
    Subdivision = ifelse(
      is.na(Subdivision),
      paste(Division, "x", sep = "_"),
      Subdivision
    ),
    Class = ifelse(is.na(Class), paste(Subdivision, "x", sep = "_"), Class),
    Order = ifelse(is.na(Order), paste(Class, "x", sep = "_"), Order)
  ) |>
  # Only keep the division Alveolata
  filter(Division == "Alveolata") |>
  # Compute the relative read abundance
  group_by(organism, survey, Order, library_ID) |>
  summarise(rra = sum(rra), .groups = "drop") |>
  group_by(library_ID) |>
  mutate(rra = rra / sum(rra)) |>
  ungroup()

## Frequency of occurrence ----
alveolata_foo <-
  alveolata_data |>
  # Ensure that when the taxa is not detected, the value is 0
  pivot_wider(names_from = Order, values_from = rra, values_fill = 0) |>
  pivot_longer(
    cols = where(is.numeric),
    values_to = "rra",
    names_to = "Order"
  ) |>
  # Logically assign if the rra is larger that 0, it is considered as detected
  mutate(detected = case_when(rra > 0 ~ 1, .default = 0)) |>
  # And summarise the data
  group_by(organism, survey, Order) |>
  summarise(foo = sum(detected) / n_distinct(library_ID), .groups = "drop")

# Plot the frequency of occurrence
alveolata_foo_plot <-
  alveolata_foo |>
  ggplot(aes(x = organism, y = Order, fill = foo * 100)) +
  geom_tile(col = "black") +
  scale_fill_viridis_c(name = "FOO (%)") +
  facet_grid(. ~ survey) +
  coord_fixed() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = NA)
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  filename = file.path("output", "figure", "order_foo.pdf"),
  plot = alveolata_foo_plot,
  width = 5,
  height = 6
)

# Combine average rra with foo for plotting the costello-plot
alveolata_summary <-
  nonmetazoan_data |>
  mutate(
    Subdivision = ifelse(
      is.na(Subdivision),
      paste(Division, "x", sep = "_"),
      Subdivision
    ),
    Class = ifelse(is.na(Class), paste(Subdivision, "x", sep = "_"), Class),
    Order = ifelse(is.na(Order), paste(Class, "x", sep = "_"), Order)
  ) |>
  filter(Division == "Alveolata") |>
  group_by(organism, survey, Order, library_ID) |>
  summarise(rra = sum(rra), .groups = "drop_last") |>
  summarise(rra = mean(rra), .groups = "drop") |>
  group_by(organism, survey) |>
  mutate(rra = rra / sum(rra)) |>
  ungroup() |>
  pivot_wider(names_from = Order, values_from = rra, values_fill = 0) |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "Order",
    values_to = "rra"
  ) |>
  left_join(alveolata_foo, by = join_by(organism, survey, Order))

# Create a data frame with labels for the taxa that have a rra > 0.5 or a foo > 0.5
label_df <-
  alveolata_summary |>
  filter(rra > 0.5 | foo > 0.5) |>
  select(Order) |>
  unique() |>
  arrange(Order) |>
  mutate(label = row_number()) |>
  right_join(
    alveolata_summary |>
      filter(rra > 0.5 | foo > 0.5),
    by = join_by(Order)
  )

load_or_install("ggrepel")

# Plot the costello-plot
costello_plot <- alveolata_summary |>
  ggplot(aes(x = foo, y = rra)) +
  geom_hex(binwidth = c(0.15, 0.15), col = "black") +
  scale_fill_gradient(
    low = "white",
    high = "#3A1772",
    limits = c(1, 25),
    name = "# Order"
  ) +
  geom_text_repel(
    data = label_df,
    aes(label = label),
    min.segment.length = 0,
    force = 5,
    seed = 100
  ) +
  geom_point(data = label_df) +
  facet_grid(organism ~ survey) +
  coord_fixed(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(color = "black", fill = NA)
  ) +
  labs(x = "Frequency of occurrence", y = "RRA")

# create the table that can be used as a legend for the label and taxa
load_or_install("gridExtra")
table_grob <- label_df |>
  select(Order, label) |>
  unique() |>
  tableGrob(rows = NULL)

# Combine the plot to the table
load_or_install("patchwork")
costello_plot_legend <- costello_plot +
  table_grob +
  plot_layout(widths = c(2, .5))
ggsave(
  plot = costello_plot_legend,
  filename = file.path("output", "figure", "costello.pdf"),
  width = 10,
  height = 5
)

# Compute the eDNA index (Winsconsin standardization)
alveolata_rra_avg <- alveolata_data |>
  pivot_wider(names_from = Order, values_from = rra, values_fill = 0) |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "Order",
    values_to = "rra"
  ) |>
  group_by(organism, survey, Order) |>
  summarise(avg_rra = mean(rra), .groups = "drop_last") |>
  mutate(avg_rra = avg_rra / sum(avg_rra)) |>
  ungroup()

edna_plot <-
  alveolata_rra_avg |>
  group_by(Order) |>
  mutate(avg_eDNA = avg_rra / max(avg_rra)) |>
  ungroup() |>
  ggplot(aes(x = organism, y = Order, fill = avg_eDNA)) +
  geom_tile(col = "black") +
  scale_fill_viridis_c(name = "eDNA index") +
  facet_grid(. ~ survey) +
  coord_fixed() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = NA)
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  filename = file.path("output", "figure", "order_edna.pdf"),
  plot = edna_plot,
  width = 5,
  height = 6
)
#rra plot
rra_plot <-
  alveolata_rra_avg |>
  ggplot(aes(x = organism, y = Order, fill = avg_rra)) +
  geom_tile(col = "black") +
  scale_fill_viridis_c(name = "RRA") +
  facet_grid(. ~ survey) +
  coord_fixed() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(color = "black", fill = NA)
  ) +
  labs(x = NULL, y = NULL)
ggsave(
  filename = file.path("output", "figure", "order_rra.pdf"),
  plot = rra_plot,
  width = 5,
  height = 6
)

## Correlation between diet and alveolata ----
diet_df <- fish_data |>
  mutate(
    prey = ifelse(
      Family == "Rotifera_XX",
      "Rotifera",
      ifelse(Family %in% c("Branchiopoda"), "Branchiopoda", Genus)
    )
  ) |>
  filter(
    prey %in%
      c(
        "Rotifera",
        "Acartia",
        "Centropages",
        "Branchiopoda",
        "Eurytemora",
        "Pseudocalanus",
        "Temora"
      ),
    library_ID %in% unique(alveolata_data$library_ID)
  ) |>
  group_by(library_ID) |>
  mutate(rra = Abundance / sum(Abundance)) |>
  ungroup() |>
  group_by(library_ID, organism, survey, prey) |>
  summarise(rra = sum(rra), .groups = "drop") |>
  # Ensure that not detected prey are = 0
  pivot_wider(names_from = prey, values_from = rra, values_fill = 0) |>
  pivot_longer(
    cols = where(is.numeric),
    values_to = "rra",
    names_to = "prey"
  ) |>
  # compute eDNA index
  group_by(prey) |>
  mutate(eDNA = rra / max(rra)) |>
  ungroup()
alveolata_df <-
  alveolata_data |>
  pivot_wider(names_from = Order, values_from = rra, values_fill = 0) |>
  pivot_longer(
    cols = where(is.numeric),
    names_to = "Order",
    values_to = "rra"
  ) |>
  group_by(Order) |>
  mutate(eDNA = rra / max(rra)) |>
  ungroup()
# correlation with rra
alveolata_rra_df <- alveolata_df |>
  select(-eDNA) |>
  pivot_wider(names_from = Order, values_from = rra) |>
  arrange(library_ID)
diet_rra_df <- diet_df |>
  select(-eDNA) |>
  pivot_wider(names_from = prey, values_from = rra) |>
  arrange(library_ID)
# Check that the row are the same
if (!isTRUE(unique(alveolata_rra_df$library_ID == diet_rra_df$library_ID)))
  message("matrices are not matching")
# transform to matrix
alveolata_rra_mat <-
  alveolata_rra_df |>
  select(-c(1:3)) |>
  as.matrix()
diet_rra_mar <-
  diet_rra_df |>
  select(-c(1:3)) |>
  as.matrix()
stopifnot(nrow(alveolata_rra_mat) == nrow(diet_rra_mar))

cor_p_rra <-
  expand.grid(
    prey = colnames(diet_rra_mar),
    alveolata = colnames(alveolata_rra_mat)
  ) |>
  rowwise() |>
  mutate(
    correlation = cor(
      diet_rra_mar[, prey],
      alveolata_rra_mat[, alveolata],
      method = "kendall"
    ),
    p_value = cor.test(
      diet_rra_mar[, prey],
      alveolata_rra_mat[, alveolata],
      method = "kendall"
    )$p.value
  ) |>
  ungroup() |>
  mutate(
    p_adj = p.adjust(p_value, method = "fdr"),
    sign = ifelse(
      p_adj <= 0.001,
      "***",
      ifelse(
        p_adj > 0.001 & p_adj <= 0.01,
        "**",
        ifelse(p_adj > 0.01 & p_adj <= 0.05, "*", "")
      )
    )
  )

cor_plot <- cor_p_rra |>
  ggplot(aes(
    y = reorder(alveolata, abs(correlation)),
    x = reorder(prey, -abs(correlation)),
    fill = correlation,
    label = sign
  )) +
  geom_tile(col = "white") +
  geom_text(col = "white") +
  scale_fill_gradient2(
    low = "#44CF6C",
    mid = "black",
    high = "#FFBB33",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(1, 10)) +
  coord_fixed() +

  theme_void() +
  theme(
    axis.text.y = element_text(color = "black", hjust = 1),
    axis.text.x = element_text(
      color = "black",
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )
ggsave(
  plot = cor_plot,
  filename = file.path("output", "figure", "kendall_correlation.pdf"),
  width = 5.5,
  height = 7
)
