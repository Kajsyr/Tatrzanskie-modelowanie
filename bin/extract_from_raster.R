#!/usr/bin/env Rscript
# -----------------------------------------------------------
# extract_named.R
# Ekstrakcja do ramki: long, lat, [nazwa_rastra]
# -----------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(terra)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="CSV z punktami (decimalLongitude, decimalLatitude)"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Ścieżka do wyjściowego CSV"),
  make_option(c("-r", "--raster"), type="character", default=NULL, 
              help="Plik TIFF")
)

opt <- parse_args(OptionParser(option_list=option_list))

# 1. Wczytanie danych i rastra
points_df <- read.csv(opt$input, header = TRUE)
r_layer <- rast(opt$raster)

# 2. Pobranie nazwy warstwy z rastra
raster_name <- names(r_layer)

# 3. Przygotowanie wektora punktów i reprojekcja[cite: 1]
points_vect <- vect(points_df, 
                    geom = c("decimalLongitude", "decimalLatitude"), 
                    crs = "EPSG:4326")

if (crs(points_vect) != crs(r_layer)) {
  points_vect <- project(points_vect, crs(r_layer))
}

extracted_data <- terra::extract(r_layer, points_vect)
extracted_values <- extracted_data[, 2]

# 5. Budowa ramki wynikowej
# Tworzymy nową ramkę, aby uniknąć problemów z indeksowaniem[cite: 1]
final_df <- data.frame(
  decimalLongitude = points_df$decimalLongitude,
  decimalLatitude = points_df$decimalLatitude,
  temp_col = extracted_values
)

# Zmiana nazwy trzeciej kolumny na nazwę rastra
colnames(final_df)[3] <- raster_name
# 6. Zapis[cite: 1]
write.csv(final_df, opt$output, row.names = FALSE)
cat(paste0("Zakończono. Kolumny: long, lat, ", raster_name, "\n"))