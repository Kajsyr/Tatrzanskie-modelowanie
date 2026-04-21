#!/usr/bin/env Rscript
# -----------------------------------------------------------
# extract_env.R
# Ekstrakcja zmiennych środowiskowych (DEM + BIOCLIM)
#
# Flagi:
#   --input   plik CSV z punktami GBIF
#   --output  plik wynikowy CSV
#   --dem     raster DEM (jeden plik .tif)
#   --bioclim katalog z plikami BIOCLIM
#
# Działanie:
#   1. Wczytuje punkty i tworzy obiekt przestrzenny.
#   2. Ekstrahuje wysokość i zmienne terenowe.
#   3. Ekstrahuje zmienne BIOCLIM.
#   4. Zapisuje wynik jako CSV.
# -----------------------------------------------------------

# Ładowanie bibliotek w trybie cichym
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(terra)
})

#Definicja argumentów
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Ścieżka do pliku CSV z GBIF", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Ścieżka do zapisu wyjściowego CSV", metavar="character"),
  make_option(c("--dem"), type="character", default=NULL, 
              help="Ścieżka do pliku DEM ", metavar="character"),
  make_option(c("--bioclim"), type="character", default=NULL, 
              help="Ścieżka do worldcim ", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


# Wczytywanie danych
gbif_alps <- read.delim(opt$input, sep = ",", header = TRUE, 
                        quote = "\"", check.names = FALSE)


#Przygotowanie punktów i rastrów
points_lonlat <- vect(gbif_alps, 
                geom = c("decimalLongitude", "decimalLatitude"), 
                crs = "EPSG:4326")


#Wczytanie dem i Bioclim
dem <- rast(opt$dem)
bioclim_stack <- rast(list.files(opt$bioclim, full.names = TRUE))
    
# Projektowanie pod punkty
points_proj <- project(points_lonlat, crs(dem))


#Wyciagamy cechy terenu
gbif_alps$elevation <- terra::extract(dem, points_proj)[,2]


#Liczymy wartośći rónych wskążników terenu
terrain_vars <- c("slope", "aspect", "TPI", "TRI", "roughness", "flowdir")
for (v in terrain_vars) {
  lyr <- terrain(dem, v = v, unit = "degrees")
  gbif_alps[[v]] <- terra::extract(lyr, points_proj)[,2]
}

#Ekstrakcja BIOCLIM
bio_ext <- terra::extract(bioclim_stack, points_lonlat)
bio_values <- bio_ext[, -1]   # usuwamy ID


#Łączenie wszystkich  i finalny zapis
final_data <- cbind(gbif_alps, bio_values)
write.csv(final_data, opt$output, row.names = FALSE)

