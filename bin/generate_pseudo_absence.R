#!/usr/bin/env Rscript

# Instalacja wymaganych paczek, jeśli ich nie ma:
# install.packages(c("sf", "optparse"))

library(sf)
library(optparse)

#Konfiguracja argumentów linii komend
option_list <- list(
  make_option(c("-p", "--shp"), type="character", default=NULL, 
              help="Ścieżka do pliku shapefile (.shp)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="output.csv", 
              help="Nazwa pliku wyjściowego [default= %default]", metavar="character"),
  make_option(c("-n", "--count"), type="integer", default=10000, 
              help="Liczba punktów do wygenerowania [default= %default]", metavar="number")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)



# Wczytywanie danych
message(paste("Wczytywanie warstwy:", opt$shp))
# Wczytywanie warstwy
shape <- st_read(opt$shp, quiet = TRUE)

# Bezpieczne sprawdzanie CRS
# Jeśli to nie jest WGS84 (EPSG 4326), wykonaj transformację
current_epsg <- st_crs(shape)$epsg

if (is.na(st_crs(shape)) || !identical(current_epsg, 4326L)) {
    message("Transformacja CRS do WGS84 (EPSG:4326)...")
    shape <- st_transform(shape, 4326)
}

# Naprawa geometrii 
shape <- st_make_valid(st_union(shape))

# Losowanie punktów
message(paste("Generowanie", opt$count, "punktów..."))
points <- st_sample(shape, size = opt$count, type = "random")

# Konwersja na ramkę danych
coords <- st_coordinates(points)
points_df <- as.data.frame(coords)

# Tworzenie finalnego zestawu
final_df <- data.frame(
    gbifID = paste0("pseudo_", seq_len(nrow(points_df))),
    decimalLongitude = points_df$X,
    decimalLatitude = points_df$Y
)

# Zapis do CSV
write.csv(final_df, opt$output, row.names = FALSE)