#!/usr/bin/env Rscript

# Ładowanie bibliotek
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(sf)
  library(terra)
  library(rnaturalearth)
  library(dplyr)
  library(ggnewscale)
})

# 1. Definicja argumentów
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Wyczyszczony plik CSV z GBIF", metavar="character"),
  make_option(c("-d", "--contours"), type="character", default=NULL, 
              help="Plik rastrowy contours (rds)", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Nazwa pliku wyjściowego (np. mapa.svg)", metavar="character"),
  make_option(c("-s", "--species"), type="character", default="Gatunek", 
              help="Nazwa gatunku do etykiety", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$contours) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("Brak wymaganych argumentów --input, --contours lub --output", call.=FALSE)
}

# 2. Wczytywanie danych
cat("Wczytywanie danych i generowanie warstwic...\n")
gbif_alps <- read.csv(opt$input, sep = "," , header = TRUE, quote = "\"")
print(colnames(gbif_alps)) # Sprawdzenie nazw kolumn
contours_sf <- readRDS(opt$contours)
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")



# 3. Tworzenie mapy
cat("Generowanie grafiki...\n")
map_plot <- ggplot() +
  # Zwiększamy grubość linii (size) do 0.4 i dodajemy tło dla czytelności
  geom_sf(data = contours_sf, aes(color = level), size = 0.4) + 
  scale_color_gradientn(
    colors = c("#2d5e39", "#e6cc80", "#b35900", "#523110"), # Bardziej kontrastowe kolory "górskie"
    name = "Wysokość (m)"
  ) +
  
  new_scale_color() +
  
  # Dodajemy granice państw pod punkty, ale nad warstwice
  geom_sf(data = world, fill = NA, color = "black", alpha = 0.5, size = 0.5) +
  
  # Punkty obecności - zwiększony rozmiar i obrys
  geom_point(data = gbif_alps, 
             aes(x = decimalLongitude, y = decimalLatitude), 
             color = "red", size = 1.2, alpha = 0.8) +
  
  coord_sf(xlim = c(4, 18), ylim = c(43.5, 48.5), expand = FALSE) +
  theme_bw() + # theme_bw lepiej odcina krawędzie mapy niż minimal
  theme(panel.background = element_rect(fill = "#f0f8ff")) + # Jasnoniebieskie tło (ocean/jeziora)

  # Stylistyka
  labs(title = paste("Rozmieszczenie:", opt$species),
       subtitle = paste("Liczba stanowisk:", nrow(gbif_alps)),
       x = "Długość geogr.", y = "Szerokość geogr.") +
  theme(legend.position = "right")



# 4. Zapis do pliku
cat("Zapisywanie mapy do:", opt$output, "\n")
ggsave(opt$output, plot = map_plot, width = 10, height = 8, device = "svg")
