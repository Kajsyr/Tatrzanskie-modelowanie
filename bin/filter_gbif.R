#!/usr/bin/env Rscript

# Ładowanie bibliotek w trybie cichym
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(sf)
  library(terra)
  library(CoordinateCleaner)
})

# Definicja argumentów
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="Ścieżka do pliku CSV z GBIF", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Ścieżka do zapisu wyjściowego CSV", metavar="character"),
  make_option(c("-s", "--shape"), type="character", default="shp_alps_simply.rds", 
              help="Ścieżka do pliku RDS z poligonem Alp", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) | is.null(opt$output)) {
  print_help(opt_parser)
  stop("Należy podać argumenty --input oraz --output", call.=FALSE)
}

# Wczytywanie danych
gbif_data <- read.delim(opt$input, sep = ",", header = TRUE, 
                        quote = "\"", check.names = FALSE)

#Tworzenie pliku z logami
log_file <- paste0(opt$output, ".log")

sink(log_file, split = TRUE)
sink(type="message" )



cat("Liczba rekordów przed czyszczeniem:", nrow(gbif_data), "\n")


#CzyszcenieCoordinateCleaner
cat("Czyszczenie danych...\n")
gbif_clean <- gbif_data %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  cc_val(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_cap(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_urb(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_aohi(lon = "decimalLongitude", lat = "decimalLatitude") %>%
  cc_cen(buffer=500) %>%
  cc_gbif()

cat("Liczba rekordów po CoordinateCleaner:", nrow(gbif_clean), "\n")


# Filtracja przestrzenna
cat("Filtracja do obszaru Alp...\n")
shp_alps <- readRDS(opt$shape)
if (inherits(shp_alps, "SpatVector")) shp_alps <- st_as_sf(shp_alps)

# Upewnij się, że shp_alps ma poprawny CRS przed transformacją
if (is.na(st_crs(shp_alps))) st_crs(shp_alps) <- 3035 

# Tworzenie obiektu przestrzennego
df_sf <- st_as_sf(gbif_clean, 
                  coords = c("decimalLongitude", "decimalLatitude"), 
                  crs = 4326)

# Ujednolicenie układów współrzędnych - rzutujemy wszystko na stopnie
shp_alps_4326 <- st_transform(shp_alps, 4326)

# Sprawdzenie punktów wewnątrz poligonu - używamy as.vector dla bezpieczeństwa
cat("Sprawdzanie przecięć...\n")
inside_matrix <- st_intersects(df_sf, shp_alps_4326, sparse = FALSE)
inside_vector <- rowSums(inside_matrix) > 0

gbif_alps <- gbif_clean[inside_vector, ]

cat("Liczba rekordów w shp:", nrow(gbif_alps), "\n")

sink(type= "message")
sink()

write.csv(gbif_alps,paste0(opt$outpu, ".csv") , col.names=TRUE, row.names = FALSE, sep = ",")

