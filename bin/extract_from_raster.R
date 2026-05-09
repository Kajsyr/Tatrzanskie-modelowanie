#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(terra)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="CSV z punktami (gbifID, decimalLongitude, decimalLatitude)"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="Ścieżka do wyjściowego CSV"),
  make_option(c("-r", "--raster"), type="character", default=NULL, 
              help="Plik TIFF (jedno- lub wielowarstwowy)")
)

opt <- parse_args(OptionParser(option_list=option_list))

#wczytanie danych
points_df <- read.csv(opt$input, header = TRUE)
r_layer <- rast(opt$raster)

#projekca punktów
points_vect <- vect(points_df, 
                    geom = c("decimalLongitude", "decimalLatitude"), 
                    crs = "EPSG:4326")

if (crs(points_vect) != crs(r_layer)) {
  points_vect <- project(points_vect, crs(r_layer))
}

#ekstrakcja warstw
extracted_data <- terra::extract(r_layer, points_vect)

# Usuwamy pierwszą kolumnę 'ID', która jest dodawana automatycznie przez terra::extract
extracted_values <- extracted_data[, -1, drop = FALSE]

#usuwamy kolumny co mają w nawie MASK
keep_cols <- !grepl("MASK", names(extracted_values), ignore.case = TRUE)

# Zostawiamy tylko wybrane kolumny
extracted_values <- extracted_values[, keep_cols, drop = FALSE]



# Łączymy gbifID z wszystkimi wyciągniętymi warstwami
final_df <- cbind(gbifID = points_df$gbifID, extracted_values)

#zapis
write.csv(final_df, opt$output, row.names = FALSE)

# Informacja w terminalu
cat(paste0("Zakończono sukcesem!\n"))
cat(paste0("Wyciągnięto warstwy: ", paste(names(r_layer), collapse = ", "), "\n"))
cat(paste0("Plik zapisany w: ", opt$output, "\n"))