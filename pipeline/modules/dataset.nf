nextflow.enable.dsl=2

// Parametry wejściowe
params.species_file = "${projectDir}/input/FLORA_UJ_TPN.xlsx"
params.alp_shapefile = "${projectDir}/input/Alpy/BiogeoRegions2016_Alpine_ExportFeatures.shp"

process EXTRACT_SPECIES {
    executor 'local' // To zadanie jest szybkie, nie wymaga klastra
    
    input:
    path xlsx

    output:
    stdout emit: species_list

    script:
    """
    #!/usr/bin/env Rscript
    library(readxl)
    df <- read_excel("$xlsx", sheet = "species")
    cat(df\$Species, sep="\\n")
    """
}



process GET_GBIF_DATA {
    conda '${baseDir}/envs/rtatry.yml'
    // Używamy nazwy gatunku w podfolderach, żeby utrzymać porządek
    storeDir "${params.output_dir}/cache_db/raw_data/${species_name.replaceAll(' ', '_')}"
    publishDir "${params.output_dir}/raw_data/${species_name.replaceAll(' ', '_')}", mode: 'copy'
    
    // Ważne: GBIF ma limity jednoczesnych pobrań (zwykle 3), 
    // ograniczamy liczbę procesów działających w tym samym czasie.
    maxForks 3 

    input:
    val species_name

    output:
    // Tworzymy plik o nazwie specyficznej dla gatunku
    path "${species_name.replaceAll(' ', '_')}_gbif.csv", emit: raw_csv

    script:
    """
    #!/usr/bin/env Rscript
    library(rgbif)
    library(dplyr)

    # 1. Pobierz klucz dla konkretnego gatunku
    key <- name_backbone(name = "${species_name}")\$usageKey
    
    if (is.null(key)) {
        stop("Nie znaleziono klucza dla: ${species_name}")
    }

    wkt_polygon <- "POLYGON((4 43.5, 18 43.5, 18 48.5, 4 43.5))"

    # 2. Wyślij zapytanie dla JEDNEGO klucza
    gbif_download_key <- occ_download(
      pred("taxonKey", key),
      pred_gte("year", 1990),
      pred_lt("coordinateUncertaintyInMeters", 20),
      pred("hasCoordinate", TRUE),
      pred_within(wkt_polygon),
      pred("hasGeospatialIssue", FALSE),
      format = "SIMPLE_CSV",
      user = "${params.gbif_user}",
      pwd  = "${params.gbif_pwd}",
      email = "${params.gbif_email}"
    )

    occ_download_wait(gbif_download_key)
    dwn <- occ_download_get(gbif_download_key)
    df <- occ_download_import(dwn)
    
    # Zapisz wynik
    filename <- paste0(gsub(" ", "_", "${species_name}"), "_gbif.csv")
    write.table(df, filename, sep = "\\t", row.names = FALSE)
    """
}

process PROCESS_DEM {
    publishDir "${params.output_dir}/spatial", mode: 'copy'

    input:
    path shp_dir // folder zawierający pliki .shp, .shx, .dbf

    output:
    path "dem_alps.tif", emit: dem_tif
    path "shp_alps_simply.rds", emit: shp_rds

    script:
    """
    #!/usr/bin/env Rscript
    library(sf)
    library(elevatr)
    library(terra)
    library(rmapshaper)

    shp_alps <- st_read("$shp_dir")
    shp_alps_simply <- ms_simplify(shp_alps, keep = 0.01, keep_shapes = TRUE)
    saveRDS(shp_alps_simply, "shp_alps_simply.rds")

    dem <- get_elev_raster(locations = shp_alps_simply, z = 10, clip = "locations") 
    writeRaster(rast(dem), "dem_alps.tif", overwrite = TRUE)
    """
}

process FILTER_AND_PLOT {
    publishDir "${params.output_dir}/plots", mode: 'copy'

    input:
    path raw_csv
    path dem_tif
    path shp_rds

    output:
    path "map_alps_all.svg"
    path "species_counts.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(sf)
    library(terra)
    library(ggplot2)
    library(CoordinateCleaner)
    library(rnaturalearth)

    # 1. Czyszczenie danych GBIF
    gbif_download <- read.delim("$raw_csv", sep = "\\t", header = TRUE, quote = "")
    colnames(gbif_download) <- tolower(colnames(gbif_download))
    
    # CoordinateCleaner wymaga konkretnych nazw
    gbif_clean <- gbif_download %>%
      rename(decimalLongitude = decimallongitude, decimalLatitude = decimallatitude) %>%
      cc_val() %>% cc_zero() %>% cc_cap() %>% cc_urb() %>% cc_gbif()

    # 2. Przestrzenne filtrowanie do obszaru Alp
    shp_alps_simply <- readRDS("$shp_rds")
    df_sf <- st_as_sf(gbif_clean, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
    df_sf <- st_transform(df_sf, st_crs(shp_alps_simply))
    
    inside_vector <- st_within(df_sf, shp_alps_simply, sparse = FALSE)
    gbif_alps <- gbif_clean[apply(inside_vector, 1, any), ]

    # 3. Statystyki
    counts <- gbif_alps %>% group_by(species) %>% summarise(count = n()) %>% arrange(desc(count))
    write.csv(counts, "species_counts.csv", row.names = FALSE)

    # 4. Mapowanie
    dem_alps <- rast("$dem_tif")
    contours <- terra::as.contour(dem_alps, levels = seq(0, 5000, 200))
    contours_sf <- st_as_sf(contours)
    world <- ne_countries(scale = "medium", returnclass = "sf")

    map_alps_all <- ggplot() +
      geom_sf(data = contours_sf, aes(color = level), size = 0.3) +
      scale_color_gradientn(colors = c("darkgreen", "yellow", "orange", "brown")) +
      ggnewscale::new_scale_color() +
      geom_sf(data = world, alpha=0.2) +
      geom_point(data = gbif_alps, aes(x = decimalLongitude, y = decimalLatitude, color = species), size = 0.5) +
      coord_sf(xlim = c(4, 18), ylim = c(43.5, 48.5), expand = FALSE) +
      theme_minimal()

    ggsave("map_alps_all.svg", plot = map_alps_all, width = 10, height = 8)
    """
}

process DOWNLOAD_BIOCLIM {
    conda '${baseDir}/envs/rtatry.yml'
    // Używamy storeDir, aby nie pobierać ciężkich plików rastrowych przy każdym uruchomieniu
    storeDir "${params.output_dir}/cache_db/worldclim"

    output:
    path "wc2.1_2.5m/*.tif", emit: bioclim_tifs

    script:
    """
    #!/usr/bin/env Rscript
    library(geodata)
    # Pobiera dane BIO (1-19) dla całego świata w rozdzielczości 2.5 minuty
    worldclim_global(var = "bio", res = 2.5, path = ".")
    """
}

process EXTRACT_ENVIRONMENTAL_DATA {
    conda '${baseDir}/envs/rtatry.yml'
    publishDir "${params.output_dir}/final_data", mode: 'copy'
    
    input:
    path "gbif_species/*.csv"
    path "bioclim_dir/*.tif"
    path dem_tif

    output:
    path "gbif_alps_env_final.csv", emit: final_data

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(terra)

    # 1. Łączenie pobranych gatunków
    files <- list.files("gbif_species", full.names = TRUE)
    gbif_all <- bind_rows(lapply(files, function(f) read.delim(f, sep="\\t")))
    colnames(gbif_all) <- tolower(colnames(gbif_all))

    # 2. Przygotowanie punktów i rastrów
    points_lonlat <- vect(gbif_all, geom = c("decimallongitude", "decimallatitude"), crs = "EPSG:4326")
    
    dem <- rast("$dem_tif")
    bioclim_stack <- rast(list.files("bioclim_dir", full.names = TRUE))
    
    # Projektowanie punktów pod DEM (jeśli DEM nie jest w stopniach)
    points_proj <- project(points_lonlat, crs(dem))

    # 3. Ekstrakcja TEREN (Elevation, Slope, etc.)
    cat("Analiza terenu...\\n")
    gbif_all\$elevation <- terra::extract(dem, points_proj)[,2]
    
    terrain_vars <- c("slope", "aspect", "TPI", "TRI", "roughness")
    for (v in terrain_vars) {
        lyr <- terrain(dem, v = v, unit = "degrees")
        gbif_all[[v]] <- terra::extract(lyr, points_proj)[,2]
    }

    # 4. Ekstrakcja BIOCLIM (BIO1-BIO19)
    cat("Ekstrakcja Bioclim...\\n")
    # Bioclim zazwyczaj jest w EPSG:4326, więc używamy points_lonlat
    bio_ext <- terra::extract(bioclim_stack, points_lonlat)
    
    # Porządkowanie nazw (upewnij się, że kolejność plików jest poprawna)
    # Pliki z geodata mają nazwy typu wc2.1_2.5m_bio_1.tif
    bio_values <- bio_ext[,-1]
    
    # 5. Łączenie i finalny zapis
    final_data <- cbind(gbif_all, bio_values)
    write.csv(final_data, "gbif_alps_env_final.csv", row.names = FALSE)
    """
}

process SOLAR_ANALYSIS_SAGA {
    conda '${baseDir}/envs/rtatry.yml'
    publishDir "${params.output_dir}/solar", mode: 'copy'

    input:
    path "gbif_all.csv"   // Połączone dane z poprzednich kroków
    path dem_tif

    output:
    path "gbif_with_solar.csv", emit: solar_data

    script:
    """
    #!/usr/bin/env Rscript
    library(terra)
    library(RSAGA)
    library(dplyr)

    # 1. Przygotowanie danych
    gbif <- read.csv("gbif_all.csv")
    dem  <- rast("$dem_tif")
    
    # SAGA wymaga formatu .sdat (Gdal raster) lub .sgrd
    # Zapisujemy DEM w formacie SAGA
    writeRaster(dem, "dem_input.sdat", overwrite=TRUE)

    # 2. Konfiguracja środowiska SAGA
    # RSAGA spróbuje automatycznie znaleźć binarki zainstalowane przez Conda
    env <- rsaga.env()

    # Definiujemy środek obszaru dla obliczeń (lub wyciągamy z danych)
    lat_center <- mean(gbif\$decimallatitude, na.rm=TRUE)
    lon_center <- mean(gbif\$decimallongitude, na.rm=TRUE)

    # 3. Uruchomienie modułu SAGA: Potential Incoming Solar Radiation
    # Uwaga: To może chwilę potrwać zależnie od rozdzielczości DEM
    rsaga.pisr2(
        in.dem      = "dem_input.sgrd",
        in.latitude  = lat_center,
        in.longitude = lon_center,
        out.direct   = "solar_direct.sgrd",
        out.diffuse  = "solar_diffuse.sgrd",
        env          = env
    )

    # 4. Odczyt wyników i ekstrakcja do punktów
    direct_rast  <- rast("solar_direct.sdat")
    diffuse_rast <- rast("solar_diffuse.sdat")
    total_solar  <- direct_rast + diffuse_rast

    points <- vect(gbif, geom = c("decimallongitude", "decimallatitude"), crs = "EPSG:4326")
    # Upewnij się, że CRS punktów zgadza się z rastrem SAGA
    points_proj <- project(points, crs(total_solar))

    solar_values <- terra::extract(total_solar, points_proj)
    
    # 5. Dodanie wyników do tabeli
    gbif\$solar_radiation <- solar_values[,2]

    write.csv(gbif, "gbif_with_solar.csv", row.names = FALSE)
    """
}

process GENERATE_PSEUDOABSENCE {
    conda '${baseDir}/envs/rtatry.yml'
    publishDir "${params.output_dir}/presence_absence", mode: 'copy'

    input:
    path "gbif_all.csv"
    path dem_tif

    output:
    path "gbif_pa_combined.csv", emit: pa_data

    script:
    """
    #!/usr/bin/env Rscript
    library(terra)
    library(dplyr)

    # 1. Wczytanie danych i rastra (jako maski)
    presences <- read.csv("gbif_all.csv")
    dem <- rast("$dem_tif")

    # 2. Tworzenie punktów tła (Background points)
    # Zazwyczaj generuje się ich tyle samo co obecności (1:1) 
    # lub stałą liczbę (np. 10 000)
    n_absences <- nrow(presences)
    
    set.seed(123) # Dla powtarzalności
    bg_points <- spatSample(dem, 
                            size = n_absences, 
                            method = "random", 
                            na.rm = TRUE, 
                            as.points = TRUE)

    # 3. Przygotowanie tabeli tła
    # Musi mieć te same kolumny co tabela obecności
    bg_df <- as.data.frame(bg_points, geom="XY") %>%
             rename(decimallongitude = x, decimallatitude = y) %>%
             mutate(species = "pseudo_absence",
                    occurrenceStatus = "ABSENT")

    # 4. Dodanie kolumny 'pa' (1 dla obecności, 0 dla tła)
    presences_final <- presences %>% mutate(pa = 1)
    absences_final <- bg_df %>% 
                      select(colnames(presences)) %>% 
                      mutate(pa = 0)

    # 5. Połączenie
    combined <- bind_rows(presences_final, absences_final)
    write.csv(combined, "gbif_pa_combined.csv", row.names = FALSE)
    """
}

process SPATIAL_THINNING {
    conda '${baseDir}/envs/rtatry.yml'
    publishDir "${params.output_dir}/thinned_data", mode: 'copy'

    input:
    path "gbif_all.csv"

    output:
    path "gbif_thinned.csv", emit: thinned_data

    script:
    """
    #!/usr/bin/env Rscript
    library(spThin)
    library(dplyr)

    df <- read.csv("gbif_all.csv")

    # Przerzedzanie dla każdego gatunku osobno
    thinned_results <- list()
    species_list <- unique(df\$species)

    for (sp in species_list) {
        sp_data <- df %>% filter(species == sp)
        
        # Wykonaj przerzedzanie (np. min. 5 km odległości)
        thin_out <- thin(
            loc.data = sp_data,
            lat.col = "decimallatitude",
            long.col = "decimallongitude",
            spec.col = "species",
            thin.dist = 5, # odległość w kilometrach
            reps = 1,
            locs.thinned.list.return = TRUE,
            write.results = FALSE
        )
        
        # Wybierz pierwszy wynik (największa możliwa liczba punktów)
        thinned_df <- as.data.frame(thin_out[[1]])
        
        # Odzyskaj oryginalne kolumny poprzez join
        sp_thinned <- inner_join(thinned_df, sp_data, 
                                 by = c("latitude" = "decimallatitude", 
                                        "longitude" = "decimallongitude"))
        
        thinned_results[[sp]] <- sp_thinned
    }

    final_df <- bind_rows(thinned_results)
    write.csv(final_df, "gbif_thinned.csv", row.names = FALSE)
    """
}

workflow DATASET {
    // 1. Nazwy gatunków z Excela
    ch_species = EXTRACT_SPECIES(params.species_file)
        .species_list.splitText().map { it.trim() }

    // 2. Pobieranie GBIF (równolegle dla każdego gatunku)
    gbif_ch = GET_GBIF_DATA(ch_species)

    // 3. Przygotowanie zasobów (DEM i BIOCLIM)
    dem_ch = PROCESS_DEM(file(params.alp_shapefile).parent)
    // 1. Zbieramy wszystkie obecności
    all_presences = gbif_ch.raw_csv.collect()

    // 2. Generujemy pseudo-absence
    pa_ch = GENERATE_PSEUDOABSENCE(all_presences, dem_ch.dem_tif)

    // 3. Przekazujemy połączone punkty (P+A) do ekstrakcji środowiskowej
    bioclim_ch = DOWNLOAD_BIOCLIM()
    env_data = EXTRACT_ENVIRONMENTAL_DATA(
        pa_ch.pa_data,
        bioclim_ch.bioclim_tifs.collect(),
        dem_ch.dem_tif
    )

    // 4. SAGA na samym końcu
    SOLAR_ANALYSIS_SAGA(env_data.final_data, dem_ch.dem_tif)
}
