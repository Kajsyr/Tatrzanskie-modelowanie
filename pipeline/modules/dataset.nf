nextflow.enable.dsl=2


process EXTRACT_SPECIES {
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
    conda '${params.project_env}/rtatry.yml'
    // Używamy nazwy gatunku w podfolderach, żeby utrzymać porządek
    storeDir "${params.output_dir}/cache_db/raw_data/${species_name.replaceAll(' ', '_')}"
    publishDir "${params.cache_dir}/raw_data/${species_name.replaceAll(' ', '_')}", mode: 'copy'
    
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

    wkt_polygon <- "POLYGON((4 43.5, 18 43.5, 18 48.5, 4 48.5, 4 43.5))"

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
    storeDir "${params.cache_dir}/dem_alps"
    conda "${params.project_envs}/rtatra.yml"

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

process FILTER {
    publishDir "${params.output_dir}/filter", mode: 'copy'
    conda "${params.project_envs}/rtatra.yml"

    input:
    path raw_csv
    path shp_rds

    output:
    path "*.csv"

    script:
    """
    filter_gbif.R \
        --input ${raw_csv} \
        --output cleaned_${raw_csv.simpleName}.csv \
        --shape ${shp_rds} 
    """
}

process PLOT {
    publishDir "${params.output_dir}/plots", mode: 'copy'
    conda "${params.project_envs}/rtatra.yml"    

    input:
    path raw_csv
    path dem_tif

    output:
    path "map_*.svg"

    script:
    def spec_name = raw_csv.simpleName.replace('_gbif', '').replace('_', ' ')
    """
    plot_gbif.R \
        --input cleaned_${raw_csv.simpleName}.csv \
        --dem dem_alps.tif \
        --species "${spec_name}" \
        --output map_${raw_csv.simpleName}.svg
    """
}


process DOWNLOAD_BIOCLIM_ALPS {
    conda "${params.project_envs}/rtatra.yml"
    // storeDir sprawi, że wycięte rastry zostaną zachowane na stałe
    storeDir "${params.cache_dir}/worldclim_alps"

    output:
    path "alps_*.tif", emit: bioclim_tifs

    script:
    """
    #!/usr/bin/env Rscript
    library(geodata)
    library(terra)

    # 1. Definicja Twojego obszaru (Bounding Box na podstawie POLYGON)
    # xmin: 4, xmax: 18, ymin: 43.5, ymax: 48.5
    roi_ext <- ext(4, 18, 43.5, 48.5)

    # 2. Pobieranie danych globalnych
    # res = 0.5 (ok. 1km) jest optymalne dla obszarów górskich (Alpy)
    cat("Pobieranie danych WorldClim (res=0.5)...\\n")
    bio_global <- worldclim_global(var = "bio", res = 0.5, path = ".")

    # 3. Przycięcie do Twojego obszaru
    cat("Przycinanie do zadanego poligonu...\\n")
    bio_alps <- crop(bio_global, roi_ext)

    # 4. Zapisywanie warstw jako osobne pliki TIFF
    # Pętla przejdzie przez wszystkie 19 zmiennych bioclimatycznych
    for(i in 1:nlyr(bio_alps)) {
        var_name <- names(bio_alps)[i]
        writeRaster(bio_alps[[i]], 
                    filename = paste0("alps_", var_name, ".tif"), 
                    overwrite = TRUE)
    }
    cat("Zakończono. Pliki zapisane w cache.\\n")
    """
}

process DOWNLOAD_CHELSA_ALPS {
    conda "${params.project_envs}/rtatra.yml"
    storeDir "${params.cache_dir}/chelsa_alps"

    output:
    path "alps_chelsa_*.tif", emit: chelsa_tifs

    script:
    """
    #!/usr/bin/env Rscript
    library(geodata)
    library(terra)

    # 1. Definicja obszaru (Bounding Box dla Alp)
    roi_ext <- ext(4, 18, 43.5, 48.5)

    # 2. Pobieranie danych CHELSA
    # var = "bio" pobierze 19 zmiennych bioclimatycznych
    # Uwaga: CHELSA v2.1 jest bardzo dokładna, ale pliki są duże.
    cat("Pobieranie danych CHELSA (bioclim)...\\n")
    
    # Pobieramy zmienne bioclimatyczne (bio)
    # geodata pobierze je do folderu tymczasowego, a my je przytniemy
    bio_global <- chelsa(var = "bio", res = 0.5, path = ".", download = TRUE)

    # 3. Przycięcie do Twojego obszaru
    cat("Przycinanie do zadanego poligonu...\\n")
    bio_alps <- crop(bio_global, roi_ext)

    # 4. Opcjonalnie: Pobranie insolacji (Solar Radiation)
    # Jeśli chcesz dodać nasłonecznienie, odkomentuj poniższe linie:
    # cat("Pobieranie insolacji...\\n")
    # sol_global <- chelsa(var = "sol", res = 0.5, path = ".", download = TRUE)
    # sol_alps <- crop(mean(sol_global), roi_ext) # średnia roczna insolacja
    # names(sol_alps) <- "solar_rad"
    # bio_alps <- c(bio_alps, sol_alps)

    # 5. Zapisywanie warstw
    for(i in 1:nlyr(bio_alps)) {
        var_name <- names(bio_alps)[i]
        writeRaster(bio_alps[[i]], 
                    filename = paste0("alps_chelsa_", var_name, ".tif"), 
                    overwrite = TRUE)
    }
    cat("Zakończono. Pliki zapisane w cache.\\n")
    """
}




process EXTRACT_ENVIRONMENTAL_DATA {
    conda "${params.project_envs}/rtatra.yml"
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
    conda "${params.project_envs}/rtatra.yml"
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
    conda "${params.project_envs}/rtatra.yml"
    publishDir "${params.output_dir}/presence_absence", mode: 'copy'

    input:
    path gbif_csv
    path dem_tif

    output:
    path "gbif_pa_combined.csv", emit: pa_data

    script:
    """
    #!/usr/bin/env Rscript
    library(terra)
    library(dplyr)

    # 1. Wczytanie danych i rastra (jako maski)
    presences <- read.csv("$gbif_csv", sep = "\\t", header = TRUE, quote = "")
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

    // Filtrowanie i mapowanie wystąpień
    filtered_gbif_ch= FILTER(gbif_ch.raw_csv, dem_ch.shp_rds)

    //Mapowanie
    PLOT(filtered_gbif_ch, dem_ch.dem_tif)

    //odfultujemy puste pliki bo cześć odpada wcześniej
    valid_gbif_ch = filtered_gbif_ch.filter { file -> 
        file.readLines().size() > 1 
    }

    // 2. Generujemy pseudo-absence
    //pa_ch = GENERATE_PSEUDOABSENCE(valid_gbif_ch, dem_ch.dem_tif)

    // 3. Przekazujemy połączone punkty (P+A) do ekstrakcji środowiskowej
    bioclim_ch = DOWNLOAD_BIOCLIM_ALPS()
    env_data = EXTRACT_ENVIRONMENTAL_DATA(
        valid_gbif_ch,
        bioclim_ch.bioclim_tifs.collect(),
        dem_ch.dem_tif
    )

    // 4. SAGA na samym końcu
    //SOLAR_ANALYSIS_SAGA(env_data.final_data, dem_ch.dem_tif)
}
