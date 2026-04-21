// Proces łączący CSV
process MERGE_CSV {
    publishDir "${params.output_dir}/features/", mode: 'copy'
    storeDir "${params.cache_dir}/features/", mode: 'copy'

    input:
        path csv_files

    output:
        path "all_species.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(data.table)

    # Pobieramy listę plików
    files <- list.files(pattern = "${csv_files}")

    # Szybkie wczytywanie i łączenie
    data_list <- lapply(files, fread)
    merged_data <- rbindlist(data_list, use.names = TRUE, fill = TRUE)

    # Zapis do pliku
    fwrite(merged_data, "all_species.csv")
    """
}

process DOWNLOAD_BIOCLIM_ALPS {
    conda "${params.project_envs}/rtatra.yml"
    storeDir "${params.cache_dir}/worldclim_alps"

    output:
    path "bioclim/"

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
                    filename = paste0("bioclim/alps_", var_name, ".tif"), 
                    overwrite = TRUE)
    }
    cat("Zakończono. Pliki zapisane w cache.\\n")
    """
}

process DOWNLOAD_CHELSA_ALPS {
    conda "${params.project_envs}/rtatra.yml"
    storeDir "${params.cache_dir}/chelsa_alps"

    output:
    path "alps_chelsa_*.tif"

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
    sol_global <- chelsa(var = "sol", res = 0.5, path = ".", download = TRUE)
    sol_alps <- crop(mean(sol_global), roi_ext) # średnia roczna insolacja
    names(sol_alps) <- "solar_rad"
    bio_alps <- c(bio_alps, sol_alps)

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
    publishDir "${params.output_dir}/features/", mode: 'copy'
    
    input:
    path gbif_species
    path dem_tiff
    path bioclim
     //path chelsa
    

    output:
    path "results/*"

    script:
    """
    mkdir results/
    extract_env.R \
      --input ${gbif_species} \
     --output results/${gbif_species.simpleName}.csv \
     --dem ${dem_tiff} \
     --bioclim ${bioclim}/

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



workflow FEATURES_EXTRACTION {
    take:
        dataset_ch
        dem_tiff
    main:

    // 3. Przekazujemy połączone punkty (P+A) do ekstrakcji środowiskowej
    bioclim_ch = DOWNLOAD_BIOCLIM_ALPS()
    //chelsa_ch = DOWNLOAD_CHELSA_ALPS()

    env_data = EXTRACT_ENVIRONMENTAL_DATA(
        dataset_ch,
        dem_tiff,
        bioclim_ch
        //chelsa_ch.collect(),
        
    )

    //merge_env_data=MERGE_CSV(env_data)
    emit:
        env_data
        //merge_env_data
        
}