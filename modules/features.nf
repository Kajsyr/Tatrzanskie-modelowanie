process PREPARE_PRESENCE_POINTS {
    conda "${params.project_envs}/rtatra.yml"
    
    input:
    path gbif_raw

    output:
    path "presence_coords.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    df <- read.csv("${gbif_raw}")
    
    # Wybieramy tylko niezbędne kolumny i ujednolicamy nazwy
    df_clean <- df %>%
      select(gbifID, decimalLongitude, decimalLatitude)
    
    write.csv(df_clean, "presence_coords.csv", row.names = FALSE)
    """
}

process BIOCLIM_CONCATENATE {
    conda "${params.project_envs}/rtatra.yml"
    publishDir "${params.output_features}/bioclim/", mode: 'copy'
    //storeDir "${params.cache_dir}/features"

    input:
    path bioclim_data

    output:
    path "bioclim_all.csv"

    script:
    """
    awk 'FNR==1 && NR!=1{next;}{print}' *.csv > bioclim_all.csv
    """
}

process BIOCLIM_EXTRACT_DATA {
    conda "${params.project_envs}/rtatra.yml"
    publishDir "${params.output_features}/bioclim/", mode: 'copy'
    //storeDir "${params.cache_dir}/features"
    
    input:
    path gbif_species
    path bioclim_dir
 

    output:
    path "results/*"

    script:
    """
    mkdir results/
    extract_env.R \
      --input ${gbif_species} \
     --output results/${gbif_species.simpleName}.csv \
     --bioclim ${bioclim_dir}/
    """
}




process GEE_POINT_EXTRACTION {
    publishDir "${params.output_features}/gee/", mode: 'copy'
    conda "${params.project_envs}/rtatra.yml"
    tag "${raster.baseName}"

    input:
        path raster
        path gatunki

    output:
    path "*.csv"

    script:
    """
    extract_from_raster.R \
        -i $gatunki \
        -r $raster \
        -o ${raster.baseName}.csv
    """
}


process GEE_JOIN_CSV {
    publishDir "${params.output_features}/gee/", mode: 'copy'
    conda "${params.project_envs}/rtatra.yml"

    input:
    path files 

    output:
    path "all_gee.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(purrr)

    # Pobranie listy plików przekazanych przez Nextflow
    file_list <- list.files(pattern = "*.csv")
    
    # Funkcja wczytująca i łącząca dane
    final_df <- file_list %>%
      map(~read.csv(.x, header = TRUE)) %>%
      reduce(inner_join, by = "gbifID")

    # Zapisanie wyniku
    write.csv(final_df, "all_gee.csv", row.names = FALSE)
    """
}

process CONCATENATE {
    publishDir "${params.output_features}", mode: 'copy'
    //storeDir "${params.cache_dir}/features"
    conda "${params.project_envs}/rtatra.yml"

    input:
    path bioclim_data
    path gee_data

    output:
    path "all_features.csv"

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)

    # Wczytywanie danych
    bioclim_df <- read.csv("${bioclim_data}")
    gee_df <- read.csv("${gee_data}")

    # Konwersja klucza na typ znakowy (character), aby umożliwić połączenie 
    # identyfikatorów numerycznych z GBIF i tekstowych z pseudoabsencji
    bioclim_df <- bioclim_df %>%
    mutate(gbifID = as.character(gbifID))

    gee_df <- gee_df %>%
    mutate(gbifID = as.character(gbifID))

    # Wykonanie złączenia
    final_df <- bioclim_df %>%
    inner_join(gee_df, by = "gbifID")

    # Zapisanie wyniku
    write.csv(final_df, "all_features.csv", row.names = FALSE)
    """
}

process STACK_RASTERS {
    publishDir "${params.output_features}/rasters/", mode: 'copy'
    conda "${params.project_envs}/rtatra.yml"

    input:
    path "gee_layers/*"
    path "bioclim_layers/*"
    path alps_shp

    output:
    path "all_layers_stack.tif"

    script:
    """
    #!/usr/bin/env Rscript
    library(terra)

    # Ustawienie limitu pamięci RAM dla terra (np. max 4GB na ten proces)
    # Zapobiega to agresywnemu alokowaniu pamięci przez R
    gdalCache(10000) 

    # 1. Pobranie listy plików TIFF
    files <- list.files(c("gee_layers", "bioclim_layers"), pattern = "\\\\.tif\$", full.names = TRUE, recursive = TRUE)
    
    if (length(files) == 0) {
        stop("Nie znaleziono żadnych plików TIFF do połączenia!")
    }

    # 2. Wczytanie i przygotowanie maski Shapefile
    mask_geom <- vect("${alps_shp.name}")

    # Stworzenie tymczasowego katalogu na docięte warstwy (zwolnienie RAMu)
    dir.create("tmp_cropped", showWarnings = FALSE)

    # 3. Przygotowanie warstwy referencyjnej (pierwszy plik)
    first_raster <- rast(files[1])
    if (crs(mask_geom) != crs(first_raster)) {
        mask_geom <- project(mask_geom, crs(first_raster))
    }
    
    # Docinanie warstwy referencyjnej i natychmiastowy zapis na dysk
    ref_raster <- crop(first_raster, mask_geom)
    ref_raster <- mask(ref_raster, mask_geom)
    
    ref_path <- "tmp_cropped/ref_layer.tif"
    writeRaster(ref_raster, ref_path, overwrite = TRUE)
    
    # Odświeżamy uchwyt do pliku, żeby R czytał go z dysku, a nie z RAMu
    ref_raster <- rast(ref_path)

    # Lista do przechowywania ŚCIEŻEK do plików, a nie samych obiektów w RAMie
    processed_paths <- list(ref_path)

    # 4. Przetwarzanie sekwencyjne (plik po pliku) z czyszczeniem pamięci
    if (length(files) > 1) {
        for (i in 2:length(files)) {
            message(paste("Przetwarzanie warstwy", i, "z", length(files), ":", files[i]))
            
            # Wczytaj bieżący raster
            current_raster <- rast(files[i])
            
            # Wstępne przycięcie do obwiedni (crop), żeby resample nie przetwarzał całego globu
            current_raster <- crop(current_raster, mask_geom)
            
            # Dopasowanie siatki (resample) i maskowanie
            if (!ext(ref_raster) == ext(current_raster) || !all(res(ref_raster) == res(current_raster))) {
                current_raster <- resample(current_raster, ref_raster, method = "bilinear")
            }
            current_raster <- mask(current_raster, mask_geom)
            
            # Zapis dociętego pliku na dysk
            out_path <- paste0("tmp_cropped/layer_", i, ".tif")
            writeRaster(current_raster, out_path, overwrite = TRUE)
            
            # Zapisujemy tylko ścieżkę tekstową
            processed_paths[[i]] <- out_path
            
            # KLUCZOWE: Usuwamy obiekty i wymuszamy Garbage Collector w R, aby zwolnić RAM
            rm(current_raster)
            gc()
        }
    }

    # 5. Łączenie na poziomie dysku (Wirtualny Stack)
    # Podając listę ścieżek tekstowych, rast() nie ładuje danych do RAMu,
    # tworzy jedynie strukturę metadanych wskazującą na pliki
    final_paths <- unlist(processed_paths)
    stacked_raster <- rast(final_paths)

    # 6. Finalny zapis (terra zapisze to blok po bloku, oszczędzając RAM)
    writeRaster(stacked_raster, "all_layers_stack.tif", overwrite = TRUE)
    
    # Czyszczenie plików tymczasowych
    unlink("tmp_cropped", recursive = TRUE)
    """
}



workflow FEATURES {
    take:
        species_dataset_ch // Surowy plik GBIF
        pseudo_absences_ch // Plik z pseudoabsencjami (już z gbifID, lat, lon)
        bioclim_ch
        gee_tif_ch
    
    main:
    // 1. Przygotowujemy czyste punkty obecności
    presence_coords_ch = PREPARE_PRESENCE_POINTS(species_dataset_ch)

    // 2. Łączymy obecności i pseudoabsencje w JEDEN plik wejściowy dla ekstrakcji
    // Używamy prostego cat (z pominięciem nagłówka drugiego pliku)
    combined_points_ch = presence_coords_ch.concat(pseudo_absences_ch)
        .collect()
        .map { files -> 
            def out = file("combined_points.csv")
            out.text = files[0].text + files[1].readLines()[1..-1].join('\n')
            return out
        }

    // 3. Ekstrakcja Bioclim (na połączonym pliku)
    bioclim_raw = BIOCLIM_EXTRACT_DATA(combined_points_ch, bioclim_ch)
    bioclim_final = BIOCLIM_CONCATENATE(bioclim_raw.collect())

    // 4. Ekstrakcja GEE (na połączonym pliku)
    gee_raw = GEE_POINT_EXTRACTION(gee_tif_ch, combined_points_ch)
    gee_final = GEE_JOIN_CSV(gee_raw.collect())

    bioclim_final.view()
    gee_final.view()

    // 5. Ostateczne połączenie cech (Join po gbifID)
    all_features = CONCATENATE(bioclim_final, gee_final)

    shapefile_ch = Channel.fromPath(params.alp_shapefile)
    
    // Łączymy wszystkie warstwy rastrów w jeden wielowarstwowy TIFF, przycinając do granic Alp

    stacked_tiff_ch = STACK_RASTERS(
        gee_tif_ch.collect(), 
        bioclim_ch.collect(), 
        shapefile_ch
    )

    emit:
    metadata = species_dataset_ch // Początkowy dataset (metadane)
    features = all_features       // Plik: gbifID + cechy (obecności + pseudo)
    stacked_raster = stacked_tiff_ch    // NOWOŚĆ: Wielowarstwowy plik TIFF
}