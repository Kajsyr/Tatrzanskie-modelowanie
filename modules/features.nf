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

    emit:
    metadata = species_dataset_ch // Początkowy dataset (metadane)
    features = all_features       // Plik: gbifID + cechy (obecności + pseudo)
}