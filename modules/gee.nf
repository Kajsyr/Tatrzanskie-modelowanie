process DOWNLOAD_SNOW {
    storeDir "${params.cache_gee}/snow"
    conda "${params.project_envs}/geedim.yml"

    output:
    path "*.tif"

    script:
    """
    geedim download \
            -i projects/lofty-complex-481222-k8/assets/snow_duration_alpy_2023 \
            --scale 20
    mv *.tif snow_duration.tif
    """
}

process DOWNLOAD_NDVI {
    storeDir "${params.cache_gee}/NDVI"
    conda "${params.project_envs}/geedim.yml"

    output:
    path "*.tif"

    script:
    """
    geedim download \
        -i projects/lofty-complex-481222-k8/assets/ndvi_alpy_2023_static \
        --scale 30

    mv *.tif NDVI.tif
    """
}

process DOWNLOAD_SOLAR_RADIATION {
    storeDir "${params.cache_gee}/solar_radiation"
    conda "${params.project_envs}/geedim.yml"

    output:
    path "*.tif"

    script:
    """
    geedim download \
        -i projects/lofty-complex-481222-k8/assets/solar_radiation_alpy_2023 \
        --scale 30

    mv *.tif solar_radiation.tif
    """
}

process DOWNLOAD_TOPOGRAPHIC_SUITE {
    //storeDir "${params.cache_gee}/NDVI"
    conda "${params.project_envs}/geedim.yml"

    output:
    path "*.tif"

    script:
    """
    geedim download \
        -i projects/lofty-complex-481222-k8/assets/topography_suite_alpy \
        --scale 30

    mv *.tif topografic_suite.tif
    """
}


process POINT_EXTRACTION {
    publishDir "${params.output_gee}", mode: 'copy'
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


process JOIN_CSV {
    publishDir "${params.output_gee}", mode: 'copy'
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


workflow GEE {
    take:
        species_ch

    main:
        snow_ch = DOWNLOAD_SNOW()
        ndvi_ch = DOWNLOAD_NDVI()
        solar_radiation_ch = DOWNLOAD_SOLAR_RADIATION()
        topographic_suite_ch = DOWNLOAD_TOPOGRAPHIC_SUITE()

        //take all tif to one chanell
        gee_tif_ch = snow_ch.mix(ndvi_ch, solar_radiation_ch, topographic_suite_ch)

        //extract point value for each individual
        gee_csv_ch = POINT_EXTRACTION(gee_tif_ch, species_ch)

        //merge all cvs base on gbif id
        merge_csv_ch =JOIN_CSV(gee_csv_ch.collect())
    emit:
        gee_feature_ch = merge_csv_ch

}