process DOWNLOAD_SNOW {
    //storeDir "${params.cache_gee}/snow"
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
    //storeDir "${params.cache_gee}/NDVI"
    conda "${params.project_envs}/geedim.yml"

    output:
    path "*.tif"

    script:
    """
    geedim download \
        -i projects/lofty-complex-481222-k8/assets/ndvi_alpy_2023_static \
        --scale 100

    mv *.tif NDVI.tif
    """
}


process POINT_EXTRACTION {
    publishDir "${params.output_dir}/gee/", mode: 'copy'
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

process MERGE_CSV {
    publishDir "${params.output_dir}/gee/", mode: 'copy'


    output:
    path "*.csv"

    script:
    """

    """
}

workflow GEE {
    take:
        species_ch

    main:
        snow_ch = DOWNLOAD_SNOW()
        ndvi_ch = DOWNLOAD_NDVI()

        //take all tif to one chanell
        gee_tif_ch = snow_ch.mix(ndvi_ch)

        //extract point value for each individual
        gee_csv_ch = POINT_EXTRACTION(gee_tif_ch, species_ch)

        //merge all cvs base on gbif id
        //merge_csv_ch =MERGE_CSV(gee_csv_ch)
    emit:
        gee_feature_ch = gee_csv_ch

}