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

workflow GEE {
    main:
        snow_ch = DOWNLOAD_SNOW()
        ndvi_ch = DOWNLOAD_NDVI()
        solar_radiation_ch = DOWNLOAD_SOLAR_RADIATION()
        topographic_suite_ch = DOWNLOAD_TOPOGRAPHIC_SUITE()

        //take all tif to one chanell
        gee_tif_ch = snow_ch.mix(ndvi_ch, solar_radiation_ch, topographic_suite_ch)


    emit:
        gee_tif_ch = gee_tif_ch

}