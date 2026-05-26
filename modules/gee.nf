process DOWNLOAD_SNOW_COVER {
    storeDir "${params.cache_gee}/snow"
    conda "${params.project_envs}/geedim.yml"

    output:
    path "*.tif"

    script:
    """
    geedim download \
            -i projects/lofty-complex-481222-k8/assets/snow_duration_alpy_2023 \
            --scale 20 \
            --max-requests ${params.max_requests} \
            --max-tile-size ${params.max_tile_size}
    mv *.tif snow_duration.tif
    """
}

process DOWNLOAD_ECO_FEATURES {
    storeDir "${params.cache_gee}/eco_features"
    conda "${params.project_envs}/geedim.yml"

    output:
    path "*.tif"

    script:
    """
    geedim download \
        -i projects/lofty-complex-481222-k8/assets/eco_features_alpy_2023 \
        --scale 30 \
        --max-requests ${params.max_requests} \
        --max-tile-size ${params.max_tile_size}

    mv *.tif eco_features.tif
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
        --scale 30 \
        --max-requests ${params.max_requests} \
        --max-tile-size ${params.max_tile_size}

    mv *.tif solar_radiation.tif
    """
}

process DOWNLOAD_TOPOGRAPHIC_SUITE {
    storeDir "${params.cache_gee}/topographic_suite"
    conda "${params.project_envs}/geedim.yml"
    tag "${band_name}"
    debug true

    input:
    val band_name

    output:
    path "${band_name}_alpy.tif"

    script:
    """
    geedim download \
        -i projects/lofty-complex-481222-k8/assets/topography_suite_alpy \
        --scale 30 \
        --band-name ${band_name} \
        --max-requests ${params.max_requests} \
        --max-tile-size ${params.max_tile_size}

    mv *.tif ${band_name}_alpy.tif
    """
}


process DOWNLOAD_SOIL {
    storeDir "${params.cache_gee}/soil"
    conda "${params.project_envs}/geedim.yml"
    tag "${soil_band}"
    debug true

    input:
    val soil_band

    output:
    path "soil_${soil_band}_alpy.tif"

    script:
    """
    geedim download \
        -i projects/lofty-complex-481222-k8/assets/soil_litho_suite_alpy \
        --scale 30 \
        --crs EPSG:32632 \
        --band-name ${soil_band} \
        --max-requests ${params.max_requests} \
        --max-tile-size ${params.max_tile_size}

    mv *.tif soil_${soil_band}_alpy.tif
    """
}


workflow GEE {
    main:
        channel
        .fromList(['DEM', 'Slope', 'Aspect', 'TRI', 'Curvature', 'TWI'])
        .set { topographic_bands_ch }

        channel
        .fromList(["Soil_pH", "Soil_Organic_Carbon", "Soil_Sand", "Soil_Clay"])
        .set { soil_bands_ch }

        snow_ch = DOWNLOAD_SNOW_COVER()
        eco_features_ch = DOWNLOAD_ECO_FEATURES()
        solar_radiation_ch = DOWNLOAD_SOLAR_RADIATION()
        topographic_suite_ch = DOWNLOAD_TOPOGRAPHIC_SUITE(topographic_bands_ch)
        soil_ch = DOWNLOAD_SOIL(soil_bands_ch)

        //take all tif to one channel
        gee_tif_ch = snow_ch.mix(eco_features_ch, 
                                solar_radiation_ch, 
                                topographic_suite_ch,
                                soil_ch)


    emit:
        //gee_tif_ch = topographic_suite_ch
        gee_tif_ch = gee_tif_ch

}