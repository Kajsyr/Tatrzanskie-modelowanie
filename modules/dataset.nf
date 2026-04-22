nextflow.enable.dsl=2


process EXTRACT_SPECIES {
    input:
    path xlsx

    output:
    path "species.txt"

    script:
    """
    #!/usr/bin/env Rscript
    library(readxl)
    df <- read_excel("$xlsx", sheet = "species")
    writeLines(df\$Species, "species.txt")
    """
}


process GET_GBIF_DATA {
    publishDir "${params.output_dir}/raw_data/"
    storeDir "${params.cache_dir}/raw_data/" 
    conda "${params.project_envs}/rtatra.yml"

    //GBIF ma limity jednoczesnych pobrań około 3 
    maxForks 3
    tag "${species_name}"
    errorStrategy 'retry'
    maxRetries 3
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
    write.table(df, filename, sep = ",", row.names = FALSE)
    """
}


process PROCESS_DEM {
    publishDir "${params.output_dir}/spatial", mode: 'copy'
    storeDir "${params.cache_dir}/dem_alps"
    conda "${params.project_envs}/rtatra.yml"

    input:
    path shp_dir // folder zawierający pliki .shp, .shx, .dbf

    output:
    path "dem_alps.tif", emit: dem_tiff
    path "shp_alps_simply.rds", emit: shp_rds
    path "contours_alps.rds", emit: contours_rds

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

    dem <- get_elev_raster(locations = shp_alps_simply, z = 8, clip = "locations") 
    writeRaster(rast(dem), "dem_alps.tif", overwrite = TRUE)
    
    # Tworzenie podkładu pod mapy
    dem <- rast("dem_alps.tif")
    # Generujemy warstwice
    contours <- terra::as.contour(dem, levels = seq(0, 4000, 250))
    contours_sf <- st_as_sf(contours)

    # KLUCZOWE: Wymuszamy układ współrzędnych taki sam jak w DEM
    st_crs(contours_sf) <- st_crs(dem)

    # Konwertujemy na stopnie (WGS84), bo w tym układzie robimy mapy
    contours_sf <- st_transform(contours_sf, 4326)

    saveRDS(contours_sf, "contours_alps.rds")

    """    
}






process FILTER {
    publishDir "${params.output_dir}", mode: 'copy'
    storeDir "${params.cache_dir}"
    conda "${params.project_envs}/rtatra.yml"
    tag "${raw_csv}"
    cpus '6'
    
    input:
    path raw_csv
    path shp_rds

    output:
    path "filtered/*.csv", emit: cleaned_data
    path "filtered/*.log",  emit: logs

    script:
    """
    mkdir filtered/
    filter_gbif.R \
        --input ${raw_csv} \
        --output filtered/${raw_csv.simpleName} \
        --shape ${shp_rds} 
    """
}

process PLOT {
    publishDir "${params.output_dir}/plots", mode: 'copy'
    conda "${params.project_envs}/rtatra.yml"    
    tag "${raw_csv}"

    input:
    path raw_csv
    path contours_rds

    output:
    path "map_*.svg"

    script:
    def spec_name = raw_csv.simpleName.replace('_gbif', '').replace('_', ' ')
    """
    plot_gbif.R \
        --input ${raw_csv.simpleName}.csv \
        --contours ${contours_rds} \
        --species "${spec_name}" \
        --output map_${raw_csv.simpleName}.svg
    """
}


workflow DATASET{
    main:
    // wczytanie listy gatunków z eksel
    ch_species = EXTRACT_SPECIES(params.species_file)
        .splitText()
        .map { it.trim() }

    //Pobieranie GBIF (równolegle dla każdego gatunku)
    download_gbif_ch = GET_GBIF_DATA(ch_species)    

    //odfultujemy puste pliki bo cześć odpada wcześniej
    gbif_ch = download_gbif_ch.filter{ file -> 
        file.readLines().size() > 1 
    }
    
    //DEM i uproszczony kształt Alp
    dem_ch = PROCESS_DEM(file(params.alp_shapefile).parent)

    // Filtrowanie i mapowanie wystąpień
    filtered_gbif_ch= FILTER(gbif_ch, dem_ch.shp_rds)

    //odfultujemy puste pliki bo cześć odpada wcześniej
    valid_gbif_ch = filtered_gbif_ch.cleaned_data
            .flatten()
            .filter { file -> file.countLines() > 1} 
            .unique()

    //Mapowanie
    PLOT(valid_gbif_ch, dem_ch.contours_rds)

    emit:
        species_ch = valid_gbif_ch
        dem_ch = dem_ch.dem_tiff

}


