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

    bioclim_names <- c(
        bio1="Annual Mean Temperature (°C)",
        bio2="Mean Diurnal Range",
        bio3="Isothermality",
        bio4="Temperature Seasonality",
        bio5="Max Temp Warmest Month",
        bio6="Min Temp Coldest Month",
        bio7="Temp Annual Range",
        bio8="Mean Temp Wettest Quarter",
        bio9="Mean Temp Driest Quarter",
        bio10="Mean Temp Warmest Quarter",
        bio11="Mean Temp Coldest Quarter",
        bio12="Annual Precipitation (mm)",
        bio13="Precipitation Wettest Month",
        bio14="Precipitation Driest Month",
        bio15="Precipitation Seasonality",
        bio16="Precipitation Wettest Quarter",
        bio17="Precipitation Driest Quarter",
        bio18="Precipitation Warmest Quarter",
        bio19="Precipitation Coldest Quarter"
    )


    # 1. Definicja Twojego obszaru (Bounding Box na podstawie POLYGON)
    # xmin: 4, xmax: 18, ymin: 43.5, ymax: 48.5
    roi_ext <- ext(4, 18, 43.5, 48.5)

    # 2. Pobieranie tylko kafla zawierającego Alpy (30x30 stopni)
    # res = 0.5 (ok. 1km) jest optymalne dla obszarów górskich (Alpy)
    # Spróbujmy tylko z wokół jednego punktu
    cat("Pobieranie kafla WorldClim (res=0.5)...\\n")
    bio_tile <- worldclim_tile(var = "bio", res = 0.5, lon = 11, lat = 46, path = ".")

    # 3. Przycięcie do Twojego obszaru
    cat("Przycinanie do zadanego poligonu...\\n")
    bio_alps <- crop(bio_tile, roi_ext)

    #Tworzenie folderu
    if (!dir.exists("bioclim")) dir.create("bioclim")

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

process PLOT_BIOCLIM {

    publishDir "${params.output_dir}/plots/", mode: 'copy'
    conda "${params.project_envs}/rtatra.yml"

    input:
    path bioclim_dir

    output:
    path "alps*.png"

    script:
    """
    #!/usr/bin/env Rscript
    library(terra)
    library(ggplot2)


    bioclim_names <- c(
        bio1="Annual Mean Temperature (°C)",
        bio2="Mean Diurnal Range",
        bio3="Isothermality",
        bio4="Temperature Seasonality",
        bio5="Max Temp Warmest Month",
        bio6="Min Temp Coldest Month",
        bio7="Temp Annual Range",
        bio8="Mean Temp Wettest Quarter",
        bio9="Mean Temp Driest Quarter",
        bio10="Mean Temp Warmest Quarter",
        bio11="Mean Temp Coldest Quarter",
        bio12="Annual Precipitation (mm)",
        bio13="Precipitation Wettest Month",
        bio14="Precipitation Driest Month",
        bio15="Precipitation Seasonality",
        bio16="Precipitation Wettest Quarter",
        bio17="Precipitation Driest Quarter",
        bio18="Precipitation Warmest Quarter",
        bio19="Precipitation Coldest Quarter"
    )

    # Wczytaj warstwy
    files <- list.files("${bioclim_dir}", pattern="tif\$", full.names=TRUE)

    #iterowanie po warstwach i zapisywanie do pliku
    for(f in files){
        r <- rast(f)
        name <- tools::file_path_sans_ext(basename(f))

        #szukamy prawdziwej nazwy
        var_match <- regmatches(name, regexpr("bio_?\\\\d+", name))
        var <- gsub("_", "", var_match)
        label <- bioclim_names[[var]]

        df <- as.data.frame(r, xy=TRUE, na.rm=TRUE)
        colnames(df)[3] <- "value"

        p <- ggplot(df, aes(x=x, y=y, fill=value)) +
            geom_raster() +
            coord_equal() +
            scale_fill_viridis_c() +
            labs(title=name) +
            theme_minimal()

        ggsave(paste0(name, ".png"), plot=p, width=8, height=6)
    }

    """
}

process GEE_DOWNSCALING_ALPS {
    conda "conda-forge::geemap conda-forge::earthengine-api"
    storeDir "${params.cache_dir}/gee_highres"

    output:
    path "niche_30m/*.tif", emit: highres_tiffs

    script:
    """
    #!/usr/bin/env python3
    import ee
    import geemap
    import os

    # Inicjalizacja (wymaga skonfigurowanego serwisu lub tokena w środowisku)
    ee.Initialize()

    roi = ee.Geometry.Rectangle([4.0, 43.5, 18.0, 48.5]) # Rozszerzony zasięg Alp
    
    # 1. Dane źródłowe
    bioclim = ee.Image("WORLDCLIM/V1/BIO").select('bio01').divide(10)
    dem_high = ee.Image("NASA/NASADEM_HGT/001").select('elevation')
    
    # 2. Downscaling (Lapse Rate Method)
    dem_low = dem_high.resample('bilinear').reproject(crs=bioclim.crs(), scale=1000)
    lapse_rate = -0.0065
    bioclim_30m = bioclim.add(dem_high.subtract(dem_low).multiply(lapse_rate)).rename('temp_30m')

    # 3. Topografia
    slope = ee.Terrain.slope(dem_high).rename('slope')
    aspect = ee.Terrain.aspect(dem_high).rename('aspect')
    
    stack = ee.Image.cat([bioclim_30m, dem_high, slope, aspect]).clip(roi)

    # 4. Eksport lokalny (przez geemap - dla małych/średnich obszarów)
    # Dla bardzo dużych obszarów lepiej użyć Export.image.toDrive i pobierać przez API
    if not os.path.exists('niche_30m'): os.makedirs('niche_30m')
    
    geemap.ee_export_image(
        stack, 
        filename="niche_30m/alps_highres_stack.tif", 
        scale=30, 
        region=roi, 
        file_per_band=True
    )
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
    storeDir "${params.cache_dir}/features"
    
    input:
    path gbif_species
    path dem_tiff
    path bioclim_dir
    //path gge_dir
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
     --bioclim ${bioclim_dir}/
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

    //Dowload BIOCLIM
    bioclim_ch = DOWNLOAD_BIOCLIM_ALPS()
    //chelsa_ch = DOWNLOAD_CHELSA_ALPS()

    //PLOT
    PLOT_BIOCLIM(bioclim_ch)

    //dowload data using GEE
    //gee_ch = GEE_DOWNSCALING_ALPS()

    //extract envirmnetal data
    env_data = EXTRACT_ENVIRONMENTAL_DATA(
        dataset_ch,
        dem_tiff,
        bioclim_ch
        //gee_ch
        //chelsa_ch.collect(),
        
    )

    //merge_env_data=MERGE_CSV(env_data)
    emit:
        env_data
        //merge_env_data
        
}