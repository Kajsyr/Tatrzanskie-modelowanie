// Proces łączący CSV
process MERGE_CSV {
    conda "${params.project_envs}/pytatra.yml"
    publishDir "${params.output_dir}/features/", mode: 'copy'
    storeDir "${params.cache_dir}/features/", mode: 'copy'

    input:
        path csv_files

    output:
        path "merged_results.csv"

    script:
    """
    #!/usr/bin/env python
    import pandas as pd
    import os

    # csv_files w bashu to lista oddzielona spacjami. 
    # Python odbierze to jako string, więc robimy z tego listę:
    files = "${csv_files}".split()

    dfs = []
    for f in files:
        df = pd.read_csv(f)
        # Wybieramy kolumny (używając Twojej logiki iloc)
        selected = pd.concat([
            df.iloc[:, [0]],
            df.iloc[:, [9]],
            df.iloc[:, 21:24],
            df.iloc[:, [25]],
            df.iloc[:, 50:]
        ], axis=1)
        dfs.append(selected)

    final_df = pd.concat(dfs, ignore_index=True)
    final_df.to_csv("merged_results.csv", index=False)
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


workflow BIO_CLIM {
    main:

    //Dowload BIOCLIM
    bioclim_ch = DOWNLOAD_BIOCLIM_ALPS()
    //chelsa_ch = DOWNLOAD_CHELSA_ALPS()

    //Plot BioCLIM
    PLOT_BIOCLIM(bioclim_ch)
    
    emit:
        bioclim_raster = bioclim_ch
        
}