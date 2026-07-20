# =============================================================================
# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 01: Limpieza y preparación de registros
# Autoras: Maria Judith Carmona


library(readxl)
library(dplyr)


# 1. Leer archivo original


raw <- read_excel("datos/registros/Copia_Puntos_georreferenciados.xlsx") %>%
  mutate(fila_excel = row_number())


# 2. Seleccionar variables de interés


registros <- raw %>%
  select(
    fila_excel,
    
    tax = Species,
    
    lat_corr = `Latitud_corregida (Y)`,
    lon_corr = `Longitud_corregida (X)`,
    
    lat_orig = Latitud,
    lon_orig = Longitud,
    
    elev_msnm = `Elevación (msnm)`,
    
    colector = Collector,
    herbario = Herbaria,
    localidad = Locality,
    
    pais = country,
    departamento = department,
    municipio = municipality,
    
    fecha = Date,
    habito = Habit,
    
    id = ID
  )


# 3. Excluir híbridos


registros <- registros %>%
  filter(
    !tax %in% c(
      "Serpocaulon x manizalense",
      "Serpocaulon x semipinnatifidum"
    )
  )


# 4. Limpiar coordenadas


registros <- registros %>%
  mutate(
    
    # eliminar espacios
    lat_corr = na_if(trimws(lat_corr), ""),
    lon_corr = na_if(trimws(lon_corr), ""),
    
    lat_orig = na_if(trimws(lat_orig), ""),
    lon_orig = na_if(trimws(lon_orig), ""),
    
    # convertir coma decimal a punto
    lat_corr = gsub(",", ".", lat_corr),
    lon_corr = gsub(",", ".", lon_corr),
    
    lat_orig = gsub(",", ".", lat_orig),
    lon_orig = gsub(",", ".", lon_orig),
    
    # convertir a numérico
    lat_corr = as.numeric(lat_corr),
    lon_corr = as.numeric(lon_corr),
    
    lat_orig = as.numeric(lat_orig),
    lon_orig = as.numeric(lon_orig),
    
    # usar coordenada corregida; si no existe usar la original
    ddlat = coalesce(lat_corr, lat_orig),
    ddlon = coalesce(lon_corr, lon_orig)
  )


# 5. Corregir longitudes positivas para Colombia


registros <- registros %>%
  mutate(
    ddlon = if_else(
      pais == "Colombia" &
        !is.na(ddlon) &
        ddlon > 0,
      -ddlon,
      ddlon
    )
  )


# 6. Conservar variables finales


registros <- registros %>%
  select(
    fila_excel,
    tax,
    ddlat,
    ddlon,
    elev_msnm,
    colector,
    herbario,
    localidad,
    pais,
    departamento,
    municipio,
    fecha,
    habito,
    id
  )


# 7. Base completa (todos los registros)


registros_totales <- registros


# 8. Base para ConR (solo registros georreferenciados)


registros_limpios <- registros_totales %>%
  filter(
    !is.na(ddlat),
    !is.na(ddlon)
  )


# 9. Control de calidad


cat("\n================ CONTROL DE CALIDAD ================\n")

cat("Registros originales:                     ",
    nrow(raw), "\n")

cat("Registros sin híbridos:                   ",
    nrow(registros_totales), "\n")

cat("Registros recuperados desde Lat/Long:     ",
    sum(
      is.na(raw$`Latitud_corregida (Y)`) &
        !is.na(raw$Latitud)
    ),
    "\n"
)

cat("Registros sin coordenadas:                ",
    sum(
      is.na(registros_totales$ddlat) |
        is.na(registros_totales$ddlon)
    ),
    "\n"
)

cat("Registros finales para ConR:              ",
    nrow(registros_limpios), "\n")

cat("====================================================\n\n")


# 10. Resumen por especie


resumen <- registros_limpios %>%
  count(tax, name = "n_registros") %>%
  arrange(tax)

print(resumen)


# 11. Guardar archivos


write.csv(
  registros_totales,
  "datos/registros/registros_totales.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  registros_limpios,
  "datos/registros/registros_limpios.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\nArchivos exportados correctamente:\n")
cat(" - datos/registros/registros_totales.csv\n")
cat(" - datos/registros/registros_limpios.csv\n")