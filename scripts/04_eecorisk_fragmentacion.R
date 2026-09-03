# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 04: Fragmentación severa, huella humana y disminución de hábitat (eecorisk)
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Referencia: GEPC - Grupo de Especialistas en Plantas de Colombia
#
# Este script implementa la metodología eecorisk del GEPC para evaluar las
# condiciones del hábitat bajo el Criterio B de la IUCN (subcriterios b(iii) y b(iv)):
#   - Fragmentación severa del hábitat disponible (AOH)
#   - Disminución continua de la calidad/cantidad de hábitat (huella humana)
#   - Disminución continua de subpoblaciones (deforestación histórica BNB)
#
# RASTERS REQUERIDOS - descargar manualmente y ubicar en datos/capas/:
#
#   coberturas_tierra/   → Capas de Bosque/No Bosque (BNB) del IDEAM
#     Qué es: producto de teledetección del IDEAM que clasifica cada píxel como
#     Bosque (1) o No Bosque (0) para Colombia continental. Se genera a partir de
#     imágenes Landsat y refleja la cobertura boscosa al año indicado.
#     Se usan tres fechas (1990, 2000, 2024) para detectar cambios históricos.
#     Resolución original: ~30 m. Se remuestrea a 300 m para este análisis (ver abajo).
#     Descargar desde: https://experience.arcgis.com/experience/568ddab184334f6b81a04d2fe9aac262
#     Buscar: "Mapa de Bosque No Bosque" para los años 1990, 2000 y 2024
#     Guardar como: BNB_1990_original.tif, BNB_2000_original.tif, BNB_2024_original.tif
#     El script los remuestrea automáticamente a 300m la primera vez.
#
#   coberturas_tierra/iheh_col.tif   → Índice de Huella Espacial Humana (IHEH) Colombia
#     Qué es: raster que cuantifica la intensidad de la presión humana acumulada
#     sobre el territorio colombiano, integrando variables como: densidad de población,
#     infraestructura vial, uso del suelo agropecuario, luces nocturnas e índice de
#     fragmentación del paisaje. Los valores van de 0 (sin huella) a 100 (máxima huella).
#     Fuente: Correa Ayram et al. (2020) "Spatiotemporal evaluation of the human
#     footprint in Colombia: Four decades of anthropic impact in highly diverse ecosystems."
#     Ecological Indicators, 117, 106630. https://doi.org/10.1016/j.ecolind.2020.106630
#     Disponible en el repositorio del proyecto (datos/capas/coberturas_tierra/iheh_col.tif).
#
# PARÁMETROS PARA Serpocaulon (helechos epífitos de bosque húmedo):
#
#   umbral = 150 km²
#     Tamaño mínimo de parche para que NO sea considerado "pequeño".
#     Para epífitas asociadas a bosque continuo, un parche < 150 km² es
#     insuficiente para mantener poblaciones viables a largo plazo.
#     El GEPC define este valor según el grupo funcional de la especie.
#
#   disper = 50 km (mayoría de especies) / 0.06 km = 60 m (10 especies)
#     Distancia máxima de dispersión efectiva para recolonización de parches.
#     Para la mayoría de especies se usa 50 km (valor GEPC para pteridófitas).
#     Para 10 especies con dispersión limitada (experimentos del experto):
#       S. antioquianum, S. attenuatum, S. biauriculatum, S. concolorum,
#       S. loriceum, S. patentissimum, S. polystichum, S. richardii,
#       S. tayronae, S. wagnerii → disper = 0.06 km (60 m).
#     Con 60 m de dispersal, cualquier par de parches separados por ≥ 300 m
#     (= 1 celda del raster) se considera aislado; el FS_score refleja entonces
#     principalmente si los parches son pequeños (< 150 km²).
#
#   umbral_HH = 40%
#     Porcentaje de huella humana promedio en el AOH a partir del cual se declara
#     disminución continua de hábitat (cod_dism_habitat = YES).
#     El índice IHEH va de 0 (sin intervención) a 100 (completamente transformado).
#     Un promedio ≥ 40% indica que una fracción importante del hábitat disponible
#     está sometida a presión humana significativa y continua.
#     El GEPC usa 40% como umbral estándar para plantas vasculares colombianas.

# Librerías ----
library(raster)
library(sf)
library(geosphere)
library(sp)
library(elevatr)   # descarga DEM automáticamente
library(dplyr)
library(writexl)

# Datos ----
# Cargar registros y agregar parámetros de especie
# umbral y disper se añaden como columnas porque AHO_fast y sfrag los leen
# directamente desde la tabla de puntos (columnas 4 y 5 respectivamente).

# Especies con dispersión efectiva limitada a 60 m (experimentos del experto)
especies_disper_bajo <- c(
  "Serpocaulon antioquianum", "Serpocaulon attenuatum",
  "Serpocaulon biauriculatum", "Serpocaulon concolorum",
  "Serpocaulon loriceum",      "Serpocaulon patentissimum",
  "Serpocaulon polystichum",   "Serpocaulon richardii",
  "Serpocaulon tayronae",      "Serpocaulon wagnerii"
)

registros <- read.csv(
  "datos/registros/registros_limpios.csv",
  encoding = "UTF-8") %>%
  filter(!is.na(ddlat), !is.na(ddlon)) %>%
  mutate(
    elev_msnm = as.numeric(elev_msnm),
    umbral = 150,
    disper = ifelse(tax %in% especies_disper_bajo, 0.06, 50)
  )

# Diagnóstico: registros sin elevación en el CSV original (se completarán con DEM)
registros_raw <- read.csv("datos/registros/registros_limpios.csv",
  encoding = "UTF-8")

registros_raw %>%
  filter(is.na(suppressWarnings(as.numeric(elev_msnm)))) %>%
  dplyr::select(tax, elev_msnm)

# DEM ----
# DEM de elevación
# Se usa para filtrar celdas del BNB fuera del rango altitudinal de la especie
# dentro de la función AHO_fast. Resolución z=6 (~1 km) es suficiente para este
# filtro grueso; la elevación precisa por registro viene del campo elev_msnm del CSV.
# z=6 descarga ~80 MB; z=9 (~30 m) descargaría ~1.2 GB y no mejora los resultados.
ruta_elv <- "datos/capas/elevacion/dem_colombia.tif"

if (!file.exists(ruta_elv)) {
  pts_sf <- st_as_sf(registros, coords = c("ddlon", "ddlat"), crs = 4326)
  elv <- get_elev_raster(pts_sf, z = 6, src = "aws")
  writeRaster(elv, ruta_elv, overwrite = TRUE)
} else {
  elv <- raster(ruta_elv)
}

# Completar elevación faltante desde el DEM

# Coordenadas de todos los registros
coords <- registros[, c("ddlon", "ddlat")]

# Extraer elevación del DEM para cada registro
elev_dem <- raster::extract(elv, coords)

# Completar únicamente los NA del CSV
registros$elev_msnm <- ifelse(
  is.na(registros$elev_msnm),
  elev_dem,
  registros$elev_msnm
)

message(sum(is.na(registros$elev_msnm)), " registros continúan sin elevación.")

# BNB remuestreo 30m → 300m ----
# Remuestreo BNB de 30 m a 300 m
# El BNB original del IDEAM tiene ~30 m de resolución (píxeles Landsat).
# Se remuestrea a 300 m (factor 10) por tres razones:
#   1. Eficiencia: clump() y las operaciones de parches sobre un raster de 30 m
#      para toda Colombia tardarían horas o días; a 300 m tarda minutos.
#   2. Escala apropiada: las subpoblaciones de Serpocaulon se definen a escala de
#      kilómetros (buffer 5 km en script 02); detectar parches a 30 m introduce
#      ruido que no corresponde a unidades ecológicas reales.
#   3. Compatibilidad con IHEH: el raster de huella humana tiene resolución similar,
#      facilitando la integración de capas.
# fun = min en aggregate: una celda de 300 m se clasifica como Bosque (1) solo si
#   TODOS (o la mayoría) sus píxeles de 30 m son bosque. Usar min es conservador:
#   una celda queda como bosque únicamente si no hay píxeles de no-bosque en ella,
#   lo que evita sobreestimar el hábitat disponible.
# SB10 = BNB 2024 a 300 m → grilla de referencia para el análisis AHO
# SB1000 = BNB 2024 a 600 m (aggregate factor 2 sobre SB10) → grilla más gruesa
#   usada solo para encontrar el extent de la especie eficientemente en el loop AHO.
# Los BNB de 1990 y 2000 se remuestrean al mismo grid que SB10 (mismo origen y
# resolución) para garantizar alineación exacta celda a celda en el stack BNBstk.
if (!file.exists("datos/capas/coberturas_tierra/BNB_300_2024.tif")) {
  r2024 <- raster("datos/capas/coberturas_tierra/BNB_2024_original.tif")
  r2024 <- aggregate(r2024, fact = 10, fun = min)
  writeRaster(r2024, "datos/capas/coberturas_tierra/BNB_300_2024.tif",
              format = "GTiff", overwrite = TRUE)
} else {
  r2024 <- raster("datos/capas/coberturas_tierra/BNB_300_2024.tif")
}

SB10   <- raster("datos/capas/coberturas_tierra/BNB_300_2024.tif")
SB1000 <- aggregate(r2024, fact = 2, fun = max)

for (año in c("1990", "2000")) {
  ruta_300 <- paste0("datos/capas/coberturas_tierra/BNB_300_", año, ".tif")
  if (!file.exists(ruta_300)) {
    r <- raster(paste0("datos/capas/coberturas_tierra/BNB_", año, "_original.tif"))
    r <- aggregate(r, fact = 10, fun = min)
    r <- resample(r, SB10)
    writeRaster(r, ruta_300, format = "GTiff", overwrite = TRUE)
  }
}

# Stack histórico BNB (1990 → 2000 → 2024)
# Las tres capas apiladas permiten comparar la cobertura boscosa en tres momentos,
# detectando si parches con registros de la especie perdieron bosque con el tiempo.
BNBstk <- stack(
  raster("datos/capas/coberturas_tierra/BNB_300_1990.tif"),
  raster("datos/capas/coberturas_tierra/BNB_300_2000.tif"),
  SB10)

# IHEH ----
# Índice de Huella Espacial Humana (IHEH)
# Se reproyecta al mismo grid que SB10 (misma resolución, extent y CRS) para poder
# enmascararlo con el AOH de cada especie y calcular el promedio de HH en el hábitat.
# projectRaster() interpola bilinealmente los valores continuos del IHEH.
hh <- raster("datos/capas/coberturas_tierra/iheh_col.tif")
hh <- projectRaster(hh, SB10)

# Tabla base eecorisk ----
# Elevación mín/máx por especie desde el CSV de registros (campo elev_msnm)
elev_min <- registros %>% filter(!is.na(elev_msnm)) %>%
  group_by(tax) %>% summarise(MinElv = min(elev_msnm), .groups = "drop")
elev_max <- registros %>% filter(!is.na(elev_msnm)) %>%
  group_by(tax) %>% summarise(MaxElv = max(elev_msnm), .groups = "drop")

# Tabla base para eecorisk

corsp <- registros %>%
  dplyr::select(
    especie = tax,
    Latitud = ddlat,
    Longitud = ddlon,
    umbral,
    disper
  ) %>%
  left_join(elev_min, by = c("especie" = "tax")) %>%
  left_join(elev_max, by = c("especie" = "tax"))

# Preparar objetos por especie

ne <- unique(corsp$especie)

dsp <- lapply(ne, function(sp) {
  corsp[corsp$especie == sp, ]
})

csp <- lapply(dsp, na.omit)
csp1 <- csp

for(i in seq_along(ne)){
  
  if(nrow(csp[[i]]) == 0) next
  
  coordinates(csp1[[i]]) <- c("Longitud","Latitud")
  
}

# Parámetros que no cambian entre especies

resx <- xres(SB1000)
resy <- yres(SB1000)

# Objeto donde se almacenará el AOH

AOOok <- vector("list", length(ne))

# Función AHO_fast ----

# AHO <- function(model,
#                 puntos,
#                 xy,
#                 bufferSize = 0.054,
#                 bufferPoints = TRUE,
#                 elv){
#   
#   if(!inherits(model,"RasterLayer"))
#     stop("model debe ser RasterLayer")
#   
#   elvg <- resample(
#     crop(elv, model),
#     model
#   )
#   
#   elvgv <- getValues(elvg)
#   
#   model[which(elvgv < as.numeric(unique(puntos[,6])) - 300)] <- 0
#   model[which(elvgv > as.numeric(unique(puntos[,7])) + 300)] <- 0
#   
#   groups <- clump(
#     model,
#     directions = 4
#   )
#   
#   pts_sf <- st_as_sf(
#     puntos,
#     coords = names(puntos)[xy],
#     crs = 4326
#   )
#   
#   if(bufferPoints){
#     
#     # Proyectar temporalmente a metros
#     pts_m <- st_transform(pts_sf, 3857)
#     
#     pbf <- st_buffer(
#       pts_m,
#       dist = bufferSize
#     )
#     
#     # Regresar a WGS84
#     pbf <- st_transform(pbf, 4326)
#     
#   }else{
#     
#     pbf <- pts_sf
#     
#   }
#   
#   pex <- extract(
#     groups,
#     as(pbf,"Spatial")
#   )
#   
#   parches <- na.omit(
#     unique(
#       pex[[1]]
#     )
#   )
#   
#   pexok <- match(
#     getValues(groups),
#     parches
#   )
#   
#   groups[which(pexok > 0)] <- 1
#   groups[which(is.na(pexok))] <- 0
#   
#   return(groups)
#   
# }


# AHO_fast - Área de Hábitat Ocupado (Area of Habitat Occupied)
# Determina qué parches de bosque (celdas BNB=1) forman parte del hábitat
# real de la especie, considerando:
#   1. Rango altitudinal: excluye celdas fuera del rango (mínElv-300, máxElv+300 m).
#      El margen de 300 m acomoda la variación microclimática y de muestreo.
#   2. Buffer de presencia: identifica parches de bosque que están dentro de
#      `bufferSize` metros de algún registro. Solo esos parches se incluyen en el AOH.
#
# bufferSize = 6000 m (6 km):
#   Equivale al 0.054° que usaba la función original en grados decimales
#   (0.054° × 111 km/° ≈ 6 km). Se usa metros porque st_buffer() en un CRS
#   proyectado (EPSG:3857) es exacto; en grados la distancia varía con la latitud.
#   El buffer captura parches de bosque adyacentes a los puntos de registro,
#   reconociendo que el espécimen pudo haber sido colectado en el borde del parche.
#
# La versión "fast" reemplaza raster::extract() por rasterize() para identificar
# los IDs de parches dentro del buffer, lo que es ~10x más rápido en rasters grandes.
AHO_fast <- function(model,
                     puntos,
                     xy,
                     bufferSize = 6000,
                     bufferPoints = TRUE,
                     elv){
  
  if(!inherits(model, "RasterLayer"))
    stop("model debe ser RasterLayer")
  
  ## Elevación
  elvg <- resample(
    crop(elv, model),
    model
  )
  
  elvgv <- getValues(elvg)
  
  model[elvgv < unique(puntos[,6]) - 300] <- 0
  model[elvgv > unique(puntos[,7]) + 300] <- 0
  
  ## Parches
  groups <- clump(
    model,
    directions = 4
  )
  
  ## Puntos
  pts_sf <- st_as_sf(
    puntos,
    coords = names(puntos)[xy],
    crs = 4326
  )
  
  if(bufferPoints){
    
    pts_sf <- st_transform(pts_sf, 3857)
    
    pbf <- st_buffer(
      pts_sf,
      dist = bufferSize
    )
    
    pbf <- st_transform(pbf, 4326)
    
  }else{
    
    pbf <- pts_sf
    
  }
  
  ## Rasterizar buffer (mucho más rápido que extract)
  pbf_sp <- as(pbf, "Spatial")
  
  mask_buf <- rasterize(
    pbf_sp,
    groups,
    field = 1,
    background = NA
  )
  
  ## IDs de parches dentro del buffer
  ids <- unique(
    getValues(
      mask(
        groups,
        mask_buf
      )
    )
  )
  
  ids <- ids[!is.na(ids)]
  
  gval <- getValues(groups)
  
  groups[] <- ifelse(
    gval %in% ids,
    1,
    0
  )
  
  return(groups)
  
}

# Función sfrag ----
# sfrag - Fragmentación Severa (Severe Fragmentation)
# Evalúa si el hábitat de la especie está severamente fragmentado según el Criterio B IUCN.
# Para cada parche de bosque en el AOH calcula:
#   - Área (km²): número de celdas × área por celda (92,106 m² a 300 m de resolución)
#   - Distancia al parche más cercano (m): usando centroides y distancias geodésicas
#   - Small: TRUE si área < umbral (150 km² para Serpocaulon)
#   - Isolated: TRUE si distancia al vecino - radio del parche > disper (50 km)
#     La resta del radio (sqrt(Area/π)) corrige el hecho de que la distancia
#     entre centroides sobreestima la distancia entre bordes en parches grandes.
#
# FS_score = % de parches que son simultáneamente pequeños Y aislados.
# Fragmentación severa (FS = TRUE) si FS_score > 50%:
#   Más de la mitad de los parches donde vive la especie son pequeños y están
#   tan aislados que la recolonización tras una extinción local es improbable.
#
# bufferSize = 20 (no se usa en la versión actual con bufferPoints = FALSE en sfrag,
#   pero se conserva por compatibilidad con la función original del GEPC).
sfrag <- function(BNB, puntos, xy = c(2, 3), bufferSize = 20, bufferPoints = TRUE) {
  if (!inherits(BNB, "RasterLayer")) stop("BNB debe ser RasterLayer")
  groups <- clump(BNB, directions = 4)
  dp     <- na.omit(as.data.frame(groups, xy = TRUE, centroids = TRUE))
  clon   <- (tapply(dp[,1], dp[,3], min) + tapply(dp[,1], dp[,3], max)) / 2
  clat   <- (tapply(dp[,2], dp[,3], min) + tapply(dp[,2], dp[,3], max)) / 2
  npix   <- tapply(dp[,1], dp[,3], length)
  Area   <- (92106 * npix) / 1e6
  corc1  <- as.data.frame(t(rbind(clon, clat))); coordinates(corc1) <- c("clon", "clat")
  if (length(Area) < 2) {
    return(list(data.frame(
        "Area km^2" = Area,
        "Dist_PMC m" = NA,
        Isolated = "Solo un parche",
        Small = Area < unique(puntos[,4])),NA,NA))
  }
  dis      <- apply(as.data.frame(distm(corc1)), 2, as.numeric); dis[dis == 0] <- NA
  minall   <- apply(dis, 2, function(x) min(x, na.rm = TRUE))
  disok    <- as.numeric(minall) - sqrt(Area / pi)
  Isolated <- disok > (unique(puntos[, 5]) * 1000)
  Small    <- Area < unique(puntos[, 4])
  FS_score <- 100 * length(na.omit(match(which(Small), which(Isolated)))) / length(Small)
  FS       <- FS_score > 50
  Tall     <- cbind("Area km^2" = Area, "Dist_PMC m" = disok,
                    "Isolated" = as.character(Isolated), "Small" = as.character(Small))
  return(list(Tall, FS_score, FS))
}

# Loop AHO por especie ----
pb <- txtProgressBar(min = 0,
                     max = length(ne),
                     style = 3)

for(i in seq_along(ne)){
  
  setTxtProgressBar(pb, i)
  
  message(ne[i])
  
  if(nrow(csp[[i]]) == 0) next
  
  cex <- raster::extract(
    SB1000,
    csp1[[i]],
    cell = TRUE
  )
  
  ext_sp <- extentFromCells(
    SB1000,
    unique(cex[,1])
  )
  
  ext_sp@xmin <- ext_sp@xmin - 2 * resx
  ext_sp@xmax <- ext_sp@xmax + 2 * resx
  ext_sp@ymin <- ext_sp@ymin - 2 * resy
  ext_sp@ymax <- ext_sp@ymax + 2 * resy
  
  SB10E <- crop(
    SB10,
    ext_sp
  )
  
  SB10E[SB10E == 2] <- 0
  SB10E[SB10E > 1]  <- 1
  SB10E[SB10E == 3] <- 0
  
  AOOok[[i]] <- AHO_fast(
    SB10E,
    csp[[i]],
    xy = c(3,2),
    bufferSize = 6000,
    bufferPoints = TRUE,
    elv = elv
  )
  
}
close(pb)

# Loop sfrag por especie ----
sg <- vector("list", length(ne))
pb2 <- txtProgressBar(min = 0, max = length(ne), style = 3)
for (i in seq_along(ne)) {
  setTxtProgressBar(pb2, i)
  message(" ", ne[i])
  if (is.null(AOOok[[i]])) next
  if (length(unique(AOOok[[i]])) == 1 && unique(AOOok[[i]]) == 0) next
  sg[[i]] <- sfrag(AOOok[[i]], dsp[[i]], xy = c(3, 2), bufferSize = 20, bufferPoints = TRUE)
}
close(pb2)

# Huella humana en AOH ----
pct_hh <- sapply(seq_along(ne), function(i) {
  
  if (is.null(AOOok[[i]]))
    return(NA)
  
  aoo_mask <- AOOok[[i]]
  
  aoo_mask[aoo_mask == 0] <- NA
  
  hh_crop <- crop(hh, aoo_mask)
  hh_mask <- mask(hh_crop, aoo_mask)
  
  round(mean(getValues(hh_mask), na.rm = TRUE), 0)
  
})

# Subpoblaciones desaparecidas (discon) ----
# discon - Discontinuidad / Subpoblaciones Desaparecidas
# Detecta si alguna subpoblación de la especie desapareció entre el pasado y el presente,
# comparando el stack histórico BNB (1990, 2000) contra la capa actual (2024).
# Lógica por especie:
#   Para cada capa histórica (1990, 2000):
#     1. Recorta BNBstk al extent de la especie + buffer de 2 celdas (eficiencia)
#     2. Identifica parches de bosque (clump) en esa capa histórica
#     3. Extrae los IDs de parches donde había registros de la especie
#     4. Verifica si esas celdas tienen bosque en 2024
#     5. Si un parche histórico con registros NO tiene bosque en 2024 → perdida = TRUE
# cod_dism_subpob = YES indica evidencia de pérdida de hábitat que implica
# posible desaparición de subpoblaciones (subcriterio B2b(iv) de la IUCN).
#
# Nota: clump() sobre el raster completo de Colombia tarda horas; sobre el extent
# recortado de una especie tarda segundos — por eso se recorta primero.
# Detectar subpoblaciones desaparecidas por especie
# Se recorta BNBstk a la extensión de cada especie antes de clumpear
# (clump sobre Colombia completa es inviablemente lento)
pb3 <- txtProgressBar(min = 0, max = length(ne), style = 3)

subpob_perdida <- sapply(seq_along(ne), function(i) {

  setTxtProgressBar(pb3, i)

  if (nrow(csp[[i]]) == 0) return(NA)

  pts    <- csp[[i]]
  coords <- as.matrix(pts[, c(3, 2)])   # Longitud, Latitud

  # Extensión con buffer de 2 celdas alrededor de los puntos
  cex    <- raster::extract(SB10, coords, cell = TRUE)
  ext_sp <- extentFromCells(SB10, unique(cex[, 1]))
  ext_sp@xmin <- ext_sp@xmin - 2 * xres(SB10)
  ext_sp@xmax <- ext_sp@xmax + 2 * xres(SB10)
  ext_sp@ymin <- ext_sp@ymin - 2 * yres(SB10)
  ext_sp@ymax <- ext_sp@ymax + 2 * yres(SB10)

  stk_sp      <- crop(BNBstk, ext_sp)
  n_layers    <- nlayers(stk_sp)
  current_v   <- getValues(stk_sp[[n_layers]])

  for (j in 1:(n_layers - 1)) {

    cl       <- clump(stk_sp[[j]], directions = 4)
    vals_cl  <- getValues(cl)
    vals_pts <- raster::extract(cl, coords)
    parches  <- unique(na.omit(vals_pts))

    if (length(parches) == 0) next

    for (p in parches) {
      idx <- which(vals_cl == p)
      if (all(current_v[idx] == 0 | is.na(current_v[idx])))
        return(TRUE)
    }
  }

  FALSE
})

close(pb3)

# Tabla de resultados ----
umbral_HH <- 40   # % para declarar disminución continua de hábitat

Tablafrag <- data.frame(
  tax               = ne,
  FS_score          = sapply(sg, function(x) if (length(x) == 0) NA else round(as.numeric(x[[2]]), 0)),
  Frag_severa       = sapply(sg, function(x) if (length(x) == 0) NA else as.logical(x[[3]])),
  pct_HH            = pct_hh,
  subpob_perdida    = subpob_perdida
) %>%
  mutate(
    cod_fragmentacion      = case_when(is.na(Frag_severa)    ~ "Unknown",
                                       Frag_severa            ~ "YES",
                                       TRUE                   ~ "NO"),
    cod_dism_habitat       = case_when(is.na(pct_HH)         ~ "Unknown",
                                       pct_HH >= umbral_HH   ~ "YES",
                                       TRUE                   ~ "NO"),
    cod_dism_subpob        = case_when(is.na(subpob_perdida) ~ "Unknown",
                                       subpob_perdida         ~ "YES",
                                       TRUE                   ~ "NO"),
    subpob_desap_sino      = case_when(is.na(subpob_perdida) ~ NA_character_,
                                       subpob_perdida         ~ "SI",
                                       TRUE                   ~ "NO"),
    fuente_dism_habitat    = ifelse(cod_dism_habitat == "Unknown", NA, "Inferred"),
    fuente_dism_subpob     = ifelse(cod_dism_subpob  == "Unknown", NA, "Inferred")
  )

dir.create("resultados/eecorisk/fragmentacion_severa", recursive = TRUE, showWarnings = FALSE)
dir.create("resultados/eecorisk/habitat_disponible",   recursive = TRUE, showWarnings = FALSE)

write.csv(Tablafrag, "resultados/eecorisk/fragmentacion_severa/resultados_eecorisk.csv",
          row.names = FALSE)

# Detalle de parches por especie
Tallfg <- Filter(Negate(is.null), lapply(seq_along(sg), function(i) {
  if (length(sg[[i]]) == 0) return(NULL)
  cbind(Especie = ne[i], as.data.frame(sg[[i]][[1]]))
}))
sapply(Tallfg, ncol)

for(i in seq_along(Tallfg)){
  names(Tallfg[[i]]) <- c(
    "Especie",
    "Area km^2",
    "Dist_PMC m",
    "Isolated",
    "Small"
  )
}

write.csv(do.call(rbind, Tallfg),
          "resultados/eecorisk/fragmentacion_severa/detalle_parches.csv",
          row.names = FALSE)

# Textos SIS ----
# Municipios y departamentos por especie (generado en script 02)
reg_mpios <- read.csv("resultados/ConR/criterioB/registros_municipios_dptos.csv",
                      encoding = "UTF-8")

mpio_dpto_sp <- function(sp) {
  d <- reg_mpios %>% filter(tax == sp, !is.na(municipio))
  if (nrow(d) == 0) return("municipios no disponibles")
  pares <- d %>% distinct(municipio, departamento) %>%
    mutate(txt = paste0(municipio, " (departamento de ", departamento, ")"))
  paste(pares$txt, collapse = ", ")
}

# Textos descriptivos SIS desde eecorisk
num_palabras <- function(n){
  
  palabras <- c(
    "una","dos","tres","cuatro","cinco",
    "seis","siete","ocho","nueve","diez"
  )
  
  sapply(n, function(x){
    
    if(is.na(x) || x == 0)
      return("ninguna")
    
    if(x >= 1 && x <= 10)
      return(palabras[x])
    
    as.character(x)
    
  })
  
}

subpop_res <- read.csv("resultados/ConR/subpoblaciones/subpoblaciones.csv")

Tablafrag <- Tablafrag %>%
  left_join(subpop_res %>% rename(n_subpop = subpop), by = "tax") %>%
  mutate(
    mpios = sapply(tax, mpio_dpto_sp),

    desc_frag = case_when(
      cod_fragmentacion != "YES" ~ NA_character_,
      TRUE ~ paste0("El ", FS_score, "% de parches de hábitat donde se encuentra la especie ",
                    "son pequeños y aislados. Estos parches se encuentran en el/los municipio(s) de ",
                    mpios, ".")
    ),

    desc_dism_hab = case_when(
      cod_dism_habitat != "YES" ~ NA_character_,
      TRUE ~ paste0(
        tools::toTitleCase(num_palabras(n_subpop)),
        " subpoblacion", ifelse(n_subpop == 1, " de", "es de"),
        " la especie se encuentran en paisajes con destrucción y degradación de su hábitat. ",
        "Estas subpoblaciones se encuentran en el/los municipio(s) de ",
        mpios, "."
      )
    ),

    desc_dism_subpob = case_when(
      cod_dism_subpob != "YES" ~ NA_character_,
      TRUE ~ paste0("Es posible que alguna(s) subpoblación(es) de la especie en ",
                    mpios, " haya(n) desaparecido por la destrucción de su hábitat.")
    ),

    tendencia         = ifelse(cod_dism_habitat == "YES", "Decreasing",
                               ifelse(cod_dism_habitat == "NO", "Stable", "Unknown")),
    fuente_tendencia  = ifelse(tendencia == "Unknown", "Unknown", "Inferred"),
    no_amenazas       = ifelse(cod_dism_habitat == "NO" & cod_dism_subpob == "NO",
                               "TRUE", "FALSE"),
    amenazas_descon   = ifelse(cod_dism_habitat == "Unknown", "TRUE", "FALSE")
  )

# Base maestra ----
Tablafrag$no_amenazas <- Tablafrag$no_amenazas == "TRUE"
Tablafrag$amenazas_descon <- Tablafrag$amenazas_descon == "TRUE"

base_maestra <- read.csv("SIS_Connect/base_maestra.csv",
                         encoding = "UTF-8", check.names = FALSE)
names(base_maestra) <- make.unique(names(base_maestra))

base_maestra <- base_maestra %>%
  left_join(dplyr::select(Tablafrag, tax, FS_score, pct_HH,
                           cod_fragmentacion, cod_dism_habitat, cod_dism_subpob,
                           subpob_desap_sino, fuente_dism_habitat, fuente_dism_subpob,
                           desc_frag, desc_dism_hab, desc_dism_subpob,
                           tendencia, fuente_tendencia, no_amenazas, amenazas_descon),
            by = c("NOMBRE CIENTÍFICO sin autor" = "tax")) %>%
  mutate(
    `% PARCHES PEQUEÑOS Y AISLADOS`                           = FS_score,
    `% HUELLA HUMANA EN LA AOO`                               = pct_HH,
    `REPORTE SUBPOBLACIONES DESAPARECIDAS`                    = subpob_desap_sino,
    `CÓDIGO SIS FRAGMENTACIÓN`                                = cod_fragmentacion,
    `DESCRIPCIÓN DE FRAGMENTACIÓN SIS`                        = desc_frag,
    `CÓDIGO SIS DISMINUCIÓN CONTINUA HÁBITAT`                 = cod_dism_habitat,
    `DESCRIPCIÓN DE DISMINUCIÓN CONTINUA HÁBITAT SIS`         = desc_dism_hab,
    `CÓDIGO SIS FUENTE DE LA DISM. CONTINUA HÁBITAT`          = fuente_dism_habitat,
    `CÓDIGO SIS DISMINUCIÓN CONTINUA SUBPOBLACIONES`          = cod_dism_subpob,
    `DESCRIPCIÓN DE DISMINUCIÓN CONTINUA SUBPOBLACIONES SIS`  = desc_dism_subpob,
    `CÓDIGO SIS FUENTE DE LA DISM. CONTINUA SUBPOBLACIONES`   = fuente_dism_subpob,
    `CÓDIGO SIS TENDENCIA POBLACIONAL`                        = tendencia,
    `CÓDIGO SIS FUENTE DE LA TENDENCIA POBLACIONAL`           = fuente_tendencia,
    `REPORTE DE "NO AMENAZAS" SIS`                            = no_amenazas,
    `REPORTE DE "AMENAZAS DESCONOCIDAS" SIS`                  = amenazas_descon
  ) %>%
  dplyr::select(-FS_score, -pct_HH, -cod_fragmentacion, -cod_dism_habitat,
                -cod_dism_subpob, -subpob_desap_sino, -fuente_dism_habitat, -fuente_dism_subpob,
                -desc_frag, -desc_dism_hab, -desc_dism_subpob,
                -tendencia, -fuente_tendencia, -no_amenazas, -amenazas_descon)


write.csv(base_maestra, "SIS_Connect/base_maestra.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# Guardar objetos intermedios para mapas de verificación (script 07)
dir.create("resultados/eecorisk/habitat_disponible", recursive = TRUE, showWarnings = FALSE)
saveRDS(AOOok, "resultados/eecorisk/habitat_disponible/AOOok.rds")
saveRDS(ne,    "resultados/eecorisk/habitat_disponible/ne.rds")

message("Script 04 completado.")
