# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 03: Fragmentación severa, huella humana y disminución de hábitat (eecorisk)
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Referencia: GEPC - Grupo de Especialistas en Plantas de Colombia
#
# RASTERS REQUERIDOS - descargar manualmente y ubicar en datos/capas/:
#
#   coberturas_tierra/   → Capas de Bosque/No Bosque del IDEAM
#     Descargar desde: https://experience.arcgis.com/experience/568ddab184334f6b81a04d2fe9aac262
#     Buscar: "Mapa de Bosque No Bosque" para los años 1990, 2000 y 2024
#     Guardar como: BNB_1990.tif, BNB_2000.tif, BNB_2024.tif
#     Luego remuestrear a 300m (ver sección "Preparar rasters BNB" abajo)
#
#   coberturas_tierra/iheh_col.tif   → Índice de Huella Espacial Humana (IHEH) Colombia
#     Cargar el archivo que tengas disponible con este nombre
#
# PARÁMETROS PARA Serpocaulon (helechos epífitos):
#   umbral = 150 m  (parche pequeño para epífitas/hierbas)
#   disper = 50 km  (dispersión por viento - esporas)
#   umbral_HH = 40  (% huella humana a partir del cual se declara disminución continua de hábitat)

library(raster)
library(sf)
library(geosphere)
library(sp)
library(elevatr)   # descarga DEM automáticamente
library(dplyr)
library(writexl)

# Cargar registros y agregar parámetros de especie
registros <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8") %>%
  filter(!is.na(ddlat), !is.na(ddlon)) %>%
  mutate(umbral = 150,
         disper = 50)

# DEM de elevación - resolución baja (z=6, ~1km) solo para filtro espacial en AHO
# Los valores mín/máx por especie se calculan desde elev_msnm del CSV (más preciso)
ruta_elv <- "datos/capas/elevacion/dem_colombia.tif"

if (!file.exists(ruta_elv)) {
  pts_sf <- st_as_sf(registros, coords = c("ddlon", "ddlat"), crs = 4326)
  elv <- get_elev_raster(pts_sf, z = 6, src = "aws")
  writeRaster(elv, ruta_elv, overwrite = TRUE)
} else {
  elv <- raster(ruta_elv)
}

# Rasters BNB del IDEAM
# Descargar manualmente desde el portal IDEAM (ver instrucciones arriba)
# y guardar en datos/capas/coberturas_tierra/ con los nombres indicados.
#
# Preparar rasters BNB a resolución 300m (correr UNA VEZ si los descargaste en resolución original):
# Preparar BNB_300 (correr UNA VEZ después de descargar los originales)
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

# Stack histórico BNB para análisis de disminución continua (1990, 2000, 2024)
BNBstk <- stack(
  raster("datos/capas/coberturas_tierra/BNB_300_1990.tif"),
  raster("datos/capas/coberturas_tierra/BNB_300_2000.tif"),
  SB10
)

# Huella humana (IHEH Colombia) - cargar el archivo disponible
hh <- raster("datos/capas/coberturas_tierra/iheh_col.tif")
hh <- projectRaster(hh, SB10)   # reproyectar a la misma resolución que SB10

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

################################################################################
# Función AHO original
################################################################################

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


# Funcion AHO optimizada
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

# Función sfrag: % parches pequeños y aislados
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
    return(list(cbind("Area km^2" = Area, "Isolated" = "Solo un parche",
                      "Small" = Area < unique(puntos[, 4])), "NULL", "NULL"))
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

# Calcular AHO por especie
pb <- txtProgressBar(min = 0,
                     max = length(ne),
                     style = 3)

for(i in seq_along(ne)){
  
  setTxtProgressBar(pb, i)
  
  message(ne[i])
  
  if(nrow(csp[[i]]) == 0) next
  
  cex <- extract(
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

# Fragmentación severa
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

# % Huella humana en el AHO por especie
pct_hh <- sapply(seq_along(ne), function(i) {
  
  if (is.null(AOOok[[i]]))
    return(NA)
  
  aoo_mask <- AOOok[[i]]
  
  aoo_mask[aoo_mask == 0] <- NA
  
  hh_crop <- crop(hh, aoo_mask)
  hh_mask <- mask(hh_crop, aoo_mask)
  
  round(mean(getValues(hh_mask), na.rm = TRUE), 0)
  
})

# Precalcular clumps históricos (solo una vez)
clumps_hist_list <- lapply(1:(nlayers(BNBstk) - 1), function(j) {
  
  clp <- clump(BNBstk[[j]], directions = 4)
  
  list(
    clump = clp,
    vals = getValues(clp)
  )
  
})

vals_current <- getValues(BNBstk[[nlayers(BNBstk)]])

# Función discon: detecta si alguna subpoblación (parche con ocurrencias) desapareció
# comparando capas históricas BNB con la actual (última capa del stack)
discon <- function(clumps_hist_list, vals_current, puntos, xy = c(3,2)) {
  
  coords <- as.matrix(puntos[, xy])
  
  for (j in seq_along(clumps_hist_list)) {
    
    clump_raster <- clumps_hist_list[[j]]$clump
    vals_clumps  <- clumps_hist_list[[j]]$vals
    
    vals_hist <- extract(clump_raster, coords)
    
    parches_hist <- unique(na.omit(vals_hist))
    
    if (length(parches_hist) == 0)
      next
    
    for (p in parches_hist) {
      
      idx <- vals_clumps == p
      
      if (all(vals_current[idx] == 0 | is.na(vals_current[idx])))
        return(TRUE)
      
    }
  }
  
  FALSE
}

# Detectar subpoblaciones desaparecidas por especie
subpob_perdida <- sapply(seq_along(ne), function(i) {
  
  if (nrow(csp[[i]]) == 0)
    return(NA)
  
  discon(
    clumps_hist_list = clumps_hist_list,
    vals_current = vals_current,
    puntos = csp[[i]],
    xy = c(3, 2)
  )
})

# Tabla de resultados
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
write.csv(do.call(rbind, Tallfg),
          "resultados/eecorisk/fragmentacion_severa/detalle_parches.csv",
          row.names = FALSE)

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
num_palabras <- function(n) {
  if (is.na(n) || n == 0) return("ninguna")
  p <- c("una","dos","tres","cuatro","cinco","seis","siete","ocho","nueve","diez")
  if (n >= 1 && n <= 10) p[n] else as.character(n)
}

subpop_res <- read.csv("resultados/ConR/subpoblaciones/subpoblaciones.csv")

Tablafrag <- Tablafrag %>%
  left_join(subpop_res %>% rename(n_subpop = subpop), by = "tax") %>%
  mutate(
    mpios = sapply(tax, mpio_dpto_sp),

    desc_frag = case_when(
      cod_fragmentacion != "YES" ~ NA_character_,
      TRUE ~ paste0("El ", FS_score, "% de parches de hábitat donde se encuentra la especie ",
                    "son pequeños y aislados. Estos parches se encuentran en ",
                    mpios, ".")
    ),

    desc_dism_hab = case_when(
      cod_dism_habitat != "YES" ~ NA_character_,
      TRUE ~ paste0(tools::toTitleCase(num_palabras(n_subpop)),
                    " subpoblacion", ifelse(n_subpop == 1, " de", "es de"),
                    " la especie se encuentran en paisajes con destrucción y degradación ",
                    "de su hábitat. Estas subpoblaciones se encuentran en ", mpios, ".")
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

# Actualizar base_maestra.csv con todos los campos derivados de eecorisk
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
    `% PARCHES PEQUEÑOS Y AISLADOS`                               = coalesce(`% PARCHES PEQUEÑOS Y AISLADOS`, FS_score),
    `% HUELLA HUMANA EN LA AOO`                                   = coalesce(`% HUELLA HUMANA EN LA AOO`, pct_HH),
    `SUBPOBLACIONES DESAPARECIDAS`                                = coalesce(`SUBPOBLACIONES DESAPARECIDAS`, subpob_desap_sino),
    `CÓDIGO SIS FRAGMENTACIÓN`                                    = coalesce(`CÓDIGO SIS FRAGMENTACIÓN`, cod_fragmentacion),
    `DESCRIPCIÓN DE FRAGMENTACIÓN SIS`                            = coalesce(`DESCRIPCIÓN DE FRAGMENTACIÓN SIS`, desc_frag),
    `CÓDIGO SIS DISMINUCIÓN CONTINUA HÁBITAT`                    = coalesce(`CÓDIGO SIS DISMINUCIÓN CONTINUA HÁBITAT`, cod_dism_habitat),
    `DESCRIPCIÓN DE DISMINUCIÓN CONTINUA HÁBITAT SIS`            = coalesce(`DESCRIPCIÓN DE DISMINUCIÓN CONTINUA HÁBITAT SIS`, desc_dism_hab),
    `CÓDIGO SIS FUENTE DE LA DISM. CONTINUA HÁBITAT`             = coalesce(`CÓDIGO SIS FUENTE DE LA DISM. CONTINUA HÁBITAT`, fuente_dism_habitat),
    `CÓDIGO SIS DISMINUCIÓN CONTINUA SUBPOBLACIONES`             = coalesce(`CÓDIGO SIS DISMINUCIÓN CONTINUA SUBPOBLACIONES`, cod_dism_subpob),
    `DESCRIPCIÓN DE DISMINUCIÓN CONTINUA SUBPOBLACIONES SIS`     = coalesce(`DESCRIPCIÓN DE DISMINUCIÓN CONTINUA SUBPOBLACIONES SIS`, desc_dism_subpob),
    `CÓDIGO SIS FUENTE DE LA DISM. CONTINUA SUBPOBLACIONES`      = coalesce(`CÓDIGO SIS FUENTE DE LA DISM. CONTINUA SUBPOBLACIONES`, fuente_dism_subpob),
    `CÓDIGO SIS TENDENCIA POBLACIONAL`                           = coalesce(`CÓDIGO SIS TENDENCIA POBLACIONAL`, tendencia),
    `CÓDIGO SIS FUENTE DE LA TENDENCIA POBLACIONAL`              = coalesce(`CÓDIGO SIS FUENTE DE LA TENDENCIA POBLACIONAL`, fuente_tendencia),
    `REPORTE DE "NO AMENAZAS" SIS`                               = coalesce(`REPORTE DE "NO AMENAZAS" SIS`, no_amenazas),
    `REPORTE DE "AMENAZAS DESCONOCIDAS" SIS`                     = coalesce(`REPORTE DE "AMENAZAS DESCONOCIDAS" SIS`, amenazas_descon)
  ) %>%
  dplyr::select(-FS_score, -pct_HH, -cod_fragmentacion, -cod_dism_habitat,
                -cod_dism_subpob, -subpob_desap_sino, -fuente_dism_habitat, -fuente_dism_subpob,
                -desc_frag, -desc_dism_hab, -desc_dism_subpob,
                -tendencia, -fuente_tendencia, -no_amenazas, -amenazas_descon)

write.csv(base_maestra, "SIS_Connect/base_maestra.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

message("Script 04 completado.")
