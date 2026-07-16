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

# Tabla base para eecorisk: especie, Latitud, Longitud, umbral, disper, MinElv, MaxElv
corsp <- registros %>%
  dplyr::select(especie = tax, Latitud = ddlat, Longitud = ddlon, umbral, disper) %>%
  left_join(elev_min, by = c("especie" = "tax")) %>%
  left_join(elev_max, by = c("especie" = "tax"))

ne   <- unique(corsp$especie)
dsp  <- lapply(ne, function(sp) corsp[corsp$especie == sp, ])
csp  <- lapply(dsp, na.omit)
csp1 <- csp
for (i in seq_along(ne)) {
  if (nrow(csp[[i]]) < 1) next
  coordinates(csp1[[i]]) <- c("Longitud", "Latitud")
}

# Función AHO: área de hábitat disponible filtrada por elevación y puntos
AHO <- function(model, puntos, xy, bufferSize = 0.054, bufferPoints = TRUE, elv) {
  if (!inherits(model, "RasterLayer")) stop("model debe ser RasterLayer")
  elvg  <- resample(crop(elv, model), model)
  elvgv <- getValues(elvg)
  model[which(elvgv < as.numeric(unique(puntos[, 6])) - 300)] <- 0
  model[which(elvgv > as.numeric(unique(puntos[, 7])) + 300)] <- 0
  groups  <- clump(model, directions = 4)
  pts_sf  <- st_as_sf(puntos, coords = names(puntos)[xy], crs = st_crs(model))
  pbf     <- if (bufferPoints) st_buffer(pts_sf, dist = bufferSize) else pts_sf
  pex     <- extract(groups, as(pbf, "Spatial"))
  parches <- na.omit(unique(pex[[1]]))
  pexok   <- match(getValues(groups), parches)
  groups[which(pexok > 0)]   <- 1
  groups[which(is.na(pexok))] <- 0
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
vs1000 <- getValues(SB1000); cells <- seq_along(vs1000)
AOOok  <- vector("list", length(ne))

pb <- txtProgressBar(min = 0, max = length(ne), style = 3)
for (i in seq_along(ne)) {
  setTxtProgressBar(pb, i)
  message(" ", ne[i])
  if (nrow(csp[[i]]) == 0) next
  cex     <- extract(SB1000, csp1[[i]], cell = TRUE)
  cells1  <- cells[na.omit(-cex[, 1])]
  SB1000p <- SB1000; SB1000p[cex[, 1]] <- 1; SB1000p[cells1] <- 0
  SP      <- rasterToPolygons(SB1000p, dissolve = TRUE, fun = function(x) x > 0)
  SB10E   <- crop(SB10, SP)
  SB10E[SB10E == 2] <- 0; SB10E[SB10E > 1] <- 1; SB10E[SB10E == 3] <- 0
  AOOok[[i]] <- AHO(SB10E, csp[[i]], xy = c(3, 2),
                    bufferSize = 6000, bufferPoints = TRUE, elv = elv)
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
  if (is.null(AOOok[[i]])) return(NA)
  aoo_mask <- AOOok[[i]]
  aoo_mask[aoo_mask == 0] <- NA
  hh_crop <- crop(hh, aoo_mask)
  hh_mask <- mask(hh_crop, aoo_mask)
  round(mean(getValues(hh_mask), na.rm = TRUE), 0)
})

# Función discon: detecta si alguna subpoblación (parche con ocurrencias) desapareció
# comparando capas históricas BNB con la actual (última capa del stack)
discon <- function(BNBstk, puntos, xy = c(3, 2)) {
  n_layers <- nlayers(BNBstk)
  current  <- BNBstk[[n_layers]]
  coords   <- puntos[, xy]
  perdida  <- FALSE
  for (j in 1:(n_layers - 1)) {
    hist_layer   <- BNBstk[[j]]
    clumps_hist  <- clump(hist_layer, directions = 4)
    vals_hist    <- extract(clumps_hist, coords)
    parches_hist <- na.omit(unique(vals_hist))
    if (length(parches_hist) == 0) next
    for (p in parches_hist) {
      cells_patch  <- which(getValues(clumps_hist) == p)
      vals_current <- getValues(current)[cells_patch]
      if (all(vals_current == 0 | is.na(vals_current))) {
        perdida <- TRUE
        break
      }
    }
    if (perdida) break
  }
  return(perdida)
}

# Detectar subpoblaciones desaparecidas por especie
subpob_perdida <- sapply(seq_along(ne), function(i) {
  if (nrow(csp[[i]]) == 0) return(NA)
  tryCatch(discon(BNBstk, csp[[i]], xy = c(3, 2)), error = function(e) NA)
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

# Actualizar base_maestra.csv con todos los campos derivados de eecorisk
base_maestra <- read.csv("SIS_Connect/base_maestra.csv",
                         encoding = "UTF-8", check.names = FALSE)
names(base_maestra) <- make.unique(names(base_maestra))

base_maestra <- base_maestra %>%
  left_join(dplyr::select(Tablafrag, tax, FS_score, pct_HH,
                           cod_fragmentacion, cod_dism_habitat, cod_dism_subpob,
                           subpob_desap_sino, fuente_dism_habitat, fuente_dism_subpob),
            by = c("NOMBRE CIENTÍFICO sin autor" = "tax")) %>%
  mutate(
    `% PARCHES PEQUEÑOS Y AISLADOS`                          = coalesce(`% PARCHES PEQUEÑOS Y AISLADOS`, FS_score),
    `% HUELLA HUMANA EN LA AOO`                              = coalesce(`% HUELLA HUMANA EN LA AOO`, pct_HH),
    `SUBPOBLACIONES DESAPARECIDAS`                           = coalesce(`SUBPOBLACIONES DESAPARECIDAS`, subpob_desap_sino),
    `CÓDIGO SIS FRAGMENTACIÓN`                               = coalesce(`CÓDIGO SIS FRAGMENTACIÓN`, cod_fragmentacion),
    `CÓDIGO SIS DISMINUCIÓN CONTINUA HÁBITAT`                = coalesce(`CÓDIGO SIS DISMINUCIÓN CONTINUA HÁBITAT`, cod_dism_habitat),
    `CÓDIGO SIS DISMINUCIÓN CONTINUA SUBPOBLACIONES`         = coalesce(`CÓDIGO SIS DISMINUCIÓN CONTINUA SUBPOBLACIONES`, cod_dism_subpob),
    `CÓDIGO SIS FUENTE DE LA DISM. CONTINUA HÁBITAT`         = coalesce(`CÓDIGO SIS FUENTE DE LA DISM. CONTINUA HÁBITAT`, fuente_dism_habitat),
    `CÓDIGO SIS FUENTE DE LA DISM. CONTINUA SUBPOBLACIONES`  = coalesce(`CÓDIGO SIS FUENTE DE LA DISM. CONTINUA SUBPOBLACIONES`, fuente_dism_subpob)
  ) %>%
  dplyr::select(-FS_score, -pct_HH, -cod_fragmentacion, -cod_dism_habitat,
                -cod_dism_subpob, -subpob_desap_sino, -fuente_dism_habitat, -fuente_dism_subpob)

write.csv(base_maestra, "SIS_Connect/base_maestra.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

message("Script 03 completado.")
