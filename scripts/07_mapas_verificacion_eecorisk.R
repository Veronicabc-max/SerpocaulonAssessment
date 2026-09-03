# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 07: Mapas de verificación AOH y huella humana por especie
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Correr después del script 04

# Librerías ----
library(raster)
library(sf)
library(ggplot2)
library(ggspatial)   # annotation_map_tile() para basemap
library(dplyr)
library(cowplot)

# Datos ----
# Cargar objetos generados por script 04
AOOok     <- readRDS("resultados/eecorisk/habitat_disponible/AOOok.rds")
ne        <- readRDS("resultados/eecorisk/habitat_disponible/ne.rds")
registros <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8") %>%
  filter(!is.na(ddlat), !is.na(ddlon))
Tablafrag <- read.csv("resultados/eecorisk/fragmentacion_severa/resultados_eecorisk.csv")

# SB10 se usa como grilla de referencia para reprojectar IHEH al mismo CRS que el AOH
SB10 <- raster("datos/capas/coberturas_tierra/BNB_300_2024.tif")
hh   <- raster("datos/capas/coberturas_tierra/iheh_col.tif")
hh   <- projectRaster(hh, SB10)   # mismo CRS/extent que los rasters AOH

umbral_HH <- 40

dir.create("resultados/mapas/verificacion_eecorisk", recursive = TRUE, showWarnings = FALSE)

# Función mapa ----
mapa_eecorisk <- function(sp, guardar = TRUE) {

  i <- which(ne == sp)

  if (is.null(AOOok[[i]])) {
    message(sp, ": sin AOH calculada, se omite.")
    return(invisible(NULL))
  }

  aoh <- AOOok[[i]]

  # Puntos de la especie
  pts <- registros %>%
    filter(tax == sp) %>%
    st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326)

  # Huella humana recortada al AOH
  aoh_mask        <- aoh
  aoh_mask[aoh_mask == 0] <- NA
  hh_crop         <- crop(hh, aoh_mask)
  hh_mask         <- mask(hh_crop, aoh_mask)
  pct_hh_val      <- round(mean(getValues(hh_mask), na.rm = TRUE), 0)

  # Parámetros eecorisk de la especie
  params <- Tablafrag %>% filter(tax == sp)
  fs_score      <- if (nrow(params) > 0) params$FS_score[1]      else NA
  cod_frag      <- if (nrow(params) > 0) params$cod_fragmentacion[1] else NA
  cod_hab       <- if (nrow(params) > 0) params$cod_dism_habitat[1]  else NA

  # Convertir rasters a data.frame para ggplot
  aoh_df <- as.data.frame(aoh, xy = TRUE) %>%
    rename(valor = 3) %>%
    mutate(panel = "AOH (área de hábitat disponible)",
           fill_aoh = ifelse(valor == 1, "Hábitat", NA))

  hh_full_df <- as.data.frame(crop(hh, extent(aoh)), xy = TRUE) %>%
    rename(hh = 3) %>%
    mutate(panel = paste0("Huella humana  [promedio en AOH: ", pct_hh_val,
                          "%  |  umbral: ", umbral_HH, "%]"))

  # Extensión con ligero margen
  ext  <- extent(aoh)
  xlim <- c(ext@xmin - 0.1, ext@xmax + 0.1)
  ylim <- c(ext@ymin - 0.1, ext@ymax + 0.1)

  # Panel izquierdo: AOH
  # annotation_map_tile: descarga teselas OpenStreetMap (requiere internet la primera vez;
  # las guarda en caché local). type = "cartolight" da fondo gris claro y legible.
  # geom_tile con alpha = 0.6 permite ver el basemap por debajo del AOH.
  p_aoh <- ggplot() +
    annotation_map_tile(type = "cartolight", zoom = NULL, quiet = TRUE) +
    geom_tile(data = aoh_df %>% filter(!is.na(fill_aoh)),
              aes(x = x, y = y), fill = "#2d8b57", alpha = 0.7) +
    geom_sf(data = pts, color = "black", fill = "#f5e642",
            shape = 21, size = 2.5, stroke = 0.8) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.3) +
    labs(title = "AOH (área de hábitat disponible)",
         subtitle = paste(nrow(pts), "registros |",
                          sum(getValues(aoh) == 1, na.rm = TRUE), "celdas"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 9, face = "bold"))

  # Panel derecho: huella humana
  p_hh <- ggplot() +
    annotation_map_tile(type = "cartolight", zoom = NULL, quiet = TRUE) +
    geom_tile(data = hh_full_df %>% filter(!is.na(hh)),
              aes(x = x, y = y, fill = hh), alpha = 0.7) +
    scale_fill_gradientn(
      colours = c("#1a9641", "#ffffbf", "#d7191c"),
      values  = scales::rescale(c(0, umbral_HH, 100)),
      limits  = c(0, 100),
      name    = "HH (%)",
      guide   = guide_colorbar(barheight = 8)
    ) +
    geom_sf(data = pts, color = "black", fill = "white",
            shape = 21, size = 2.5, stroke = 0.8) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    labs(title = paste0("Huella humana  [promedio AOH: ", pct_hh_val,
                        "%  |  umbral: ", umbral_HH, "%]"),
         subtitle = paste("Dism. hábitat:", cod_hab,
                          " | Frag. severa:", cod_frag,
                          " (FS score:", fs_score, "%)"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(size = 9, face = "bold"))

  p <- plot_grid(p_aoh, p_hh, ncol = 2) +
    theme(plot.title = element_text(face = "italic")) +
    labs(title = sp)

  p <- ggdraw(p) +
    draw_label(sp, x = 0.5, y = 0.98, vjust = 1,
               fontface = "italic", size = 13)

  if (guardar) {
    nombre <- gsub(" ", "_", sp)
    ggsave(
      filename = paste0("resultados/mapas/verificacion_eecorisk/", nombre, "_eecorisk.png"),
      plot = p, width = 14, height = 7, dpi = 150
    )
    message("Guardado: ", nombre, "_eecorisk.png")
  }

  return(p)
}

# Prueba ----
sp_prueba <- ne[!sapply(AOOok, is.null)][1]
mapa_eecorisk(sp_prueba, guardar = FALSE)

# Loop mapas ----
for (sp in ne) {
  tryCatch(
    mapa_eecorisk(sp, guardar = TRUE),
    error = function(e) message("Error en ", sp, ": ", e$message)
  )
}

message("Mapas de verificación guardados en resultados/mapas/verificacion_eecorisk/")
