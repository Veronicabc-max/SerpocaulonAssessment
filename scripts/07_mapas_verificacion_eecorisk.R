# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 07: Mapas de verificación AOH y huella humana por especie
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Correr después del script 04

library(raster)
library(sf)
library(ggplot2)
library(ggspatial)
library(dplyr)
library(cowplot)

# Cargar objetos generados por script 04
AOOok     <- readRDS("resultados/eecorisk/habitat_disponible/AOOok.rds")
ne        <- readRDS("resultados/eecorisk/habitat_disponible/ne.rds")
hh        <- raster("datos/capas/coberturas_tierra/iheh_col.tif")
registros <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8") %>%
  filter(!is.na(ddlat), !is.na(ddlon))
Tablafrag <- read.csv("resultados/eecorisk/fragmentacion_severa/resultados_eecorisk.csv")

umbral_HH <- 40

dir.create("resultados/mapas/verificacion_eecorisk", recursive = TRUE, showWarnings = FALSE)

# Función: mapa AOH + huella humana para una especie
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
  aoh_df <- as.data.frame(aoh,  xy = TRUE) %>% filter(!is.na(layer))
  hh_df  <- as.data.frame(hh_mask, xy = TRUE) %>%
    rename(hh = 3) %>% filter(!is.na(hh))

  # Extensión con ligero margen
  ext  <- extent(aoh)
  xlim <- c(ext@xmin - 0.1, ext@xmax + 0.1)
  ylim <- c(ext@ymin - 0.1, ext@ymax + 0.1)

  # Texto de parámetros
  label_txt <- paste0(
    "HH promedio: ", pct_hh_val, "%  (umbral: ", umbral_HH, "%)\n",
    "Dism. hábitat: ", cod_hab, "\n",
    "Frag. severa: ", cod_frag, "  (FS score: ", fs_score, "%)"
  )

  p <- ggplot() +
    # AOH (verde)
    geom_raster(data = aoh_df %>% filter(layer == 1),
                aes(x = x, y = y), fill = "#2d8b57", alpha = 0.6) +
    # Huella humana sobre AOH (gradiente rojo)
    geom_raster(data = hh_df, aes(x = x, y = y, fill = hh), alpha = 0.8) +
    scale_fill_gradient(low = "#ffffcc", high = "#d73027",
                        name = "Huella\nhumana (%)") +
    # Umbral de huella humana como línea de referencia en la leyenda (visual)
    geom_point(data = pts, aes(geometry = geometry), stat = "sf_coordinates",
               color = "black", fill = "#f5e642", shape = 21, size = 3, stroke = 1) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.3) +
    annotation_north_arrow(location = "tr",
                           style = north_arrow_fancy_orienteering(
                             text_size = 8)) +
    annotate("label", x = xlim[1] + 0.05, y = ylim[2] - 0.05,
             label = label_txt, hjust = 0, vjust = 1, size = 3,
             fill = "white", alpha = 0.85, label.size = 0.3) +
    labs(title = sp,
         subtitle = paste(nrow(pts), "registros |",
                          sum(getValues(aoh) == 1, na.rm = TRUE), "celdas AOH"),
         x = NULL, y = NULL) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "italic"),
          legend.position = "right")

  if (guardar) {
    nombre <- gsub(" ", "_", sp)
    ggsave(
      filename = paste0("resultados/mapas/verificacion_eecorisk/", nombre, "_eecorisk.png"),
      plot = p, width = 10, height = 7, dpi = 150
    )
    message("Guardado: ", nombre, "_eecorisk.png")
  }

  return(p)
}

# Mapa de prueba: primera especie con AOH disponible
sp_prueba <- ne[!sapply(AOOok, is.null)][1]
mapa_eecorisk(sp_prueba, guardar = FALSE)

# Loop: generar mapas para todas las especies
for (sp in ne) {
  tryCatch(
    mapa_eecorisk(sp, guardar = TRUE),
    error = function(e) message("Error en ", sp, ": ", e$message)
  )
}

message("Mapas de verificación guardados en resultados/mapas/verificacion_eecorisk/")
