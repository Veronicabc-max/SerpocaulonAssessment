# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 03: Mapas de registros por especie (verificación y reporte)
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Correr DESPUÉS del script 02 (necesita resultados ConR) y ANTES del script 04

library(sf)
library(ggplot2)
library(dplyr)
library(geodata)
library(ggspatial)   # escala y norte
library(cowplot)     # mapa inset
library(readxl)
library(mapview)     # mapa interactivo para verificación rápida

dir.create("resultados/mapas/por_especie", recursive = TRUE, showWarnings = FALSE)

registros <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8") %>%
  filter(!is.na(ddlat), !is.na(ddlon))

puntos_sf <- st_as_sf(registros, coords = c("ddlon", "ddlat"), crs = 4326)

# Capas base
ruta_colombia <- "datos/capas/pais/Colombia.gpkg"
if (!file.exists(ruta_colombia)) {
  colombia <- gadm("Colombia", level = 0, path = tempdir(), version = "latest") %>% st_as_sf()
  st_write(colombia, ruta_colombia)
} else {
  colombia <- st_read(ruta_colombia, quiet = TRUE)
}

ruta_deptos <- "datos/capas/pais/Colombia_deptos.gpkg"
if (!file.exists(ruta_deptos)) {
  deptos <- gadm("Colombia", level = 1, path = tempdir(), version = "latest") %>% st_as_sf()
  st_write(deptos, ruta_deptos)
} else {
  deptos <- st_read(ruta_deptos, quiet = TRUE)
}

# Resultados ConR para la caja de texto
eoo_res   <- read.csv("resultados/ConR/EOO/EOO_resultados.csv")
aoo_res   <- read.csv("resultados/ConR/AOO/AOO_resultados.csv")
loc_res   <- read.csv("resultados/ConR/criterioB/localidades.csv")
subp_res  <- read.csv("resultados/ConR/subpoblaciones/subpoblaciones.csv")
pct_ap    <- read.csv("SIS_Connect/base_maestra.csv", check.names = FALSE) %>%
  dplyr::select(tax = `NOMBRE CIENTÍFICO sin autor`,
                pct_ap = `% OCURRENCIAS EN AREAS PROTEGIDAS`,
                n_ap   = `# LOCALIDADES "locations"`)

params_conr <- eoo_res %>%
  rename(EOO_km2 = eoo) %>%
  left_join(aoo_res  %>% rename(AOO_km2   = aoo),       by = "tax") %>%
  left_join(loc_res  %>% rename(n_loc     = locations),  by = "tax") %>%
  left_join(subp_res %>% rename(n_subpop  = subpop),     by = "tax") %>%
  left_join(pct_ap,                                       by = "tax")

# Función para construir un mapa de especie
mapa_especie <- function(sp, guardar = TRUE) {
  pts_sp <- puntos_sf %>% filter(tax == sp)
  if (nrow(pts_sp) == 0) return(invisible(NULL))

  # Extent con buffer
  bbox   <- st_bbox(pts_sp)
  buf    <- 0.8
  xlim   <- c(bbox["xmin"] - buf, bbox["xmax"] + buf)
  ylim   <- c(bbox["ymin"] - buf, bbox["ymax"] + buf)

  # Parámetros ConR
  par_sp <- params_conr %>% filter(tax == sp)
  eoo_v  <- if (nrow(par_sp) > 0 && !is.na(par_sp$EOO_km2))  round(par_sp$EOO_km2, 0)  else "N/D"
  aoo_v  <- if (nrow(par_sp) > 0 && !is.na(par_sp$AOO_km2))  round(par_sp$AOO_km2, 0)  else "N/D"
  loc_v  <- if (nrow(par_sp) > 0 && !is.na(par_sp$n_loc))    par_sp$n_loc              else "N/D"
  sub_v  <- if (nrow(par_sp) > 0 && !is.na(par_sp$n_subpop)) par_sp$n_subpop           else "N/D"
  pct_v  <- if (nrow(par_sp) > 0 && !is.na(par_sp$pct_ap))   par_sp$pct_ap             else "N/D"

  texto_params <- paste0(
    "EOO = ", eoo_v, " km²\n",
    "AOO (2 km) = ", aoo_v, " km²\n",
    "Registros únicos = ", nrow(pts_sp), "\n",
    "Subpoblaciones (5 km) = ", sub_v, "\n",
    "Localidades (10 km) = ", loc_v, "\n",
    "% registros en APs = ", pct_v
  )

  # Mapa principal
  p_main <- ggplot() +
    geom_sf(data = deptos,   fill = "grey96", color = "grey75", linewidth = 0.2) +
    geom_sf(data = colombia, fill = NA,       color = "grey30", linewidth = 0.5)

  # EOO si existe
  nombre_shp <- file.path("resultados/ConR/formas_EOO",
                           paste0(gsub(" ", "_", sp), ".shp"))
  if (file.exists(nombre_shp)) {
    eoo_shp <- st_read(nombre_shp, quiet = TRUE)
    p_main <- p_main +
      geom_sf(data = eoo_shp, fill = "grey60", alpha = 0.3,
              color = "grey40", linewidth = 0.5, linetype = "dashed")
  }

  p_main <- p_main +
    geom_sf(data = pts_sp, shape = 21, fill = "#d73027",
            color = "black", size = 2, stroke = 0.4, alpha = 0.85) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.3,
                     pad_x = unit(0.3, "cm"), pad_y = unit(0.3, "cm")) +
    annotate("label", x = xlim[1] + 0.1, y = ylim[2] - 0.1,
             label = texto_params, hjust = 0, vjust = 1, size = 2.8,
             label.size = 0.3, fill = "white", alpha = 0.85, family = "mono") +
    labs(title = sp, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title    = element_text(face = "italic", size = 12),
      panel.grid    = element_line(color = "grey90", linewidth = 0.2),
      panel.border  = element_rect(color = "grey50", fill = NA, linewidth = 0.5)
    )

  # Mapa inset de Colombia con ubicación de la especie
  p_inset <- ggplot() +
    geom_sf(data = colombia, fill = "grey90", color = "grey50", linewidth = 0.3) +
    geom_sf(data = pts_sp,   color = "#d73027", size = 0.8) +
    theme_void() +
    theme(panel.border = element_rect(color = "grey50", fill = NA, linewidth = 0.5))

  # Combinar
  p_final <- ggdraw(p_main) +
    draw_plot(p_inset, x = 0.72, y = 0.02, width = 0.26, height = 0.28)

  if (guardar) {
    nombre_out <- file.path("resultados/mapas/por_especie",
                            paste0(gsub(" ", "_", sp), ".png"))
    ggsave(nombre_out, p_final, width = 8, height = 7, dpi = 150)
  }

  return(p_final)
}

# Mapa interactivo de todos los registros para verificación de coordenadas
mapview(puntos_sf, zcol = "tax", layer.name = "Especie",
        map.types = "Esri.WorldShadedRelief")

# Mapa de prueba - verificar que se ve bien antes de correr el loop
sp_prueba <- sort(unique(registros$tax))[1]
print(mapa_especie(sp_prueba, guardar = FALSE))

# Loop completo - correr cuando el mapa de prueba se vea bien
especies <- sort(unique(registros$tax))
pb <- txtProgressBar(min = 0, max = length(especies), style = 3)
for (i in seq_along(especies)) {
  setTxtProgressBar(pb, i)
  mapa_especie(especies[i], guardar = TRUE)
}
close(pb)
message("Mapas guardados en resultados/mapas/por_especie/")
