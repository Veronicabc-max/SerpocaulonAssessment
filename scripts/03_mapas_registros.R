# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 03: Mapas de registros por especie (verificación y reporte)
# Autora: Maria Judith Carmona | 2026
# Correr ANTES del script 04 para verificar que los registros sean correctos

library(sf)
library(ggplot2)
library(dplyr)
library(geodata)
library(mapview)   # mapa interactivo para verificación

dir.create("resultados/mapas/por_especie", recursive = TRUE, showWarnings = FALSE)

registros <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8") %>%
  filter(!is.na(ddlat), !is.na(ddlon))

puntos_sf <- st_as_sf(registros, coords = c("ddlon", "ddlat"), crs = 4326)

# Capa de Colombia
ruta_colombia <- "datos/capas/pais/Colombia.gpkg"
if (!file.exists(ruta_colombia)) {
  colombia <- gadm("Colombia", level = 0, path = tempdir(), version = "latest") %>% st_as_sf()
  st_write(colombia, ruta_colombia)
} else {
  colombia <- st_read(ruta_colombia, quiet = TRUE)
}

# Departamentos para contexto
ruta_deptos <- "datos/capas/pais/Colombia_deptos.gpkg"
if (!file.exists(ruta_deptos)) {
  deptos <- gadm("Colombia", level = 1, path = tempdir(), version = "latest") %>% st_as_sf()
  st_write(deptos, ruta_deptos)
} else {
  deptos <- st_read(ruta_deptos, quiet = TRUE)
}

# Mapa interactivo de todos los registros para verificación
# Colorea por especie; hacer zoom para revisar puntos sospechosos
mapview(puntos_sf, zcol = "tax", layer.name = "Especie",
        map.types = "Esri.WorldShadedRelief")

# Mapas estáticos por especie con EOO (si ya existe del script 02)
especies <- sort(unique(registros$tax))

# Mapa de prueba con la primera especie - verificar que se ve bien antes de correr el loop
sp_prueba <- especies[1]
pts_prueba <- puntos_sf %>% filter(tax == sp_prueba)
p_prueba <- ggplot() +
  geom_sf(data = deptos, fill = "grey95", color = "grey70", linewidth = 0.2) +
  geom_sf(data = colombia, fill = NA, color = "grey30", linewidth = 0.4) +
  geom_sf(data = pts_prueba, color = "#d73027", size = 1.8, alpha = 0.8) +
  coord_sf(xlim = c(-80, -66), ylim = c(-5, 13)) +
  labs(title    = sp_prueba,
       subtitle = paste0("n = ", nrow(pts_prueba), " registros — MAPA DE PRUEBA"),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "italic"))
print(p_prueba)

# Si el mapa se ve bien, correr el loop completo
# Si hay algo que ajustar (escala, colores, etc.) hacerlo aquí antes de continuar
pb <- txtProgressBar(min = 0, max = length(especies), style = 3)

for (i in seq_along(especies)) {
  sp <- especies[i]
  setTxtProgressBar(pb, i)

  pts_sp <- puntos_sf %>% filter(tax == sp)
  n_pts  <- nrow(pts_sp)

  p <- ggplot() +
    geom_sf(data = deptos, fill = "grey95", color = "grey70", linewidth = 0.2) +
    geom_sf(data = colombia, fill = NA, color = "grey30", linewidth = 0.4)

  # Agregar EOO si existe
  nombre_shp <- file.path("resultados/ConR/formas_EOO",
                          paste0(gsub(" ", "_", sp), ".shp"))
  if (file.exists(nombre_shp)) {
    eoo_shp <- st_read(nombre_shp, quiet = TRUE)
    p <- p + geom_sf(data = eoo_shp, fill = "#2166ac", alpha = 0.15,
                     color = "#2166ac", linewidth = 0.5)
  }

  p <- p +
    geom_sf(data = pts_sp, color = "#d73027", size = 1.8, alpha = 0.8) +
    coord_sf(xlim = c(-80, -66), ylim = c(-5, 13)) +
    labs(title    = sp,
         subtitle = paste0("n = ", n_pts, " registros"),
         x = NULL, y = NULL,
         caption  = "Polígono azul = EOO | Puntos rojos = registros") +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "italic"),
          panel.grid    = element_line(color = "grey90", linewidth = 0.2))

  nombre_out <- file.path("resultados/mapas/por_especie",
                          paste0(gsub(" ", "_", sp), ".png"))
  ggsave(nombre_out, p, width = 7, height = 8, dpi = 150)
}

close(pb)
message("Mapas guardados en resultados/mapas/por_especie/")
