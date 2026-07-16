# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 02: Cálculo de parámetros ConR y evaluación Criterio B
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Referencia: GEPC - Grupo de Especialistas en Plantas de Colombia

library(sf)
library(sp)
library(terra)
library(raster)
library(lwgeom)
library(ConR)
library(writexl)
library(geodata)   # GADM - límites administrativos detallados
library(wdpar)
library(dplyr)

# Crear carpetas de output si no existen
dirs_out <- c("resultados/ConR/EOO",
              "resultados/ConR/formas_EOO",
              "resultados/ConR/AOO",
              "resultados/ConR/subpoblaciones",
              "resultados/ConR/criterioB")
invisible(lapply(dirs_out, dir.create, recursive = TRUE, showWarnings = FALSE))

# Cargar registros limpios
# ConR requiere columnas en este orden: ddlat, ddlon, tax
registros <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8")
MyData <- registros[, c("ddlat", "ddlon", "tax")]

# Capa de Colombia (GADM nivel 0 = límite nacional)
# Más precisa que rnaturalearth; se descarga una sola vez
ruta_colombia <- "datos/capas/pais/Colombia.gpkg"

if (!file.exists(ruta_colombia)) {
  colombia <- gadm("Colombia", level = 0, path = tempdir(), version = "latest") %>%
    st_as_sf()
  st_write(colombia, ruta_colombia)
} else {
  colombia <- st_read(ruta_colombia)
}

# Áreas protegidas Colombia (fuente: Protected Planet / WDPA)
# wdpar descarga el archivo oficial; puede tardar unos minutos la primera vez
# Se recortan al extent de Colombia para reducir el peso del archivo
ruta_ap <- "datos/capas/areas_protegidas/wdpa_colombia.gpkg"

if (!file.exists(ruta_ap)) {
  sf_use_s2(FALSE)
  raw_pa <- wdpa_fetch("Colombia", wait = TRUE)

  # wdpa_clean() falla en Colombia por un polígono con geometría degenerada.
  # Se replica manualmente el filtrado esencial y se usa st_buffer(0)
  # como alternativa al paso de reparación que produce el error.
  areas_prot <- raw_pa %>%
    filter(STATUS %in% c("Designated", "Inscribed", "Established", "Adopted")) %>%
    filter(!grepl("Biosphere Reserve", DESIG, ignore.case = TRUE)) %>%
    filter(REP_AREA > 0 | GIS_AREA > 0) %>%
    st_buffer(0) %>%
    st_transform(st_crs(colombia)) %>%
    st_crop(st_bbox(colombia))

  sf_use_s2(TRUE)
  st_write(areas_prot, ruta_ap)
} else {
  areas_prot <- st_read(ruta_ap)
}

# Alinear CRS (la versión nueva de ConR usa sf directamente, no Spatial)
areas_prot <- st_transform(areas_prot, st_crs(colombia))

# Calcular EOO por especie
# Los shapefiles de EOO se guardan en el directorio de trabajo; luego se mueven
eoo <- EOO.computing(MyData,
                     country_map   = colombia,
                     export_shp    = TRUE,
                     write_shp     = TRUE,
                     show_progress = TRUE)

write.csv(eoo$results, "resultados/ConR/EOO/EOO_resultados.csv", row.names = FALSE)

# ConR guarda los shapefiles en shapesIUCN/ — moverlos a su carpeta
shp_eoo <- list.files("shapesIUCN", full.names = TRUE)
if (length(shp_eoo) > 0)
  file.rename(shp_eoo, file.path("resultados/ConR/formas_EOO", basename(shp_eoo)))

# Calcular AOO por especie (celda 2x2 km, estándar IUCN)
aoo <- AOO.computing(MyData,
                     cell_size_AOO    = 2,
                     nbe.rep.rast.AOO = 30,
                     show_progress    = TRUE)

write.csv(aoo, "resultados/ConR/AOO/AOO_resultados.csv", row.names = FALSE)

# Calcular número de localidades
loc <- locations.comp(MyData,
                      method              = "fixed_grid",
                      nbe_rep             = 30,
                      cell_size_locations = 10,
                      method_polygons     = "no_more_than_one",
                      show_progress       = TRUE)

write.csv(loc$locations, "resultados/ConR/criterioB/localidades.csv", row.names = FALSE)

# Calcular subpoblaciones (resolución 5 km)
subpop <- subpop.comp(MyData,
                      resol_sub_pop = 5,
                      show_progress = TRUE)

write.csv(subpop, "resultados/ConR/subpoblaciones/subpoblaciones.csv", row.names = FALSE)

# Evaluación Criterio B completa
# criterion_B calcula todo internamente para garantizar consistencia entre parámetros.
# Los cálculos individuales (EOO, AOO, loc, subpop) ya están guardados arriba como referencia.
# DrawMap = FALSE porque no es compatible con subpops; los mapas se generan en script 04.
# Combinar resultados individuales en una tabla resumen para SIS Connect
# La categoría IUCN se asigna manualmente por los evaluadores; criterion_B no es necesario.
resumen_conr <- eoo$results %>%
  rename(EOO_km2 = eoo) %>%
  full_join(aoo %>% rename(AOO_km2 = aoo), by = "tax") %>%
  full_join(loc$locations %>% rename(n_localidades = locations), by = "tax") %>%
  full_join(subpop %>% rename(n_subpoblaciones = subpop), by = "tax")

write_xlsx(resumen_conr, "resultados/ConR/criterioB/resumen_parametros_ConR.xlsx")

# Intersección de registros con áreas protegidas
# Para reportar en qué áreas se encuentra cada especie (útil para SIS Connect)
sf_use_s2(FALSE)
areas_prot <- st_make_valid(areas_prot)

puntos_sf <- registros %>%
  filter(!is.na(ddlat), !is.na(ddlon)) %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326)

interseccion <- st_intersection(puntos_sf, areas_prot) %>%
  st_drop_geometry()

write_xlsx(interseccion, "resultados/ConR/criterioB/registros_en_areas_protegidas.xlsx")

# Porcentaje de registros en áreas protegidas por especie
pct_ap <- interseccion %>%
  count(tax, name = "n_en_ap") %>%
  left_join(count(registros, tax, name = "n_total"), by = "tax") %>%
  mutate(pct_en_ap = round(n_en_ap / n_total * 100, 1))

# Listado de áreas protegidas por especie (nombres únicos separados por coma)
listado_ap <- interseccion %>%
  group_by(tax) %>%
  summarise(areas_protegidas = paste(sort(unique(NAME)), collapse = "; "), .groups = "drop")

# Intersección de registros con municipios y departamentos (GADM nivel 2)
# Se guarda para reusar en scripts 04, 05 y 06
ruta_mpios <- "datos/capas/pais/Colombia_mpios.gpkg"
if (!file.exists(ruta_mpios)) {
  mpios <- gadm("Colombia", level = 2, path = tempdir(), version = "latest") %>% st_as_sf()
  st_write(mpios, ruta_mpios)
} else {
  mpios <- st_read(ruta_mpios, quiet = TRUE)
}

reg_mpios <- puntos_sf %>%
  st_join(mpios %>% dplyr::select(municipio = NAME_2, departamento = NAME_1)) %>%
  st_drop_geometry() %>%
  dplyr::select(tax, id, municipio, departamento)

write.csv(reg_mpios, "resultados/ConR/criterioB/registros_municipios_dptos.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# Texto tamaño poblacional por especie (plantilla SIS)
num_palabras <- function(n) {
  if (is.na(n) || n == 0) return("ninguna")
  palabras <- c("una","dos","tres","cuatro","cinco","seis","siete","ocho","nueve","diez")
  if (n >= 1 && n <= 10) palabras[n] else as.character(n)
}

desc_tamano <- resumen_conr %>%
  mutate(desc_tamano_pob = paste0(
    "La especie tiene registradas hasta el momento ",
    sapply(n_subpoblaciones, num_palabras),
    " subpoblacion", ifelse(n_subpoblaciones == 1, ".", "es."),
    " No se conoce nada sobre la abundancia o tendencia poblacional de la especie."
  )) %>%
  dplyr::select(tax, desc_tamano_pob)

# Actualizar base_maestra.csv con resultados de ConR
base_maestra <- read.csv("SIS_Connect/base_maestra.csv",
                         encoding = "UTF-8", check.names = FALSE)
# La BASE MAESTRA tiene columnas duplicadas; make.unique las hace únicas para el join
names(base_maestra) <- make.unique(names(base_maestra))

# Los nombres de columna exactos del output de criterion_B pueden variar según versión;
# ajustar si es necesario revisando names(criterioB)
base_maestra <- base_maestra %>%
  left_join(dplyr::select(resumen_conr, tax, EOO_km2, AOO_km2,
                          n_loc = n_localidades, n_subpop = n_subpoblaciones),
            by = c("NOMBRE CIENTÍFICO sin autor" = "tax")) %>%
  left_join(dplyr::select(pct_ap, tax, pct_en_ap),
            by = c("NOMBRE CIENTÍFICO sin autor" = "tax")) %>%
  left_join(dplyr::select(listado_ap, tax, areas_protegidas),
            by = c("NOMBRE CIENTÍFICO sin autor" = "tax")) %>%
  left_join(desc_tamano, by = c("NOMBRE CIENTÍFICO sin autor" = "tax")) %>%
  mutate(
    `EOO (km2)`                            = coalesce(`EOO (km2)`, EOO_km2),
    `AOO (km2)`                            = coalesce(`AOO (km2)`, AOO_km2),
    `# LOCALIDADES "locations"`            = coalesce(`# LOCALIDADES "locations"`, n_loc),
    `# SUBPOBLACIONES`                     = coalesce(`# SUBPOBLACIONES`, n_subpop),
    `% OCURRENCIAS EN AREAS PROTEGIDAS`    = coalesce(`% OCURRENCIAS EN AREAS PROTEGIDAS`, pct_en_ap),
    `LISTADO DE AREAS PROTEGIDAS CON OCURRENCIAS` = coalesce(
      `LISTADO DE AREAS PROTEGIDAS CON OCURRENCIAS`, areas_protegidas),
    `DESCRIPCIÓN TAMAÑO POBLACIONAL Y DEMOGRAFÍA` = coalesce(
      `DESCRIPCIÓN TAMAÑO POBLACIONAL Y DEMOGRAFÍA`, desc_tamano_pob),
    `REPORTE DE PRESENCIA EN AREAS PROTEGIDAS`    = coalesce(
      `REPORTE DE PRESENCIA EN AREAS PROTEGIDAS`,
      ifelse(is.na(pct_en_ap) | pct_en_ap == 0, "NO", "YES"))
  ) %>%
  dplyr::select(-EOO_km2, -AOO_km2, -n_loc, -n_subpop, -pct_en_ap,
                -areas_protegidas, -desc_tamano_pob)

write.csv(base_maestra, "SIS_Connect/base_maestra.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
