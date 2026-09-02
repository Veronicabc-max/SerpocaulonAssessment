# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 02: Cálculo de parámetros ConR y evaluación Criterio B
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Referencia: GEPC - Grupo de Especialistas en Plantas de Colombia
#
# Este script calcula los tres parámetros del Criterio B de la IUCN:
#   EOO  - Extensión de la presencia (Extent of Occurrence)
#   AOO  - Área de occupación (Area of Occupancy)
#   Loc  - Número de localidades
#   Subpob - Número de subpoblaciones
# Todos los resultados se guardan en resultados/ConR/ y se transfieren
# automáticamente a SIS_Connect/base_maestra.csv

# Librerías ----
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

# Carpetas de output ----
# Crear carpetas de output si no existen
dirs_out <- c("resultados/ConR/EOO",
              "resultados/ConR/formas_EOO",
              "resultados/ConR/AOO",
              "resultados/ConR/subpoblaciones",
              "resultados/ConR/criterioB")
invisible(lapply(dirs_out, dir.create, recursive = TRUE, showWarnings = FALSE))

# Función subpop personalizada ----
# Función personalizada para calcular subpoblaciones asignando un ID a cada una.
# ConR::subpop.comp() solo retorna el conteo; esta versión retorna además
# a qué subpoblación pertenece cada registro individual, lo que permite
# intersectar con amenazas en los scripts 04 y 05.
# El radio (resol_sub_pop) define el tamaño del buffer alrededor de cada punto;
# dos puntos en la misma subpoblación son aquellos cuyos buffers se solapan.
subpop.comp.records <- function(
    XY,
    resol_sub_pop = NULL,
    proj_type = "cea"
){

  if(is.null(resol_sub_pop))
    stop("Debe indicar resol_sub_pop")

  proj_type <- ConR:::proj_crs(proj_type)

  if(is.data.frame(resol_sub_pop)){

    XY <- merge(
      XY,
      resol_sub_pop,
      by = "tax",
      all.x = TRUE,
      sort = FALSE
    )

    XY <- XY[,c("ddlat","ddlon","tax","radius")]

  }else{

    XY$radius <- resol_sub_pop

  }

  lista <- ConR:::coord.check(
    XY = XY,
    listing = TRUE,
    proj_type = proj_type
  )[[1]]

  salida <- vector("list", length(lista))

  for(i in seq_along(lista)){

    datos <- lista[[i]]

    puntos <- sf::st_as_sf(
      datos,
      coords = c("ddlon","ddlat"),
      crs = 4326
    )

    puntos <- sf::st_transform(
      puntos,
      proj_type
    )

    buffers <- sf::st_buffer(
      puntos,
      dist = unique(datos$radius) * 1000
    )

    buffers <- sf::st_union(buffers)

    buffers <- sf::st_cast(
      buffers,
      "POLYGON"
    )

    subpop <- sf::st_as_sf(
      data.frame(
        geometry = buffers
      )
    )

    subpop$subpop_id <- seq_len(nrow(subpop))

    relacion <- sf::st_intersects(
      puntos,
      subpop
    )

    puntos$subpop_id <-
      sapply(relacion, `[`, 1)

    puntos <- sf::st_drop_geometry(puntos)

    salida[[i]] <- puntos

  }

  registros <- dplyr::bind_rows(salida)

  resumen <-

    registros %>%

    dplyr::group_by(tax) %>%

    dplyr::summarise(
      subpop = dplyr::n_distinct(subpop_id),
      .groups = "drop"
    )

  list(
    number_subpop = resumen,
    registros = registros
  )

}

# Datos y capas base ----
# Cargar registros limpios
# ConR requiere columnas en este orden: ddlat, ddlon, tax
registros <- read.csv("datos/registros/registros_Amplia_Distrib.csv", encoding = "UTF-8")
MyData <- registros[, c("ddlat", "ddlon", "tax")]

# Límite nacional de Colombia (GADM nivel 0)
# Se usa GADM en lugar de rnaturalearth porque tiene mayor precisión costera,
# lo que mejora el recorte del EOO en especies costeras o de islas.
# Se descarga una sola vez y se guarda localmente.
ruta_colombia <- "datos/capas/pais/Colombia.gpkg"

if (!file.exists(ruta_colombia)) {
  colombia <- gadm("Colombia", level = 0, path = tempdir(), version = "latest") %>%
    st_as_sf()
  st_write(colombia, ruta_colombia)
} else {
  colombia <- st_read(ruta_colombia)
}

# Áreas protegidas de Colombia (fuente: Protected Planet / WDPA)
# wdpar descarga el shapefile oficial de UNEP-WCMC con todas las APs de Colombia.
# Se filtran solo las APs con estatus legal activo (Designated, Inscribed, etc.)
# y se excluyen las Reservas de Biosfera (no son APs de protección estricta).
# wdpa_clean() falla en Colombia por un polígono con geometría degenerada;
# se replica su lógica manualmente usando st_buffer(0) para reparar geometrías.
ruta_ap <- "datos/capas/areas_protegidas/wdpa_colombia.gpkg"

if (!file.exists(ruta_ap)) {
  sf_use_s2(FALSE)
  raw_pa <- wdpa_fetch("Colombia", wait = TRUE)

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

areas_prot <- st_transform(areas_prot, st_crs(colombia))

# EOO ----
# EOO - Extensión de la Presencia (Extent of Occurrence)
# La EOO es el área del polígono convexo mínimo (convex hull) que encierra
# todos los registros de la especie. Representa el rango geográfico potencial
# total de la especie, sin considerar hábitat o discontinuidades.
# Umbrales IUCN Criterio B1:
#   CR  < 100 km²
#   EN  < 5,000 km²
#   VU  < 20,000 km²
# country_map recorta el EOO al territorio colombiano, evitando que el polígono
# se extienda sobre el mar o países vecinos cuando hay puntos cerca de fronteras.
# export_shp = TRUE guarda el polígono EOO como shapefile para usar en los mapas.
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

# AOO ----
# AOO - Área de Ocupación (Area of Occupancy)
# La AOO cuenta el número de celdas de una grilla fija que contienen al menos
# un registro de la especie, multiplicado por el área de cada celda.
# Tamaño de celda = 2x2 km, estándar obligatorio de la IUCN.
# Umbrales IUCN Criterio B2:
#   CR  < 10 km²
#   EN  < 500 km²
#   VU  < 2,000 km²
# nbe.rep.rast.AOO = 30: la grilla se desplaza aleatoriamente 30 veces y se
# reporta el mínimo, cumpliendo el requisito IUCN de usar la posición más
# conservadora (menor AOO posible). 30 repeticiones son suficientes para
# estabilizar el estimado sin costo computacional excesivo.
aoo <- AOO.computing(MyData,
                     cell_size_AOO    = 2,
                     nbe.rep.rast.AOO = 30,
                     show_progress    = TRUE)

write.csv(aoo, "resultados/ConR/AOO/AOO_resultados.csv", row.names = FALSE)

# Localidades ----
# Una "localidad" es un área geográficamente distinta donde la especie puede
# ser afectada por una misma amenaza en un momento dado. Es un concepto
# cualitativo de la IUCN, pero ConR lo aproxima con una grilla.
# method = "fixed_grid": usa una grilla de celdas fijas (más reproducible que
#   el método de polígonos de Voronoi).
# cell_size_locations = 10 km: tamaño de celda para agrupar puntos en la misma
#   localidad. Se usa 10 km porque es la distancia a la que una amenaza típica
#   (e.g., un frente de deforestación) afectaría a todos los individuos del área.
#   El equipo GEPC recomienda este valor para plantas vasculares colombianas.
# nbe_rep = 30: igual que AOO, se usa el mínimo de 30 posiciones de grilla.
# method_polygons = "no_more_than_one": cada polígono de amenaza puede contener
#   a lo sumo una localidad, evitando sobreestimación.
# Umbrales IUCN Criterio B (localidades):
#   CR  ≤ 1
#   EN  ≤ 5
#   VU  ≤ 10
loc <- locations.comp(MyData,
                      method              = "fixed_grid",
                      nbe_rep             = 30,
                      cell_size_locations = 10,
                      method_polygons     = "no_more_than_one",
                      show_progress       = TRUE)

write.csv(loc$locations, "resultados/ConR/criterioB/localidades.csv", row.names = FALSE)

# Subpoblaciones ----
# Una subpoblación es un grupo de individuos separados geográficamente de otros
# grupos, con poco o ningún intercambio demográfico o genético.
# ConR las aproxima creando buffers circulares alrededor de cada registro y
# uniendo los que se solapan; cada componente conectado = una subpoblación.
# resol_sub_pop = 5 km: radio del buffer por registro. Dos registros a menos de
#   10 km de distancia (2 x 5 km) quedan en la misma subpoblación.
#   Para Serpocaulon (helechos epífitos dispersados por esporas) se usa 5 km
#   porque la dispersión por viento de esporas es efectiva a escala local
#   (~1-10 km), pero el flujo génico a mayor distancia es incierto.
#   Valores menores (ej. 2 km) fragmentarían demasiado; mayores (ej. 20 km)
#   unirían poblaciones que probablemente son independientes.
# export_shp = TRUE: guarda los polígonos de subpoblaciones para visualización.
subpop <- subpop.comp(
  MyData,
  resol_sub_pop = 5,
  export_shp = TRUE,
  show_progress = TRUE)

# Asignar cada registro a su subpoblación (necesario para scripts 04 y 05)
puntos_sf <- registros %>%
  filter(!is.na(ddlat), !is.na(ddlon)) %>%
  st_as_sf(
    coords = c("ddlon", "ddlat"),
    crs = 4326)

puntos_sf <- st_transform(
  puntos_sf,
  st_crs(subpop$poly_subpop))

subpop$poly_subpop$subpop_id <-
  seq_len(nrow(subpop$poly_subpop))

subpop$poly_subpop$subpop_id <- ave(
  seq_len(nrow(subpop$poly_subpop)),
  subpop$poly_subpop$tax,
  FUN = seq_along)

lista_registros <- split(puntos_sf, puntos_sf$tax)
lista_poligonos <- split(subpop$poly_subpop, subpop$poly_subpop$tax)

registros_subpop <-

  purrr::map_dfr(

    names(lista_registros),

    function(sp){

      puntos <- lista_registros[[sp]]

      poligonos <- lista_poligonos[[sp]]

      if(is.null(poligonos))
        return(NULL)

      salida <-

        sf::st_join(
          puntos,
          poligonos[, "subpop_id"],
          left = FALSE,
          join = sf::st_intersects
        )

      sf::st_drop_geometry(salida)

    }

  )

# Verificación: el conteo de subpoblaciones por especie debe coincidir
# entre la función original de ConR y la asignación por intersección.
# Si hay diferencias, revisar registros en los bordes de los polígonos.
a <- registros_subpop %>%
       group_by(tax) %>%
       summarise(
           subpop_calculadas = n_distinct(subpop_id)) %>%
       left_join(
           subpop$number_subpop,
           by = "tax" ) %>%
       mutate(coincide = subpop_calculadas == subpop)
View(a)

write.csv(registros_subpop,
  "resultados/ConR/subpoblaciones/subpoblaciones_registros.csv",
  row.names = FALSE)

write.csv(subpop$number_subpop,
  "resultados/ConR/subpoblaciones/subpoblaciones.csv",
  row.names = FALSE)

# Tabla resumen Criterio B ----
# Tabla resumen con todos los parámetros del Criterio B
# La categoría IUCN final la asignan los evaluadores manualmente considerando
# el conjunto de subcriterios (B1, B2) y condiciones adicionales (a, b, c).
# Este script no asigna categoría automáticamente.
resumen_conr <- eoo$results %>%
  rename(EOO_km2 = eoo) %>%
  full_join(
    aoo %>% rename(AOO_km2 = aoo),
    by = "tax") %>%
  full_join(
    loc$locations %>% rename(n_localidades = locations),
    by = "tax") %>%
  full_join(
    subpop$number_subpop %>%
      rename(n_subpoblaciones = subpop),
    by = "tax")

write_xlsx(resumen_conr, "resultados/ConR/criterioB/resumen_parametros_ConR.xlsx")

# Áreas protegidas ----
# Intersección con áreas protegidas
# Se usa sf_use_s2(FALSE) porque las geometrías WDPA tienen ocasionalmente
# vértices que cruzan el antimeridiano o son inválidos bajo geometría esférica (S2).
# st_make_valid() repara geometrías degeneradas antes de la intersección.
sf_use_s2(FALSE)
areas_prot <- st_make_valid(areas_prot)

puntos_sf <- registros %>%
  filter(!is.na(ddlat), !is.na(ddlon)) %>%
  st_as_sf(coords = c("ddlon", "ddlat"), crs = 4326)

interseccion <- st_intersection(puntos_sf, areas_prot) %>%
  st_drop_geometry()

write_xlsx(interseccion, "resultados/ConR/criterioB/registros_en_areas_protegidas.xlsx")

pct_ap <- interseccion %>%
  count(tax, name = "n_en_ap") %>%
  left_join(count(registros, tax, name = "n_total"), by = "tax") %>%
  mutate(pct_en_ap = round(n_en_ap / n_total * 100, 1))

listado_ap <- interseccion %>%
  group_by(tax) %>%
  summarise(areas_protegidas = paste(sort(unique(NAME)), collapse = "; "), .groups = "drop")

# Municipios y departamentos ----
# Intersección con municipios y departamentos (GADM nivel 2)
# Se guarda en disco para reutilizar en scripts 04, 05 y 06 sin recalcular.
ruta_mpios <- "datos/capas/pais/Colombia_mpios.gpkg"
if (!file.exists(ruta_mpios)) {
  mpios <- gadm("Colombia", level = 2, path = tempdir(), version = "latest") %>% st_as_sf()
  st_write(mpios, ruta_mpios)
} else {
  mpios <- st_read(ruta_mpios, quiet = TRUE)
}

reg_mpios <- puntos_sf %>%
  st_join(mpios, left = TRUE) %>%
  st_drop_geometry() %>%
  transmute(
    tax,
    id,
    municipio = NAME_2,
    departamento = NAME_1)

write.csv(reg_mpios, "resultados/ConR/criterioB/registros_municipios_dptos.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

# Texto tamaño poblacional SIS ----
# Texto de tamaño poblacional para SIS Connect
# Plantilla estándar del GEPC: reporta el número de subpoblaciones conocidas
# y declara que no se conoce abundancia ni tendencia (típico para herbario-based assessments).
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

# Base maestra ----
# Actualizar base_maestra.csv con resultados de ConR
# Todos los campos calculados se sobreescriben en cada corrida para garantizar
# que la base_maestra siempre refleje los resultados más recientes.
base_maestra <- read.csv("SIS_Connect/base_maestra.csv",
                         encoding = "UTF-8", check.names = FALSE)
names(base_maestra) <- make.unique(names(base_maestra))

# Verificar que todas las especies de registros estén en base_maestra; agregar las faltantes
especies_reg    <- unique(registros$tax)
especies_bm     <- base_maestra[["NOMBRE CIENTÍFICO sin autor"]]
especies_faltan <- setdiff(especies_reg, especies_bm)

if (length(especies_faltan) > 0) {
  message("Especies en registros pero no en base_maestra — se agregan automáticamente:")
  message(paste(" -", especies_faltan, collapse = "\n"))
  filas_nuevas <- data.frame(
    REINO   = "Plantae", PHYLLUM = "Tracheophyta", CLASE = "Polypodiopsida",
    ORDEN   = "Polypodiales", FAMILIA = "Polypodiaceae",
    `GÉNERO` = "Serpocaulon",
    ESPECIE = sub("Serpocaulon ", "", especies_faltan),
    `NOMBRE CIENTÍFICO sin autor` = especies_faltan,
    check.names = FALSE
  )
  base_maestra <- bind_rows(base_maestra, filas_nuevas) %>%
    arrange(`NOMBRE CIENTÍFICO sin autor`)
  write.csv(base_maestra, "SIS_Connect/base_maestra.csv",
            row.names = FALSE, fileEncoding = "UTF-8")
}

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
    `EOO (km2)`                                   = EOO_km2,
    `AOO (km2)`                                   = AOO_km2,
    `# LOCALIDADES "locations"`                   = n_loc,
    `# SUBPOBLACIONES`                            = n_subpop,
    `% OCURRENCIAS EN AREAS PROTEGIDAS`           = pct_en_ap,
    `LISTADO DE AREAS PROTEGIDAS CON OCURRENCIAS` = areas_protegidas,
    `DESCRIPCIÓN TAMAÑO POBLACIONAL Y DEMOGRAFÍA` = desc_tamano_pob,
    `REPORTE DE PRESENCIA EN AREAS PROTEGIDAS`    = ifelse(
      is.na(pct_en_ap) | pct_en_ap == 0, "NO", "YES")
  ) %>%
  dplyr::select(-EOO_km2, -AOO_km2, -n_loc, -n_subpop, -pct_en_ap,
                -areas_protegidas, -desc_tamano_pob)

write.csv(base_maestra, "SIS_Connect/base_maestra.csv",
          row.names = FALSE, fileEncoding = "UTF-8")
