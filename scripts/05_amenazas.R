# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 05: Análisis de amenazas y códigos SIS
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
#
# CAPAS DE AMENAZAS INCLUIDAS:
#   ✓ Petróleo y gas (ANH): datos/capas/amenazas/Tierras_Junio_170625.shp
#     Fuente: Agencia Nacional de Hidrocarburos (ANH), junio 2025
#     Código SIS: 3.1 - Oil & gas drilling
#
#   ✓ Minería de metales (ANM): datos/capas/amenazas/Titulo_vigente_030122.shp
#     Fuente: Agencia Nacional de Minería (ANM), enero 2022
#     NOTA: Actualizar con datos más recientes desde https://www.anm.gov.co
#     Código SIS: 3.2 - Mining & quarrying
#
#   ✓ Vías (OpenStreetMap): descargado automáticamente con osmdata
#     Código SIS: 4.1 - Roads & railroads
#
#   ⏳ Cobertura de la tierra (IDEAM 2022): datos/capas/coberturas_tierra/cobertura_2022/
#     Descargar desde el portal IDEAM: https://experience.arcgis.com/experience/568ddab184334f6b81a04d2fe9aac262
#     Buscar: "Cobertura de la Tierra 100K Periodo 2022 limite administrativo"
#     Códigos SIS derivados: 2.1.4 (cultivos), 2.3.4 (ganadería), 1.1 (urbano)
#     PENDIENTE: agregar cuando esté disponible la capa
#
#   ✗ Incendios (FIRMS/NASA): no se incluyó porque los datos de FIRMS solo están
#     disponibles para los últimos 7 días de forma directa. Para datos históricos
#     (últimos 10 años) se requiere registro en NASA Earthdata y descarga del archivo:
#     https://firms.modaps.eosdis.nasa.gov/country/
#     Alternativa: MapBiomas Fuego https://mapbiomas.org/fire (cicatrices anuales)
#     Esta amenaza es de baja relevancia para helechos epífitos de bosque húmedo.

library(sf)
library(dplyr)
library(osmdata)

sf_use_s2(FALSE)

registros  <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8") %>%
  filter(!is.na(ddlat), !is.na(ddlon))
puntos_sf  <- st_as_sf(registros, coords = c("ddlon", "ddlat"), crs = 4326)
reg_mpios  <- read.csv("resultados/ConR/criterioB/registros_municipios_dptos.csv",
                       encoding = "UTF-8")
subpop_res <- read.csv("resultados/ConR/subpoblaciones/subpoblaciones.csv")
base_eeco  <- read.csv("resultados/eecorisk/fragmentacion_severa/resultados_eecorisk.csv",
                       encoding = "UTF-8")

# Función: % subpoblaciones afectadas y código de alcance
alcance_sis <- function(pct_afectadas) {
  case_when(
    is.na(pct_afectadas)    ~ "Unknown",
    pct_afectadas > 90      ~ "Whole (>90%)",
    pct_afectadas >= 50     ~ "Majority (50-90%)",
    TRUE                    ~ "Minority (<50%)"
  )
}

# Función: intersección de puntos con una capa de amenaza
pct_subpob_amenaza <- function(puntos, amenaza_sf, subpop_res) {
  inter <- st_intersection(puntos, amenaza_sf) %>%
    st_drop_geometry() %>%
    distinct(tax, .keep_all = FALSE) %>%
    count(tax, name = "n_subpob_amenazada")
  subpop_res %>%
    left_join(inter, by = "tax") %>%
    mutate(
      n_subpob_amenazada = replace_na(n_subpob_amenazada, 0),
      pct_afectadas      = round(100 * n_subpob_amenazada / subpop, 1)
    ) %>%
    dplyr::select(tax, n_subpob_amenazada, pct_afectadas)
}

# Petróleo y gas (ANH)
petroleo <- st_read("datos/capas/amenazas/Tierras_Junio_170625.shp", quiet = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  filter(CLASIFICAC == "ASIGNADA")   # solo contratos activos

af_petroleo <- pct_subpob_amenaza(puntos_sf, petroleo, subpop_res) %>%
  rename(pct_petroleo = pct_afectadas, n_petroleo = n_subpob_amenazada)

# Minería de metales (ANM)
# NOTA: capa de 2022, actualizar cuando esté disponible versión más reciente
mineria <- st_read("datos/capas/amenazas/Titulo_vigente_030122.shp", quiet = TRUE) %>%
  st_transform(4326) %>%
  st_make_valid() %>%
  filter(ESTADO == "Activo")

af_mineria <- pct_subpob_amenaza(puntos_sf, mineria, subpop_res) %>%
  rename(pct_mineria = pct_afectadas, n_mineria = n_subpob_amenazada)

# Vías (OpenStreetMap - descarga automática)
ruta_vias <- "datos/capas/amenazas/vias_colombia.gpkg"
if (!file.exists(ruta_vias)) {
  message("Descargando vías de OpenStreetMap...")
  bbox_col <- c(-79.0, -4.2, -66.8, 12.6)
  vias_osm <- opq(bbox = bbox_col) %>%
    add_osm_feature(key = "highway",
                    value = c("motorway","trunk","primary","secondary")) %>%
    osmdata_sf()
  vias <- vias_osm$osm_lines %>%
    dplyr::select(osm_id, highway) %>%
    st_transform(4326)
  st_write(vias, ruta_vias, quiet = TRUE)
} else {
  vias <- st_read(ruta_vias, quiet = TRUE)
}

# Buffer de 1 km alrededor de vías para capturar impacto
vias_buf <- st_buffer(vias, dist = 0.009)   # ~1 km en grados decimales

af_vias <- pct_subpob_amenaza(puntos_sf, vias_buf, subpop_res) %>%
  rename(pct_vias = pct_afectadas, n_vias = n_subpob_amenazada)

# Tabla de amenazas por especie
amenazas_sp <- subpop_res %>%
  left_join(af_petroleo, by = "tax") %>%
  left_join(af_mineria,  by = "tax") %>%
  left_join(af_vias,     by = "tax") %>%
  left_join(base_eeco %>% dplyr::select(tax, cod_dism_habitat), by = "tax")

# Construir códigos SIS
amenazas_sp <- amenazas_sp %>%
  mutate(
    tiene_petroleo = pct_petroleo > 0 & !is.na(pct_petroleo),
    tiene_mineria  = pct_mineria  > 0 & !is.na(pct_mineria),
    tiene_vias     = pct_vias     > 0 & !is.na(pct_vias),
    tiene_hab      = cod_dism_habitat == "YES"
  ) %>%
  rowwise() %>%
  mutate(
    codigos_lista = list(c(
      if (tiene_petroleo) "3.1",
      if (tiene_mineria)  "3.2",
      if (tiene_vias)     "4.1"
    )),
    alcances_lista = list(c(
      if (tiene_petroleo) alcance_sis(pct_petroleo),
      if (tiene_mineria)  alcance_sis(pct_mineria),
      if (tiene_vias)     alcance_sis(pct_vias)
    )),

    cod_amenazas     = if (length(codigos_lista) == 0) "N/A" else
                         paste(codigos_lista,  collapse = ", "),
    cod_tiempo       = if (length(codigos_lista) == 0) "N/A" else
                         paste(rep("Ongoing", length(codigos_lista)), collapse = ", "),
    cod_alcance      = if (length(alcances_lista) == 0) "N/A" else
                         paste(alcances_lista, collapse = ", "),
    cod_severidad    = if (length(codigos_lista) == 0) "N/A" else
                         paste(rep("Very rapid declines", length(codigos_lista)), collapse = ", "),
    cod_presiones    = if (length(codigos_lista) == 0) "N/A" else
                         paste(rep("1.1 | 1.2 | 2.1 | 2.2", length(codigos_lista)), collapse = " // ")
  ) %>%
  ungroup()

# Texto descripción de amenazas
mpio_dpto_sp <- function(sp) {
  d <- reg_mpios %>% filter(tax == sp, !is.na(municipio))
  if (nrow(d) == 0) return("municipios no disponibles")
  pares <- d %>% distinct(municipio, departamento) %>%
    mutate(txt = paste0(municipio, " (departamento de ", departamento, ")"))
  paste(pares$txt, collapse = ", ")
}

num_palabras <- function(n) {
  if (is.na(n) || n == 0) return("ninguna")
  p <- c("una","dos","tres","cuatro","cinco","seis","siete","ocho","nueve","diez")
  if (n >= 1 && n <= 10) p[n] else as.character(n)
}

amenazas_sp <- amenazas_sp %>%
  mutate(mpios = sapply(tax, mpio_dpto_sp)) %>%
  mutate(desc_amenazas = case_when(
    !tiene_hab & !tiene_petroleo & !tiene_mineria & !tiene_vias ~
      "Las subpoblaciones conocidas de la especie se encuentran en hábitats poco perturbados por actividades humanas.",
    TRUE ~ paste0(
      tools::toTitleCase(num_palabras(subpop)), " subpoblacion",
      ifelse(subpop == 1, " conocida", "es conocidas"),
      " de la especie se encuentran en sitios con destrucción y degradación de su hábitat. ",
      "Las alteraciones del hábitat se deben principalmente a",
      ifelse(tiene_petroleo, " la minería de hidrocarburos (petróleo y/o gas),", ""),
      ifelse(tiene_mineria,  " la minería de metales,", ""),
      ifelse(tiene_vias,     " la construcción de infraestructura (vías),", ""),
      " en los municipios de ", mpios, "."
    )
  ))

write.csv(amenazas_sp, "resultados/amenazas/tabla_amenazas.csv", row.names = FALSE)
dir.create("resultados/amenazas", recursive = TRUE, showWarnings = FALSE)
write.csv(amenazas_sp, "resultados/amenazas/tabla_amenazas.csv", row.names = FALSE)

# Actualizar base_maestra.csv
base_maestra <- read.csv("SIS_Connect/base_maestra.csv",
                         encoding = "UTF-8", check.names = FALSE)
names(base_maestra) <- make.unique(names(base_maestra))

base_maestra <- base_maestra %>%
  left_join(dplyr::select(amenazas_sp, tax, desc_amenazas, cod_amenazas,
                           cod_tiempo, cod_alcance, cod_severidad, cod_presiones),
            by = c("NOMBRE CIENTÍFICO sin autor" = "tax")) %>%
  mutate(
    `DESCRIPCIÓN AMENAZAS`                     = coalesce(`DESCRIPCIÓN AMENAZAS`, desc_amenazas),
    `CÓDIGO SIS AMENAZAS`                      = coalesce(`CÓDIGO SIS AMENAZAS`, cod_amenazas),
    `CÓDIGO SIS TIEMPO DE LAS AMENAZAS`        = coalesce(`CÓDIGO SIS TIEMPO DE LAS AMENAZAS`, cod_tiempo),
    `CÓDIGO SIS ALCANCE DE LAS AMENAZAS`       = coalesce(`CÓDIGO SIS ALCANCE DE LAS AMENAZAS`, cod_alcance),
    `CÓDIGO SIS SEVERIDAD DE LAS AMENAZAS`     = coalesce(`CÓDIGO SIS SEVERIDAD DE LAS AMENAZAS`, cod_severidad),
    `CÓDIGO SIS PRESIONES RESULTADO DE AMENAZAS` = coalesce(
      `CÓDIGO SIS PRESIONES RESULTADO DE AMENAZAS`, cod_presiones)
  ) %>%
  dplyr::select(-desc_amenazas, -cod_amenazas, -cod_tiempo,
                -cod_alcance, -cod_severidad, -cod_presiones)

write.csv(base_maestra, "SIS_Connect/base_maestra.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

message("Script 05 completado.")
