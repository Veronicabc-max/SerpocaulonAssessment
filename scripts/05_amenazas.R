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
#   ✓ Cobertura de la tierra (IDEAM 2022): datos/capas/coberturas_tierra/ECOSISTEMAS_18062025.gdb
#     Fuente: IDEAM – "Cobertura de la Tierra 100K Periodo 2022 limite administrativo"
#     Descargable en: https://experience.arcgis.com/experience/568ddab184334f6b81a04d2fe9aac262
#     Layer usada: e_cobertura_tierra_2022_admin
#     Campos: nivel_1 (categoría CLC principal), nivel_2 (subcategoría)
#       nivel_1 == "1" → territorios artificializados → SIS 1.1 (urbano)
#       nivel_2 %in% c("21","22","24") → cultivos y áreas agrícolas heterg. → SIS 2.1.4
#       nivel_2 == "23" → pastos → SIS 2.3.4
#     NOTA: archivo ~4 GB, no está en GitHub. Estudiante debe descargar y guardar en:
#       datos/capas/coberturas_tierra/ECOSISTEMAS_18062025.gdb
#
#   ✗ Incendios (FIRMS/NASA): no se incluyó porque los datos de FIRMS solo están
#     disponibles para los últimos 7 días de forma directa. Para datos históricos
#     (últimos 10 años) se requiere registro en NASA Earthdata y descarga del archivo:
#     https://firms.modaps.eosdis.nasa.gov/country/
#     Alternativa: MapBiomas Fuego https://mapbiomas.org/fire (cicatrices anuales)
#     Esta amenaza es de baja relevancia para helechos epífitos de bosque húmedo.

library(sf)
library(dplyr)
library(osmextract)
library(sf)

ruta_vias <- "datos/capas/amenazas/vias_colombia.gpkg"

if (!file.exists(ruta_vias)) {
  
  message("Descargando vías principales de Colombia...")
  
  vias <- oe_get(
    place = "Colombia",
    layer = "lines",
    force_download = TRUE,
    quiet = FALSE
  )
  
  vias <- vias[vias$highway %in%
                 c("motorway",
                   "trunk",
                   "primary",
                   "secondary"), ]
  
  st_write(vias, ruta_vias, delete_dsn = TRUE, quiet = TRUE)
  
} else {
  
  vias <- st_read(ruta_vias, quiet = TRUE)
  
}

sf_use_s2(FALSE)

registros  <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8") %>%
  filter(!is.na(ddlat), !is.na(ddlon))
puntos_sf  <- st_as_sf(registros, coords = c("ddlon", "ddlat"), crs = 4326)
ruta_mpios_csv <- "resultados/ConR/criterioB/registros_municipios_dptos.csv"
if (!file.exists(ruta_mpios_csv))
  stop("Archivo no encontrado: ", ruta_mpios_csv,
       "\nAsegúrate de haber corrido el script 02_ConR_criterioB.R primero.")
reg_mpios  <- read.csv(ruta_mpios_csv, encoding = "UTF-8")
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
      n_subpob_amenazada = dplyr::coalesce(n_subpob_amenazada, 0L),
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

# Buffer de 1 km alrededor de vías para capturar impacto
vias_buf <- st_buffer(vias, dist = 0.009)   # ~1 km en grados decimales

af_vias <- pct_subpob_amenaza(puntos_sf, vias_buf, subpop_res) %>%
  rename(pct_vias = pct_afectadas, n_vias = n_subpob_amenazada)

# Cobertura de la tierra (IDEAM 2022) - CLC Corine Land Cover
# Archivo ~4 GB, no está en GitHub; descargar desde portal IDEAM y guardar en:
# datos/capas/coberturas_tierra/ECOSISTEMAS_18062025.gdb
ruta_cob <- "datos/capas/coberturas_tierra/ECOSISTEMAS_18062025.gdb"

tiene_cobertura <- file.exists(ruta_cob)

if (tiene_cobertura) {
  cob <- st_read(ruta_cob, layer = "e_cobertura_tierra_2022_admin", quiet = TRUE) %>%
    st_transform(4326) %>%
    st_make_valid()

  cultivos <- cob %>% filter(nivel_2 %in% c("21", "22", "24"))
  ganaderia <- cob %>% filter(nivel_2 == "23")
  urbano    <- cob %>% filter(nivel_1 == "1")

  af_cultivos  <- pct_subpob_amenaza(puntos_sf, cultivos,  subpop_res) %>%
    rename(pct_cultivos  = pct_afectadas, n_cultivos  = n_subpob_amenazada)
  af_ganaderia <- pct_subpob_amenaza(puntos_sf, ganaderia, subpop_res) %>%
    rename(pct_ganaderia = pct_afectadas, n_ganaderia = n_subpob_amenazada)
  af_urbano    <- pct_subpob_amenaza(puntos_sf, urbano,    subpop_res) %>%
    rename(pct_urbano    = pct_afectadas, n_urbano    = n_subpob_amenazada)
} else {
  message("Capa de cobertura 2022 no encontrada en ", ruta_cob)
  message("Las amenazas de cobertura (2.1.4, 2.3.4, 1.1) no se calcularán.")
  especies <- unique(subpop_res$tax)
  af_cultivos  <- data.frame(tax = especies, pct_cultivos  = NA_real_, n_cultivos  = 0L)
  af_ganaderia <- data.frame(tax = especies, pct_ganaderia = NA_real_, n_ganaderia = 0L)
  af_urbano    <- data.frame(tax = especies, pct_urbano    = NA_real_, n_urbano    = 0L)
}

# Tabla de amenazas por especie
amenazas_sp <- subpop_res %>%
  left_join(af_petroleo,  by = "tax") %>%
  left_join(af_mineria,   by = "tax") %>%
  left_join(af_vias,      by = "tax") %>%
  left_join(af_cultivos,  by = "tax") %>%
  left_join(af_ganaderia, by = "tax") %>%
  left_join(af_urbano,    by = "tax") %>%
  left_join(base_eeco %>% dplyr::select(tax, cod_dism_habitat), by = "tax")

# Construir códigos SIS
amenazas_sp <- amenazas_sp %>%
  mutate(
    tiene_petroleo  = !is.na(pct_petroleo)  & pct_petroleo  > 0,
    tiene_mineria   = !is.na(pct_mineria)   & pct_mineria   > 0,
    tiene_vias      = !is.na(pct_vias)      & pct_vias      > 0,
    tiene_cultivos  = !is.na(pct_cultivos)  & pct_cultivos  > 0,
    tiene_ganaderia = !is.na(pct_ganaderia) & pct_ganaderia > 0,
    tiene_urbano    = !is.na(pct_urbano)    & pct_urbano    > 0,
    tiene_hab       = cod_dism_habitat == "YES"
  ) %>%
  rowwise() %>%
  mutate(
    codigos_lista = list(c(
      if (tiene_urbano)    "1.1",
      if (tiene_cultivos)  "2.1.4",
      if (tiene_ganaderia) "2.3.4",
      if (tiene_petroleo)  "3.1",
      if (tiene_mineria)   "3.2",
      if (tiene_vias)      "4.1"
    )),
    alcances_lista = list(c(
      if (tiene_urbano)    alcance_sis(pct_urbano),
      if (tiene_cultivos)  alcance_sis(pct_cultivos),
      if (tiene_ganaderia) alcance_sis(pct_ganaderia),
      if (tiene_petroleo)  alcance_sis(pct_petroleo),
      if (tiene_mineria)   alcance_sis(pct_mineria),
      if (tiene_vias)      alcance_sis(pct_vias)
    )),

    cod_amenazas  = if (length(codigos_lista) == 0) "N/A" else
                      paste(codigos_lista,  collapse = ", "),
    cod_tiempo    = if (length(codigos_lista) == 0) "N/A" else
                      paste(rep("Ongoing", length(codigos_lista)), collapse = ", "),
    cod_alcance   = if (length(alcances_lista) == 0) "N/A" else
                      paste(alcances_lista, collapse = ", "),
    cod_severidad = if (length(codigos_lista) == 0) "N/A" else
                      paste(rep("Very rapid declines", length(codigos_lista)), collapse = ", "),
    cod_presiones = if (length(codigos_lista) == 0) "N/A" else
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
  
  p <- c("una","dos","tres","cuatro","cinco",
         "seis","siete","ocho","nueve","diez")
  
  res <- ifelse(
    is.na(n) | n == 0,
    "ninguna",
    ifelse(
      n >= 1 & n <= 10,
      p[n],
      as.character(n)
    )
  )
  
  return(res)
}

amenazas_sp <- amenazas_sp %>%
  mutate(mpios = sapply(tax, mpio_dpto_sp)) %>%
  mutate(desc_amenazas = case_when(
    !tiene_hab & !tiene_petroleo & !tiene_mineria & !tiene_vias &
      !tiene_cultivos & !tiene_ganaderia & !tiene_urbano ~
      "Las subpoblaciones conocidas de la especie se encuentran en hábitats poco perturbados por actividades humanas.",
    TRUE ~ paste0(
      tools::toTitleCase(num_palabras(subpop)), " subpoblacion",
      ifelse(subpop == 1, " conocida", "es conocidas"),
      " de la especie se encuentran en sitios con destrucción y degradación de su hábitat. ",
      "Las alteraciones del hábitat se deben principalmente a",
      ifelse(tiene_urbano,    " la expansión urbana,", ""),
      ifelse(tiene_cultivos,  " la agricultura (cultivos transitorios y permanentes),", ""),
      ifelse(tiene_ganaderia, " la ganadería,", ""),
      ifelse(tiene_petroleo,  " la extracción de hidrocarburos (petróleo y/o gas),", ""),
      ifelse(tiene_mineria,   " la minería de metales,", ""),
      ifelse(tiene_vias,      " la construcción de infraestructura vial,", ""),
      " en los municipios de ", mpios, "."
    )
  ))

dir.create("resultados/amenazas", recursive = TRUE, showWarnings = FALSE)
amenazas_sp <- amenazas_sp %>%
  dplyr::select(-codigos_lista, -alcances_lista)
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
