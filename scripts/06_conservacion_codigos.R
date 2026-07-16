# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 06: Acciones de conservación, investigación y campos finales SIS
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Correr después del script 05

library(dplyr)
library(sf)

registros  <- read.csv("datos/registros/registros_limpios.csv", encoding = "UTF-8")
subpop_res <- read.csv("resultados/ConR/subpoblaciones/subpoblaciones.csv")
base_eeco  <- read.csv("resultados/eecorisk/fragmentacion_severa/resultados_eecorisk.csv",
                       encoding = "UTF-8")

base_maestra <- read.csv("SIS_Connect/base_maestra.csv",
                         encoding = "UTF-8", check.names = FALSE)
names(base_maestra) <- make.unique(names(base_maestra))

# Figuras de conservación: fuerte vs débil
# Figuras fuertes: PNN, SFF, Santuario de Flora, Via Parque, Reserva Natural
# Figuras débiles: DMI, DRMI, Reserva Forestal, Area de Conservacion de Suelos
figuras_fuertes <- c("Parque Nacional Natural", "Parque Natural Regional",
                     "Santuario de Fauna y Flora", "Santuario de Flora",
                     "Vía Parque", "Reserva Natural de la Sociedad Civil")
figuras_debiles <- c("Distrito de Manejo Integrado", "Distrito Regional de Manejo Integrado",
                     "Reserva Forestal Protectora", "Reserva Forestal",
                     "Area de Conservacion de Suelos", "Área de Conservación de Suelos")

tipo_ap <- function(listado) {
  if (is.na(listado) || listado == "") return("ninguna")
  tiene_fuerte <- any(sapply(figuras_fuertes, function(f) grepl(f, listado, ignore.case = TRUE)))
  tiene_debil  <- any(sapply(figuras_debiles, function(f) grepl(f, listado, ignore.case = TRUE)))
  if (tiene_fuerte) return("fuerte")
  if (tiene_debil)  return("debil")
  return("ninguna")
}

num_palabras <- function(n) {
  if (is.na(n) || n == 0) return("ninguna")
  p <- c("una","dos","tres","cuatro","cinco","seis","siete","ocho","nueve","diez")
  if (n >= 1 && n <= 10) p[n] else as.character(n)
}

# Construir campos de conservación
conservacion <- base_maestra %>%
  dplyr::select(
    tax        = `NOMBRE CIENTÍFICO sin autor`,
    pct_ap     = `% OCURRENCIAS EN AREAS PROTEGIDAS`,
    listado_ap = `LISTADO DE AREAS PROTEGIDAS CON OCURRENCIAS`,
    dism_hab   = `CÓDIGO SIS DISMINUCIÓN CONTINUA HÁBITAT`
  ) %>%
  left_join(subpop_res, by = "tax") %>%
  left_join(base_eeco %>% dplyr::select(tax, pct_HH), by = "tax") %>%
  mutate(
    amenazada    = dism_hab == "YES",
    tipo_ap_sp   = sapply(listado_ap, tipo_ap),
    tiene_ap     = !is.na(pct_ap) & pct_ap > 0,

    # Descripción acciones de conservación
    desc_conservacion = case_when(
      !amenazada & !tiene_ap ~
        "La especie no está presente en áreas del \"Sistema Nacional de Áreas Protegidas\". No se conocen acciones de conservación para la especie.",
      !amenazada & tipo_ap_sp == "fuerte" ~
        "La especie está presente en áreas del \"Sistema Nacional de Áreas Protegidas\". No se conocen otras acciones de conservación para la especie.",
      !amenazada & tipo_ap_sp == "debil" ~
        "La especie está presente en áreas del \"Sistema Nacional de Áreas Protegidas\"; sin embargo, las subpoblaciones se encuentran en áreas bajo figuras legales que son poco efectivas para la conservación de la biodiversidad. No se conocen otras acciones de conservación para la especie.",
      amenazada & !tiene_ap ~
        "La especie no está presente en áreas del \"Sistema Nacional de Áreas Protegidas\". No se conocen acciones de conservación para la especie.",
      amenazada & tipo_ap_sp == "fuerte" ~
        paste0(tools::toTitleCase(num_palabras(subpop)),
               " subpoblacion", ifelse(subpop == 1, " conocida está presente", "es conocidas están presentes"),
               " en áreas del \"Sistema Nacional de Áreas Protegidas\". No se conocen otras acciones de conservación para la especie."),
      amenazada & tipo_ap_sp == "debil" ~
        paste0(tools::toTitleCase(num_palabras(subpop)),
               " subpoblacion", ifelse(subpop == 1, " conocida está presente", "es conocidas están presentes"),
               " en áreas del \"Sistema Nacional de Áreas Protegidas\"; sin embargo, las subpoblaciones se encuentran en áreas bajo figuras legales que son poco efectivas para la conservación de la biodiversidad. No se conocen otras acciones de conservación para la especie."),
      TRUE ~ NA_character_
    ),

    # Código conservación requerida
    # 1.1 siempre si amenazada y solo APs débiles; 3.4.2 y 4.3 siempre si amenazada
    # 2.3 si HH > 90%; 3.3.1 si subpobs < 5 (muy pocas)
    cod_conservacion = case_when(
      !amenazada ~ "N/A",
      amenazada & tipo_ap_sp != "fuerte" & !is.na(pct_HH) & pct_HH > 90 ~
        "1.1, 2.3, 3.4.2, 4.3",
      amenazada & tipo_ap_sp != "fuerte" ~
        "1.1, 3.4.2, 4.3",
      amenazada & !is.na(pct_HH) & pct_HH > 90 ~
        "2.3, 3.4.2, 4.3",
      TRUE ~ "3.4.2, 4.3"
    ),

    # Código investigación requerida
    cod_investigacion = case_when(
      !amenazada ~ "N/A",
      TRUE       ~ "1.2, 1.3, 1.5, 1.6, 2.2, 3.1, 3.4"
    ),

    # Exsitu: Unknown por defecto (requiere consulta manual a jardines botánicos)
    reporte_exsitu = "Unknown",
    listado_exsitu = "N/A"
  )

base_maestra <- base_maestra %>%
  left_join(dplyr::select(conservacion, tax, desc_conservacion, cod_conservacion,
                           cod_investigacion, reporte_exsitu, listado_exsitu),
            by = c("NOMBRE CIENTÍFICO sin autor" = "tax")) %>%
  mutate(
    `DESCRIPCIÓN ACCIONES DE CONSERVACIÓN`  = coalesce(`DESCRIPCIÓN ACCIONES DE CONSERVACIÓN`, desc_conservacion),
    `CÓDIGO SIS CONSERVACIÓN REQUERIDA`     = coalesce(`CÓDIGO SIS CONSERVACIÓN REQUERIDA`, cod_conservacion),
    `CÓDIGO SIS INVESTIGACION REQUERIDA`    = coalesce(`CÓDIGO SIS INVESTIGACION REQUERIDA`, cod_investigacion),
    `REPORTE DE CONSERVACION EXSITU`        = coalesce(`REPORTE DE CONSERVACION EXSITU`, reporte_exsitu),
    `LISTADO DE COLECCIONES EXSITU CON INDIVIDUOS` = coalesce(
      `LISTADO DE COLECCIONES EXSITU CON INDIVIDUOS`, listado_exsitu)
  ) %>%
  dplyr::select(-desc_conservacion, -cod_conservacion, -cod_investigacion,
                -reporte_exsitu, -listado_exsitu)

write.csv(base_maestra, "SIS_Connect/base_maestra.csv",
          row.names = FALSE, fileEncoding = "UTF-8")

message("Script 06 completado. Base maestra actualizada con todos los campos disponibles.")
message("Campos pendientes de revisión manual:")
message("  - DESCRIPCIÓN TAMAÑO POBLACIONAL Y DEMOGRAFÍA: verificar si hay literatura")
message("  - LISTADO DE COLECCIONES EXSITU: consultar herbarios y jardines botánicos")
message("  - DESCRIPCIÓN AMENAZAS: revisar y ajustar texto generado automáticamente")
message("  - Agregar amenazas de cobertura de tierra cuando esté disponible la capa (script 05)")
