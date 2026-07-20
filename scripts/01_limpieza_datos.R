# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 01: Limpieza de datos y preparación para ConR y eecorisk
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026

library(readxl)
library(dplyr)

# Cargar datos originales
raw <- read_excel("datos/registros/Copia_Puntos_georreferenciados.xlsx")

# Seleccionar y renombrar columnas de interés
# ddlat y ddlon son los nombres que espera ConR
registros <- raw %>%
  select(
    tax          = Species,
    ddlat        = `Latitud_corregida (Y)`,
    ddlon        = `Longitud_corregida (X)`,
    elev_msnm    = `Elevación (msnm)`,
    colector     = Collector,
    herbario     = Herbaria,
    localidad    = Locality,
    pais         = country,
    departamento = department,
    municipio    = municipality,
    fecha        = Date,
    habito       = Habit,
    id           = ID
  ) %>%
  # Excluir híbridos
  filter(!tax %in% c("Serpocaulon x manizalense", "Serpocaulon x semipinnatifidum")) %>%
  # Excluir registros sin coordenadas
  filter(!is.na(ddlat), !is.na(ddlon)) %>%
  mutate(
    ddlat     = as.numeric(ddlat),
    ddlon     = as.numeric(ddlon),
    elev_msnm = as.numeric(elev_msnm)
  )

# Registros por especie
count(registros, tax, name = "n_registros") %>% arrange(tax)

# Guardar CSV limpio para usar en los siguientes scripts
write.csv(registros,
          "datos/registros/registros_limpios.csv",
          row.names = FALSE,
          fileEncoding = "UTF-8")

# corregir registros
library(dplyr)

registros <- read.csv(
  "datos/registros/registros_limpios.csv",
  encoding = "UTF-8"
)

registros <- registros %>%
  mutate(
    ddlon = if_else(
      pais == "Colombia" & ddlon > 0,
      -ddlon,
      ddlon
    )
  )

write.csv(
  registros,
  "datos/registros/registros_limpios.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
