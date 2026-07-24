################################################################################
# Evaluación IUCN - Serpocaulon spp.
# Script 02b
# Genera la asignación de cada registro a una subpoblación
################################################################################

library(sf)
library(dplyr)
library(igraph)

dir.create("resultados/ConR/subpoblaciones",
           recursive = TRUE,
           showWarnings = FALSE)

################################################################################
# Leer registros
################################################################################

registros <- read.csv(
  "datos/registros/registros_limpios.csv",
  encoding = "UTF-8"
)

registros <- registros %>%
  filter(!is.na(ddlat),
         !is.na(ddlon))

################################################################################
# Convertir a sf
################################################################################

puntos <- st_as_sf(
  registros,
  coords = c("ddlon","ddlat"),
  crs = 4326
)

################################################################################
# CRS métrico para medir 5 km
################################################################################

puntos <- st_transform(puntos, 9377)

################################################################################
# Función para una especie
################################################################################

asignar_subpoblaciones <- function(datos_sp,
                                   distancia = 5000){
  
  n <- nrow(datos_sp)
  
  if(n == 1){
    
    datos_sp$subpop_id <- 1
    
    return(datos_sp)
    
  }
  
  dist <- st_distance(datos_sp)
  
  vecinos <- which(
    dist <= distancia,
    arr.ind = TRUE
  )
  
  vecinos <- vecinos[
    vecinos[,1] < vecinos[,2],
  ]
  
  g <- graph_from_edgelist(
    as.matrix(vecinos),
    directed = FALSE
  )
  
  if(vcount(g) < n){
    
    g <- add_vertices(
      g,
      n - vcount(g)
    )
    
  }
  
  comp <- components(g)
  
  datos_sp$subpop_id <- comp$membership
  
  datos_sp
  
}

################################################################################
# Aplicar por especie
################################################################################

subpoblaciones_registros <-
  
  puntos %>%
  
  group_split(tax) %>%
  
  lapply(asignar_subpoblaciones) %>%
  
  bind_rows()

################################################################################
# Resumen
################################################################################

resumen <- subpoblaciones_registros %>%
  st_drop_geometry() %>%
  group_by(tax) %>%
  summarise(
    subpop = n_distinct(subpop_id),
    .groups="drop"
  )

################################################################################
# Guardar
################################################################################

write.csv(
  st_drop_geometry(subpoblaciones_registros),
  "resultados/ConR/subpoblaciones/subpoblaciones_registros.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  resumen,
  "resultados/ConR/subpoblaciones/subpoblaciones_recalculadas.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("\n-------------------------------------\n")
cat("Proceso terminado\n")
cat("-------------------------------------\n")

cat("\nNúmero de especies:",
    nrow(resumen))

cat("\nSubpoblaciones totales:",
    sum(resumen$subpop))

cat("\n\nArchivo generado:\n")
cat("resultados/ConR/subpoblaciones/subpoblaciones_registros.csv\n")