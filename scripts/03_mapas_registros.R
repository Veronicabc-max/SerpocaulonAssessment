# Evaluación IUCN - Serpocaulon spp. Colombia
# Script 03: Mapas de registros por especie (verificación y reporte)
# Autora: Verónica Bedoya; Maria Judith Carmona | 2026
# Correr DESPUÉS del script 02 (necesita resultados ConR) y ANTES del script 04

library(sf)
library(ggplot2)
library(dplyr)
library(geodata)
library(ggspatial)
library(cowplot)
library(readxl)
library(mapview)

dir.create("resultados/mapas/por_especie",
           recursive = TRUE,
           showWarnings = FALSE)

registros <- read.csv(
  "datos/registros/registros_limpios.csv",
  encoding = "UTF-8"
) |>
  filter(!is.na(ddlat), !is.na(ddlon))

puntos_sf <- st_as_sf(
  registros,
  coords = c("ddlon","ddlat"),
  crs = 4326
)

colombia <- st_read("datos/capas/pais/Colombia.gpkg", quiet = TRUE)
deptos   <- st_read("datos/capas/pais/Colombia_deptos.gpkg", quiet = TRUE)

eoo_res  <- read.csv("resultados/ConR/EOO/EOO_resultados.csv")
aoo_res  <- read.csv("resultados/ConR/AOO/AOO_resultados.csv")
loc_res  <- read.csv("resultados/ConR/criterioB/localidades.csv")
subp_res <- read.csv("resultados/ConR/subpoblaciones/subpoblaciones.csv")

pct_ap <- read.csv(
  "SIS_Connect/base_maestra.csv",
  check.names = FALSE
) |>
  dplyr::select(
    tax = `NOMBRE CIENTÍFICO sin autor`,
    pct_ap = `% OCURRENCIAS EN AREAS PROTEGIDAS`
  )

params_conr <- eoo_res |>
  rename(EOO_km2 = eoo) |>
  left_join(rename(aoo_res, AOO_km2 = aoo), by="tax") |>
  left_join(rename(loc_res, n_loc = locations), by="tax") |>
  left_join(rename(subp_res, n_subpop = subpop), by="tax") |>
  left_join(pct_ap, by="tax")

eoo_shp <- st_read(
  "resultados/ConR/formas_EOO/EOO_poly.shp",
  quiet = TRUE
)

mapa_especie <- function(sp, guardar=TRUE){
  
  pts_sp <- puntos_sf |> filter(tax==sp)
  if(nrow(pts_sp)==0) return(invisible(NULL))
  
  bbox <- st_bbox(pts_sp)
  
  buf <- 0.8
  
  dx <- max(as.numeric(bbox["xmax"]-bbox["xmin"]),1)
  dy <- max(as.numeric(bbox["ymax"]-bbox["ymin"]),1)
  
  xmid <- mean(c(bbox["xmin"],bbox["xmax"]))
  ymid <- mean(c(bbox["ymin"],bbox["ymax"]))
  
  xlim <- c(xmid-dx/2-buf,xmid+dx/2+buf)
  ylim <- c(ymid-dy/2-buf,ymid+dy/2+buf)
  
  par_sp <- params_conr |> filter(tax==sp)
  
  valor <- function(x){
    if(length(x)==0 || is.na(x)) "N/D" else x
  }
  
  texto <- paste0(
    "EOO = ", valor(round(par_sp$EOO_km2,0)), " km²\n",
    "AOO (2 km) = ", valor(round(par_sp$AOO_km2,0)), " km²\n",
    "Número de registros = ", nrow(pts_sp), "\n",
    "Subpoblaciones (5 km) = ", valor(par_sp$n_subpop), "\n",
    "Localidades (10 km) = ", valor(par_sp$n_loc), "\n",
    "% registros en APs = ", valor(par_sp$pct_ap)
  )
  
  eoo_sp <- eoo_shp |> filter(tax==sp)
  
  p_main <- ggplot() +
    geom_sf(data=deptos,
            fill="grey96",
            color="grey75",
            linewidth=0.2) +
    geom_sf(data=colombia,
            fill=NA,
            color="grey30",
            linewidth=0.5)
  
  if(nrow(eoo_sp)>0 && all(!st_is_empty(eoo_sp))){
    p_main <- p_main +
      geom_sf(
        data=eoo_sp,
        fill="#6baed6",
        alpha=0.25,
        color="#2171b5",
        linewidth=0.7
      )
  }
  
  p_main <- p_main +
    geom_sf(
      data=pts_sp,
      shape=21,
      fill="#d73027",
      color="black",
      size=2,
      stroke=0.4
    )+
    coord_sf(
      xlim=xlim,
      ylim=ylim,
      expand=FALSE
    )+
    annotation_scale(
      location="bl",
      width_hint=0.30
    )+
    annotate(
      "label",
      x=xlim[1]+0.1,
      y=ylim[2]-0.1,
      label=texto,
      hjust=0,
      vjust=1,
      size=2.8,
      linewidth=0.3,
      fill="white",
      alpha=0.85,
      family="mono"
    )+
    labs(title=sp)+
    theme_minimal()+
    theme(
      plot.title=element_text(face="italic"),
      panel.border=element_rect(
        colour="grey50",
        fill=NA
      )
    )
  
  p_inset <- ggplot() +
    geom_sf(
      data=colombia,
      fill="grey90",
      color="grey50",
      linewidth=0.3
    )
  
  if(nrow(eoo_sp)>0 && all(!st_is_empty(eoo_sp))){
    p_inset <- p_inset +
      geom_sf(
        data=eoo_sp,
        fill="#6baed6",
        alpha=0.25,
        color="#2171b5",
        linewidth=0.4
      )
  }
  
  p_inset <- p_inset +
    geom_sf(
      data=pts_sp,
      color="#d73027",
      size=0.8
    )+
    theme_void()+
    theme(
      panel.border=element_rect(
        colour="grey50",
        fill=NA
      )
    )
  
  p_final <- ggdraw(p_main) +
    draw_plot(
      p_inset,
      x=0.72,
      y=0.02,
      width=0.26,
      height=0.28
    )
  
  if(guardar){
    ggsave(
      file.path(
        "resultados/mapas/por_especie",
        paste0(gsub(" ","_",sp),".png")
      ),
      p_final,
      width=8,
      height=7,
      dpi=150
    )
  }
  
  p_final
}

mapview(
  puntos_sf,
  zcol="tax",
  layer.name="Especie",
  map.types="Esri.WorldShadedRelief"
)

sp_prueba <- sort(unique(registros$tax))[1]
print(mapa_especie(sp_prueba,FALSE))

especies <- sort(unique(registros$tax))

pb <- txtProgressBar(
  min=0,
  max=length(especies),
  style=3
)

for(i in seq_along(especies)){
  setTxtProgressBar(pb,i)
  mapa_especie(especies[i],TRUE)
}

close(pb)

message("Mapas guardados en resultados/mapas/por_especie/")
