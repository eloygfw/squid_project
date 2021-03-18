library(sf)
library(ggplot2)
library(dplyr)


readHDF5 <- function(files) {
  
  rad <-
    purrr::map(files, function(file) {
      rhdf5::h5read(file, name = "/All_Data/VIIRS-DNB-SDR_All/Radiance")
    })
  rad_matrix <- do.call("cbind", rad)
  
  lon <-
    purrr::map(files, function(file) {
      rhdf5::h5read(file, name = "/All_Data/VIIRS-DNB-GEO_All/Longitude")
    })
  lon_matrix <- do.call("cbind", lon)
  
  lat <-
    purrr::map(files, function(file) {
      rhdf5::h5read(file, name = "/All_Data/VIIRS-DNB-GEO_All/Latitude")
    })
  lat_matrix <- do.call("cbind", lat)
  
  out <-
    list(rad_matrix = rad_matrix,
         lon_matrix = lon_matrix,
         lat_matrix = lat_matrix)
  
  return(out)
}


readHDF5_df <- function(files) {
  
  list_of_matrix <- readHDF5(files)
  
  df <- data.frame(
    rad = as.vector(list_of_matrix$rad_matrix),
    lon = as.vector(list_of_matrix$lon_matrix),
    lat = as.vector(list_of_matrix$lat_matrix)
  )
  
  return(df)
}


make_raster <- function(df, bin = 100, fun = base::mean){
  
  df_raster <-
    df %>%
    dplyr::mutate(lon_bin = round(bin * lon), lat_bin = round(bin * lat)) %>%
    dplyr::group_by(lon_bin, lat_bin) %>%
    dplyr::summarise(rad = mean(rad)) %>%
    dplyr::mutate(lon = lon_bin * (1.0/bin), lat = lat_bin * (1.0/bin)) %>% 
    dplyr::select(lon, lat, rad)
  
  return(df_raster)
}


plotDNBonly <- function(
  dnb_raster,
  land_sf = NULL,
  eez_sf = NULL,
  file = NULL,
  lower = 0.1, # unit of nano watt
  upper = 10   # unit of nano watt
){
  
  
  
  # cropping data with rad value
  lower <- lower * 1e-9
  upper <- upper * 1e-9
  dnb_raster <-
    dnb_raster %>%
    dplyr::mutate(rad = ifelse(rad < lower, lower, rad)) %>%
    dplyr::mutate(rad = ifelse(rad > upper, upper, rad))
  
  
  
  # plot
  p <-
    ggplot2::ggplot() +
    
    # DNB
    ggplot2::geom_raster(data = dnb_raster,  aes(lon, lat, fill = rad)) +
    ggplot2::scale_colour_gradient(low = "black",
                                   high = "white",
                                   aesthetics = "fill")
  
  if(!is.null(eez_sf)){
    p <- 
      p +
      # plot eez
      ggplot2::geom_sf(
        data = eez_sf,
        color = "yellow",
        fill = NA,
        size = 0.05
      )
  }
  
  if(!is.null(land_sf)){
    p <- 
      p +
      # plot eez
      ggplot2::geom_sf(
        data = land_sf,
        color = "yellow",
        fill = NA,
        size = 0.05
      )
  }
  
  p <- p +
    # restrict plot area to the dnb data
    ggplot2::coord_sf(xlim = c(min(dnb_raster$lon),
                               max(dnb_raster$lon)),
                      ylim = c(min(dnb_raster$lat),
                               max(dnb_raster$lat))) +
    ggplot2::theme_bw()
  
  
  
  if(!is.null(file)){
    # save plot as pdf
    ggplot2::ggsave(
      filename = file,
      plot = p,
      #path = "plot",
      units = "cm",
      width = 20,
      height = 20
    )
  }
  
  
  return(p)
}

######################################################

source("./code/eloy_functions.R")


file1 <- "./input/GDNBO-SVDNB_npp_d20190828_t0701573_e0707377_b40589_c20200214081403538036_nobc_ops.h5"

data <- readHDF5_df(file1)

eez_sf <-
  read_sf('./input/eez_200NM/') %>%
  filter(Territory1 %in% c("Peru", "Ecuador", "Chile", "Colombia", "Panama"))

land_sf <- 
  fishwatchr::land_sf %>% 
  filter(ID %in% c("Peru", "Ecuador", "Chile", "Colombia", "Panama"))


data_raster <- make_raster(data)

p <- plotDNBonly(
  file="./output/peru_20190828.pdf",
  dnb_raster = data_raster,
  eez_sf = eez_sf,
  land_sf = land_sf,
  lower = 0.1,
  upper = 10)



# geotiff -------


library(raster)

#data <- readHDF5_df(c(file1)) %>% dplyr::select(lat, lon, rad) %>% filter(rad>0)

data <-
  readHDF5_df(c(file1)) %>%
  dplyr::select(lat, lon, rad) %>%
  dplyr::mutate(rad = ifelse(rad > 0.0, rad, 0.0)) # replace minus value (detection error) into zero.


test <-
  raster::rasterize(
    x = data[c("lon", "lat")],
    y = raster::raster(ncols = 32512, nrows = 24576),
    field = data$rad,
    fun = mean,
    filename = "./output/peru_20190828_32512_24576.tif",
    overwrite=TRUE
  )
