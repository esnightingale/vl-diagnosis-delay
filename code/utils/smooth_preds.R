#' Kernel weighted smoothing with arbitrary bounding area
#'
#' @param df sf object (points)
#' @param field weigth field in sf
#' @param bandwith kernel bandwidth (map units)
#' @param resolution output grid resolution (map units)
#' @param zone sf study zone (polygon)
#' @param out_crs EPSG (should be an equal-area projection)
#'
#' @return a raster object
#' @import btb, raster, fasterize, dplyr, plyr, sf
lissage <- function(df, field, bandwidth, resolution, zone, out_crs = 3035) {
  
  if (st_crs(zone)$epsg != out_crs) {
    message("reprojecting data...")
    zone <- st_transform(zone, out_crs)
  }
  
  if (st_crs(df)$epsg != out_crs) {
    message("reprojecting study zone...")
    df <- st_transform(df, out_crs)
  }
  
  zone_bbox <- st_bbox(zone)
  
  # grid generation
  message("generating reference grid...")
  zone_xy <- zone %>% 
    dplyr::select(geometry) %>% 
    st_make_grid(cellsize = resolution,
                 offset = c(plyr::round_any(zone_bbox[1] - bandwidth, resolution, f = floor),
                            plyr::round_any(zone_bbox[2] - bandwidth, resolution, f = floor)),
                 what = "centers") %>%
    st_sf() %>% 
    st_join(zone, join = st_intersects, left = FALSE) %>% 
    st_coordinates() %>% 
    as_tibble() %>% 
    dplyr::select(x = X, y = Y)
  
  # kernel
  message("computing kernel...")
  kernel <- df %>% 
    cbind(., st_coordinates(.)) %>%
    st_set_geometry(NULL) %>% 
    dplyr::select(x = X, y = Y, field) %>% 
    btb::kernelSmoothing(dfObservations = .,
                         sEPSG = out_crs,
                         iCellSize = resolution,
                         iBandwidth = bandwidth,
                         vQuantiles = NULL,
                         dfCentroids = zone_xy)
  
  # rasterization
  message("\nrasterizing...")
  raster::raster(xmn = plyr::round_any(zone_bbox[1] - bandwidth, resolution, f = floor),
                 ymn = plyr::round_any(zone_bbox[2] - bandwidth, resolution, f = floor),
                 xmx = plyr::round_any(zone_bbox[3] + bandwidth, resolution, f = ceiling),
                 ymx = plyr::round_any(zone_bbox[4] + bandwidth, resolution, f = ceiling),
                 resolution = resolution) %>% 
    fasterize::fasterize(kernel, ., field = field)
}

# ---------------------------------------------------------------------------- #
# library(raster)
# 
# st_geometry(boundary) <- st_geometry(boundary)*1000
# st_geometry(dat) <- st_geometry(dat)*1000
# 
# boundary <- boundary %>%
#   st_set_crs(7759)
# dat <- dat %>%
#   st_set_crs(7759)
# 
# dat %>%
#   mutate(pred = pred2$Mean,
#          field = 1) -> dat
#   lissage("field", 20000, 2000, boundary, 7759) -> test
# raster::writeRaster(test, "test.tif")
# 
# test <- raster('test.tif',sep="")
# raster::plot(test, box = FALSE)
# 
# #convert the raster to points for plotting
# test.p <- rasterToPoints(test)
# 
# #Make the points a dataframe for ggplot
# df <- data.frame(test.p)
# #Make appropriate column headings
# colnames(df) <- c("longitude", "latitude", "value")
# 
# ggplot() +
#   geom_raster(data=df, aes(y=latitude, x=longitude, fill=value)) +
#   geom_sf(data=boundary, fill = NA) + 
#   scale_fill_viridis_c(option = "viridis", direction = -1, end = 0.9)
