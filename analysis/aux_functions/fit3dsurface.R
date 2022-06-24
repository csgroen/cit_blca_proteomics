.fit3dsurface <- function(data, x, y, z, dim_size = 100) {
    #-- Get points as spacial coordinates
    coords <- data[,c(x,y)]
    s <- SpatialPointsDataFrame(coords, data)
    
    #-- Make thin-plate spline
    tps <- fields::Tps(coordinates(s), data[,z])
    raster_p   <- raster::raster(s, nrow = dim_size, ncol = dim_size) %>% 
        raster::interpolate(tps)
    
    #-- Return coordinates
    p_df <- raster::as.data.frame(raster_p, xy = TRUE)
    colnames(p_df) <- c(x, y, z)
    return(p_df)
}