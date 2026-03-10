library(fields)

simCov <- function(x_seq,y_seq){
  grid_list <- list(x = x_seq, y = y_seq)
  
  setup <- fields::matern.image.cov(
    setup = TRUE, 
    grid = list(x = x_seq, y = y_seq), 
    theta = 50,       # identical to scale = 50
    smoothness = 0.6, # identical to nu = 0.6
    M = 1024, N = 1024 
  )
  
  sim_val <- sqrt(0.1) * fields::sim.rf(setup) # identical to var = 0.1
  
  
  spatialCov <- list(cov = raster::flip(
    raster::raster(
      t(sim_val), 
      xmn = min(x_seq), xmx = max(x_seq),
      ymn = min(y_seq), ymx = max(y_seq)
    ),
    direction = 'y' 
  ))
  
  return(spatialCov)
}
