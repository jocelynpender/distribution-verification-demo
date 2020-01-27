

extract_clim <- function(gbif_data, clim) {
  bioclim_data <- extract(clim, gbif_data, method = "simple", buffer = NULL, small = FALSE, 
    cellnumbers = FALSE, fun = NULL, na.rm = TRUE)
  return(bioclim_data)
}

run_extract_clim <- function(gbif_data, clim) {
  clim_data <- vector("list", length(gbif_data)) # init list with length of num of spp 
  for (i in 1:length(gbif_data)) {
    print(names(gbif_data)[[i]])
    clim_spp <- extract_clim(gbif_data[[i]], clim)
    clim_spp_data <- cbind(gbif_data[[i]], clim_spp)
    clim_data[[i]] <- clim_spp_data
  }
  return(clim_data)
}