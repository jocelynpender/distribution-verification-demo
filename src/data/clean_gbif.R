require(dplyr)
require(ggplot2)
require(rgbif)
require(sp)
require(countrycode)
require(CoordinateCleaner)
require(ggplot2)
require(BIEN)
require(maps)  #Useful for making quick maps of occurrences
require(gtools)
require(maps)
require(sp)
require(maptools)

download_gbif <- function(sp_name) {
  print(sp_name)
  dat <- occ_search(scientificName = sp_name, limit = 200000, return = "data", hasCoordinate = T)  # obtain data from GBIF via rgbif
  if (is.character(dat) == FALSE) { # GBIF will return a character if there is no data found
    dat <- dat %>% dplyr::select(one_of("species", "decimalLongitude", "decimalLatitude", 
      "countryCode", "gbifID", "family", "taxonRank", "coordinateUncertaintyInMeters", 
      "year", "basisOfRecord", "institutionCode", "datasetName"))  # select columns of interest
    dat <- dat %>% filter(!is.na(decimalLongitude)) %>% filter(!is.na(decimalLatitude))  # remove records without coordinates
    dat$countryCode <- countrycode(dat$countryCode, origin = "iso2c", destination = "iso3c")  # convert country code from ISO2c to ISO3c
    dat <- data.frame(dat)
    return(dat)
  }
}

clean_gbif <- function(dat) {
  print(dat$species %>% unique)
  flag <- try(clean_coordinates(x = dat, lon = "decimalLongitude", lat = "decimalLatitude", 
    countries = "countryCode", species = "species", tests = c("capitals", "centroids", 
      "equal", "gbif", "institutions", "zeros", "duplicates", "seas", "validity")))  # most test are on by default
  return(flag)
}

get_range <- function(dat) {
  sp_name <- dat$species %>% unique
  range <- try(BIEN_ranges_load_species(species = sp_name))
  return(range)
}

get_range_name <- function(sp_name) {
  print(sp_name)
  range <- try(BIEN_ranges_load_species(species = sp_name))
  return(range)
}

clean_ranges <- function(ranges) {
  failed <- c()
  for (i in 1:length(ranges)) {
    if (is.character(ranges[[i]])) {
      failed <- c(failed, i)
    } else {
      ranges[[i]]@data[, ] <- names(ranges)[i]  # Fix BIEN name, returns species name separated by an underscore
    }
  }
  ranges <- ranges[-failed]
  return(ranges)
}

range_comparison <- function(dat, range) {
  range_flags <- try(cc_iucn(x = dat, range = range, lon = "decimalLongitude", lat = "decimalLatitude", 
    value = "flagged", buffer = 1))
  dat_range_flags <- cbind(dat, range_flags)
  return(dat_range_flags)
}

plot_flag <- function(flag, device) {
  sp_name <- flag$species %>% unique
  date <- format(Sys.time(), "%Y-%m-%d")
  file_name <- paste(sp_name, "_map_", date, ".", device, sep = "")
  plot(flag, lon = "decimalLongitude", lat = "decimalLatitude", details = TRUE)
  ggsave(file_name, device = device)
}

plot_range_comparison <- function(dat_range_flags, sp_range) {
  sp_name <- dat_range_flags$species %>% unique
  date <- format(Sys.time(), "%Y-%m-%d")
  file_name <- paste(sp_name, "_range_map_", date, ".pdf", sep = "")
  exclude <- dat_range_flags[dat_range_flags$range_flags == FALSE, ]  # Plot the range test results
  include <- dat_range_flags[dat_range_flags$range_flags == TRUE, ]
  pdf(file = file_name, height = 20, width = 20)
  map("world", fill = TRUE, col = "grey")  #, xlim = c(-180,-20), ylim = c(-60,80))
  plot(sp_range, col = "green", add = TRUE)
  points(exclude$decimalLongitude, exclude$decimalLatitude, col = "red", cex = 0.6)
  points(include$decimalLongitude, include$decimalLatitude, col = "purple", cex = 0.6)
  dev.off()
}

get_poly_area <- function(range) {
  poly <- sf::st_as_sf(range)  # calculate area of polygon using sf package
  poly$area <- st_area(poly)  #Take care of units
  area <- poly$area %>% as.numeric
  return(area)
}

merge_flags <- function(flags, dat_range_flags, cutoff, areas) {
  merged_flags <- list()
  for (i in 1:length(flags)) {
    print(to_merge_flags[[i]]$species %>% unique)
    if (areas[[i]] > cutoff) {
      merged_flags[[i]] <- cbind(flags[[i]], dat_range_flags[[i]]$range_flags)
      # add range comparison flags to flags list as new column
    } else {
      merged_flags[[i]] <- rep(NA, nrow(flags[[i]])) %>% cbind(flags[[i]], .)
    }
    names(merged_flags[[i]])[ncol(merged_flags[[i]])] <- "range_flags"
  }
  return(merged_flags)
}

summarize_flags <- function(flags) {
  summaries <- lapply(flags, summary)
  names(summaries) <- names(flags)
  cc_tests <- lapply(summaries, length) %>% as.numeric
  few_tests <- which(cc_tests < 10)
  bind_first <- summaries[-few_tests]
  bind_last <- summaries[few_tests]
  summary <- bind_rows(bind_first) %>% t # summary function from cc
  summary_few <- bind_rows(bind_last) %>% t 
  #summaries_df <- lapply(summaries, as.data.frame)
  #matched_summaries <- do.call(smartbind, summaries_df)
  colnames(summary) <- names(summaries[[1]])  # Apply summary result names to data frame
  colnames(summary_few) <- summaries[few_tests][[1]] %>% names
  summary_bound <- smartbind(summary, summary_few)
  rownames(summary_bound) <- c(rownames(summary), rownames(summary_few))
  summary_bound <- summary_bound[order(rownames(summary_bound)),]
  return(summary)
}

# function to compare number of records returned after cleaning of first download
# set in 2014 and modern download set
count_clean <- function(dat, clean_dat, clean_range_dat) {
  counts <- matrix(nrow = 5) %>% data.frame
  for (i in 1:length(dat)) {
    before <- nrow(dat[[i]])
    after <- nrow(clean_dat[[i]])
    range <- nrow(clean_range_dat[[i]])
    prop <- (after/before) * 100
    range_prop <- (range/after) * 100
    counts <- rbind(before, after, range, prop, range_prop) %>% cbind(counts, 
      .)
  }
  counts <- counts[, -1]  # Remove initliazed first column
  colnames(counts) <- names(clean_dat)
  return(counts)
}

decimal_places <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    num_places <- nchar(strsplit(sub("0+$", "", as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    tens <- 1 * 10^num_places
    return(tens)
  } else {
    return(0)
  }
}

precise_enough <- function(latitude) {
  calc_error <- decimal_places(latitude)
  precise_enough <- ifelse(calc_error > 10, TRUE, FALSE)
  return(precise_enough)
}

# The function below has been abandoned
run_subset <- function(fun, size, input_names = NULL, input_data = NULL) {
  groups <- split(input_names, ceiling(seq_along(input_names)/size))
  subset_dat <- list()
  dat <- list()
  
  for (i in 1:length(groups)) {
    print(paste("==== running on group", i, "of", length(groups), "groups ===="))
    if (is.null(input_data) == TRUE) {
      names <- unlist(groups[i])
      subset_dat <- lapply(names, fun)
    } else {
      subset <- groups[i] %>% unlist(., use.names = FALSE)
      subset_dat <- lapply(input_data[subset], fun)  # run coordinate cleaner
    }  # download gbif data using occ_search and only return columns we want. Clean up country codes
    dat <- c(dat, subset_dat)
  }
  
  return(dat)
}

check_usa_polygon <- function(data, usa) {
  points <- data[,2:3] # Long and lat in this order
  test_points <- lapply(usa@polygons, function(x) pnt.in.poly(points, 
    x@Polygons[[1]]@coords))
  usa_points <- lapply(test_points, function(x) x[x$pip == TRUE,]) %>% bind_rows
  return(usa_points)
}

# gbifdata[[i]]$calc_error <- ifelse(gbifdata[[i]]$decimalLatitude ==
# as.integer(gbifdata[[i]]$decimalLatitude), 100, ifelse((10 *
# gbifdata[[i]]$decimalLatitude) == as.integer(10 *
# gbifdata[[i]]$decimalLatitude), 10, ifelse((100 *
# gbifdata[[i]]$decimalLatitude) == as.integer(100 *
# gbifdata[[i]]$decimalLatitude), 1, ifelse((1000 *
# gbifdata[[i]]$decimalLatitude) == as.integer(1000 *
# gbifdata[[i]]$decimalLatitude), 0.1, ifelse((10000 *
# gbifdata[[i]]$decimalLatitude) == as.integer(10000 *
# gbifdata[[i]]$decimalLatitude), 0.01, ifelse((1e+05 *
# gbifdata[[i]]$decimalLatitude) == as.integer(1e+05 *
# gbifdata[[i]]$decimalLatitude), 0.001, 1e-04)))))) gbifdata[[i]]$precise_enough
# <- ifelse(gbifdata[[i]]$calc_error < 10, TRUE, FALSE)
