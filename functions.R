veg_scale <- scale_fill_manual(values = 
                                 c("firebrick",
                                   "darkorange1",
                                   "darkgreen", 
                                   "cadetblue", 
                                   "burlywood4"), 
                               #labels = c("Very frequent fire \nforest", 
                               #           "Frequent fire \nforest", 
                               labels=c("Frequent fire forest\n(FRI < 20)",
                                        "Frequent fire forest\n(FRI 20-40)",
                                        "Infrequent fire forest",
                                        "Grass",
                                        "Shrub"),
                               name = "",
                               na.translate=F)

cat3_scale <- scale_fill_manual(#values = 
                                # viridisLite::viridis(n = 5),
                               values = c(
                                 "indianred4",
                                 "indianred1",
                                 "lightyellow",
                                 "cadetblue1",
                                 "cadetblue4"),
                               name = "",
                               na.translate=F)

cat5_scale <- scale_fill_manual(#values = 
                                #  viridisLite::viridis(n = 7),
                                values = c(
                                  "indianred4",
                                  "indianred3",
                                  "indianred1",
                                  "lightyellow",
                                  "cadetblue1",
                                  "cadetblue3",
                                  "cadetblue4"),
                                name = "",
                                na.translate=F)

blank_axes <- theme(axis.title.x=element_blank(),
                    axis.text.x=element_blank(),
                    axis.ticks.x=element_blank(),
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    panel.border = element_blank())

do_a_region <- function(region_name, region_shape, make_plots = F) {
  
  #----- Clip FRI ------------------------------------------------------------#
  fri_crop <- terra::crop(fri     , region_shape)
  fri_crop <- terra::mask(fri_crop, region_shape)
  
  #----- Clip veg ------------------------------------------------------------#
  veg_crop <- terra::crop(veg     , region_shape)
  veg_crop <- terra::mask(veg_crop, region_shape)
  
  #----- Clip WUI ------------------------------------------------------------#
  wui_crop <- terra::crop(wui     , region_shape)
  wui_crop <- terra::mask(wui_crop, region_shape)
  #---------------------------------------------------------------------------#
  
  
  
  
  
  #---------------------------------------------------------------------------#
  # Forests.
  # Frequent fire forests FRI 0-20 years
  # Frequent fire forest FRI 20-40
  # Infrequent fire forests (FRI > 40)
  #---------------------------------------------------------------------------#
  
  #----- ID forests by fire frequency ----------------------------------------#
  forest_mask <- terra::ifel(veg_crop == "Forest", 1, NA)
  forest_fri  <- terra::mask(fri_crop, forest_mask)  
  forest_mat  <- cbind(c(0,20,40), c(20, 40, 100000), 1:3)
  forest_type <- terra::classify(forest_fri, forest_mat)
  levels(forest_type) <- data.frame(ID = 1:3, 
                                    category = c(
                                      "Frequent fire (FRI < 20)", 
                                      "Frequent fire (FRI 20-40)", 
                                      "Infrequent fire"))
  
  
  
  
  
  #---------------------------------------------------------------------------#
  # If plotting: veg map, with 5 vegetation types
  #---------------------------------------------------------------------------#
  if (make_plots) {
    
    all_veg_map <- veg_crop + 2
    all_veg_map <- terra::classify(all_veg_map, cbind(3,0))
    all_veg_map <- terra::classify(all_veg_map, cbind(NA,0))
    ff <- terra::classify(forest_fri, forest_mat)
    ff <- terra::classify(ff, cbind(NA, 0))
    
    all_veg_map <- all_veg_map + ff
    all_veg_map <- terra::classify(all_veg_map, cbind(0, NA))
    levels(all_veg_map) <- data.frame(ID = 1:5,
                                      category = c(
                                        "Frequent fire (FRI < 20)", 
                                        "Frequent fire (FRI 20-40)", 
                                        "Infrequent fire",
                                        "Grass",
                                        "Shrub"))
    
    suppressMessages(veg_map_plot <- ggplot() + 
                       geom_spatraster(data = all_veg_map) +
                       geom_spatvector(data = region_shape, fill=NA) +
                       theme_bw() +
                       blank_axes +
                       ggtitle("Vegetation map") +
                       veg_scale)
    
    rm(ff, all_veg_map)
  }
  rm(forest_fri)
  invisible(gc())
  
  
  
  #----- Make a total fire deficit map, including long deficit ---------------#
  long_crop <- terra::crop(long_forst_def, region_shape)
  long_crop <- terra::mask(long_crop     , region_shape)
  
  short_crop <- terra::crop(short_def , region_shape)
  short_crop <- terra::mask(short_crop, region_shape)
  
  #----- I want to only substitute eco-level deficit where there isn't short -#
  long_crop <- terra::mask(long_crop, short_crop, inverse = T)
  all_def <- terra::subst(long_crop, NA, 0) + terra::subst(short_crop, NA, 0)
  rm(short_crop, long_crop)
  invisible(gc())
  #----- Mask back to FRI, no values if FRI was NA ---------------------------#
  all_def <- terra::mask(all_def, fri_crop)
  all_def <- terra::mask(all_def, forest_type)
  
  if (make_plots) {
    forest_all_def <- all_def
  }
  
  #---------------------------------------------------------------------------#
  # Extract table stats
  #---------------------------------------------------------------------------#
  
  #----- Overall area --------------------------------------------------------#
  vv <- values(forest_mask)
  forest_cells <- sum(vv, na.rm=T)
  rm(vv)
  invisible(gc())
  
  vv <- values(forest_type)
  tab <- table(vv)
  # Make sure all categories are present
  valu <- rep(0, nrow(levels(forest_type)[[1]]))
  valu[as.numeric(names(tab))] <- as.numeric(unlist(tab))
  
  summary_dat <- data.frame(Vegetation = "All forest",
                            Numcells = forest_cells)

  summary_dat <- rbind(summary_dat,
                       data.frame(Vegetation = levels(forest_type)[[1]]$category,
                                  Numcells = valu))
  rm(vv, tab, valu)
  invisible(gc())
  
  #----- Area in surplus -----------------------------------------------------#
  summary_dat$Numsurpluscells <- 0
  
  surplus <- terra::ifel(all_def > 1.2, 1, 0)
  # Surplus - All forest
  vv <- values(surplus)
  surplus_cells <- sum(vv, na.rm=T)
  x <- which(summary_dat$Vegetation == "All forest")
  summary_dat$Numsurpluscells[x] <- surplus_cells
  rm(vv, surplus_cells)
  
  # Surplus - By forest type
  datum <- terra::zonal(surplus, forest_type, fun="sum")
  x <- match(datum$category, summary_dat$Vegetation)
  summary_dat$Numsurpluscells[x] <- datum$forest_def
  
  rm(surplus, datum)
  invisible(gc())
  
  #----- Area in deficit -----------------------------------------------------#
  summary_dat$Numdeficitcells <- 0
  deficit <- terra::ifel(all_def < -1.2, 1, 0)
  # All forest
  vv <- values(deficit)
  deficit_cells <- sum(vv, na.rm=T)
  x <- which(summary_dat$Vegetation == "All forest")
  summary_dat$Numdeficitcells[x] <- deficit_cells
  rm(vv, deficit_cells)
  
  # By forest type
  datum <- terra::zonal(deficit, forest_type, fun="sum")
  x <- match(datum$category, summary_dat$Vegetation)
  summary_dat$Numdeficitcells[x] <- datum$forest_def
  rm(deficit, datum)
  invisible(gc())
  
  #----- Just right - -1.2 to 1.2 --------------------------------------------#
  summary_dat$Numgoldilockscells <- 0
  mat <- matrix(c(-100000000, -1.2, 0,
           -1.2, 1.2, 1,
           1.2, 1000000000, 0), byrow=T, nrow=3)
  goldilocks <- terra::classify(all_def, mat)
  x <- which(summary_dat$Vegetation == "All forest")
  summary_dat$Numgoldilockscells[x] <- sum(values(goldilocks), na.rm=T)
  
  # By forest type
  datum <- terra::zonal(goldilocks, forest_type, fun="sum")
  x <- match(datum$category, summary_dat$Vegetation)
  summary_dat$Numgoldilockscells[x] <- datum$forest_def
  rm(goldilocks)
  invisible(gc())
  
  #----- Area for which we have no data --------------------------------------#
  summary_dat$Nodatacells <- NA
  no_dat <- !is.na(forest_mask) & is.na(all_def)
  x <- which(summary_dat$Vegetation == "All forest")
  summary_dat$Nodatacells[x] <- sum(values(no_dat))
  rm(no_dat)
  invisible(gc())
  
  #----- Mean and sd of deficit ----------------------------------------------#
  summary_dat$mean_deficit <- NA
  summary_dat$sd_deficit   <- NA
  
  # All forest
  vv <- values(all_def)
  x <- which(summary_dat$Vegetation == "All forest")
  summary_dat$mean_deficit[x] <- mean(vv, na.rm=T)
  summary_dat$sd_deficit[x]   <- sd  (vv, na.rm=T)
  rm(vv)
  
  datum <- terra::zonal(all_def, forest_type, fun="mean")
  x <- match(datum$category, summary_dat$Vegetation)
  summary_dat$mean_deficit[x] <- unlist(datum$forest_def)
  
  datum <- terra::zonal(all_def, forest_type, fun="sd")
  x <- match(datum$category, summary_dat$Vegetation)
  summary_dat$sd_deficit[x] <- unlist(datum$forest_def)
  
  #---------------------------------------------------------------------------#
  # Forest WUI 
  #---------------------------------------------------------------------------#
  summary_dat <- rbind(summary_dat, 
                       data.frame(Vegetation = c("Forest WUI", "Forest non-WUI"),
                                  Numcells = 0, Numsurpluscells = 0, 
                                  Numdeficitcells = 0, Numgoldilockscells = NA, 
                                  Nodatacells = NA,
                                  mean_deficit = 0, sd_deficit = 0))

  #----- Overall forest WUI area ---------------------------------------------#
  datum <- mask(wui_crop, forest_type)
  summary_dat$Numcells[summary_dat$Vegetation == "Forest WUI"] <- sum(!is.na(values(datum)), na.rm=T)
  datum <- mask(forest_type, wui_crop, inverse=T)
  summary_dat$Numcells[summary_dat$Vegetation == "Forest non-WUI"] <- sum(!is.na(values(datum)), na.rm=T)
  rm(datum)
  invisible(gc())
  
  #----- Area in surplus and deficit -----------------------------------------#
  surplus <- terra::ifel(all_def > 1.2, 1, NA)
  surp <- terra::mask(surplus, wui_crop)
  summary_dat$Numsurpluscells[summary_dat$Vegetation == "Forest WUI"] <- sum(values(surp), na.rm=T)
  surp <- terra::mask(surplus, wui_crop, inverse = T)
  summary_dat$Numsurpluscells[summary_dat$Vegetation == "Forest non-WUI"] <- sum(values(surp), na.rm=T)
  
  rm(surplus, surp)
  invisible(gc())
  
  deficit <- terra::ifel(all_def < -1.2, 1, NA)
  defr <- terra::mask(deficit, wui_crop)
  summary_dat$Numdeficitcells[summary_dat$Vegetation == "Forest WUI"] <- sum(values(defr), na.rm=T)
  defr <- terra::mask(deficit, wui_crop, inverse = T)
  summary_dat$Numdeficitcells[summary_dat$Vegetation == "Forest non-WUI"] <- sum(values(defr), na.rm=T)
  rm(deficit, defr)
  invisible(gc())
  
  goldilocks <- terra::classify(all_def, mat)
  goldf <- terra::mask(goldilocks, wui_crop)
  summary_dat$Numgoldilockscells[summary_dat$Vegetation == "Forest WUI"] <- sum(values(goldf), na.rm=T)
  goldf <- terra::mask(goldilocks, wui_crop, inverse = T)
  summary_dat$Numgoldilockscells[summary_dat$Vegetation == "Forest non-WUI"] <- sum(values(goldf), na.rm=T)
  rm(goldilocks, goldf)
  invisible(gc())
  
  #----- Mean and sd of deficit ----------------------------------------------#
  defr <- terra::mask(all_def, wui_crop)
  summary_dat$mean_deficit[summary_dat$Vegetation == "Forest WUI"] <- mean(values(defr), na.rm=T)
  summary_dat$sd_deficit  [summary_dat$Vegetation == "Forest WUI"] <- sd(values(defr), na.rm=T)
  
  defr <- terra::mask(all_def, wui_crop, inverse=T)
  summary_dat$mean_deficit[summary_dat$Vegetation == "Forest non-WUI"] <- mean(values(defr), na.rm=T)
  summary_dat$sd_deficit  [summary_dat$Vegetation == "Forest non-WUI"] <- sd(values(defr), na.rm=T)
  
  rm(forest_type, all_def)
  invisible(gc())
  
  
  
  #---------------------------------------------------------------------------#
  # Grass
  #---------------------------------------------------------------------------#
  grass_mask <- terra::ifel(veg_crop == "Grass", 1, NA)

  #----- Make a total fire deficit map, including long deficit ---------------#
  long_crop <- terra::crop(long_grass_def, region_shape)
  long_crop <- terra::mask(long_crop     , region_shape)
  
  short_crop <- terra::crop(short_def , region_shape)
  short_crop <- terra::mask(short_crop, region_shape)

  #----- I want to only substitute eco-level deficit where there isn't short -#
  long_crop <- terra::mask(long_crop, short_crop, inverse = T)
  all_def <- terra::subst(long_crop, NA, 0) + terra::subst(short_crop, NA, 0)
  rm(short_crop, long_crop)
  invisible(gc())
  #----- Mask back to FRI, no values if FRI was NA ---------------------------#
  all_def <- terra::mask(all_def, fri_crop)
  all_def <- terra::mask(all_def, grass_mask)
  
  if (make_plots) {
    grass_all_def <- all_def
  }
  #---------------------------------------------------------------------------#
  
  #----- Extract table stats -------------------------------------------------#
  
  # Overall area
  vv <- values(grass_mask)
  grass_cells <- sum(vv, na.rm=T)
  
  rm(vv)
  invisible(gc())
  
  #----- Area in surplus and deficit -----------------------------------------#
  surplus <- terra::ifel(all_def > 1.2, 1, 0)
  vv <- values(surplus)
  surplus_cells <- sum(vv, na.rm=T)
  rm(surplus, vv)
  
  deficit <- terra::ifel(all_def < -1.2, 1, 0)
  vv <- values(deficit)
  deficit_cells <- sum(vv, na.rm=T)
  rm(deficit, vv)
  invisible(gc())
  
  #----- Just right ----------------------------------------------------------#
  goldilocks <- terra::classify(all_def, mat)
  goldilocks_cells <- sum(values(goldilocks), na.rm=T)
  rm(goldilocks)
  invisible(gc())
  
  #----- Area for which we have no data --------------------------------------#
  no_dat <- !is.na(grass_mask) & is.na(all_def)
  no_data_cells <- sum(values(no_dat))
  rm(no_dat)
  invisible(gc())
  
  #----- Mean and sd of deficit ----------------------------------------------#
  vv <- values(all_def)
  mean_deficit <- mean(vv, na.rm=T)
  sd_deficit   <- sd  (vv, na.rm=T)
  
  summary_dat <- rbind(summary_dat, data.frame(
    Vegetation         = "Grassland",
    Numcells           = grass_cells,
    Numsurpluscells    = surplus_cells,
    Numdeficitcells    = deficit_cells,
    Numgoldilockscells = goldilocks_cells,
    Nodatacells        = no_data_cells,
    mean_deficit       = mean_deficit,
    sd_deficit         = sd_deficit
  ))
  
  #---------------------------------------------------------------------------#
  # Grass WUI 
  #---------------------------------------------------------------------------#
  summary_dat <- rbind(summary_dat, 
                       data.frame(Vegetation = c("Grass WUI", "Grass non-WUI"),
                                  Numcells = 0, Numsurpluscells = 0, 
                                  Numdeficitcells = 0, Numgoldilockscells = NA,
                                  Nodatacells = NA,
                                  mean_deficit = 0, sd_deficit = 0))

  #----- Overall grass WUI area ---------------------------------------------#
  datum <- mask(grass_mask, wui_crop)
  numcells <- sum(!is.na(values(datum)), na.rm=T)
  summary_dat$Numcells[summary_dat$Vegetation == "Grass WUI"] <- numcells
  datum <- mask(grass_mask, wui_crop, inverse=T)
  numcells <- sum(!is.na(values(datum)), na.rm=T)
  summary_dat$Numcells[summary_dat$Vegetation == "Grass non-WUI"] <- numcells
  rm(datum)
  invisible(gc())
  
  #----- Area in surplus and deficit -----------------------------------------#
  surplus <- terra::ifel(all_def > 1.2, 1, NA)
  surp <- terra::mask(surplus, wui_crop)
  numcells <- sum(values(surp), na.rm=T)
  summary_dat$Numsurpluscells[summary_dat$Vegetation == "Grass WUI"] <- numcells
  surp <- terra::mask(surplus, wui_crop, inverse = T)
  numcells <- sum(values(surp), na.rm=T)
  summary_dat$Numsurpluscells[summary_dat$Vegetation == "Grass non-WUI"] <- numcells

  rm(surplus, surp)
  invisible(gc())
  
  deficit <- terra::ifel(all_def < -1.2, 1, NA)
  defr <- terra::mask(deficit, wui_crop)
  numcells <- sum(values(defr), na.rm=T)
  summary_dat$Numdeficitcells[summary_dat$Vegetation == "Grass WUI"] <- numcells
  defr <- terra::mask(deficit, wui_crop, inverse = T)
  numcells <- sum(values(defr), na.rm=T)
  summary_dat$Numdeficitcells[summary_dat$Vegetation == "Grass non-WUI"] <- numcells
  rm(deficit, defr)
  invisible(gc())
  
  goldilocks <- terra::classify(all_def, mat)
  goldf <- terra::mask(goldilocks, wui_crop)
  summary_dat$Numgoldilockscells[summary_dat$Vegetation == "Grass WUI"] <- sum(values(goldf), na.rm=T)
  goldf <- terra::mask(goldilocks, wui_crop, inverse = T)
  summary_dat$Numgoldilockscells[summary_dat$Vegetation == "Grass non-WUI"] <- sum(values(goldf), na.rm=T)
  rm(goldilocks, goldf)
  invisible(gc())
  
  #----- Mean and sd of deficit ----------------------------------------------#
  defr <- terra::mask(all_def, wui_crop)
  summary_dat$mean_deficit[summary_dat$Vegetation == "Grass WUI"] <- mean(values(defr), na.rm=T)
  summary_dat$sd_deficit  [summary_dat$Vegetation == "Grass WUI"] <- sd(values(defr), na.rm=T)
  
  defr <- terra::mask(all_def, wui_crop, inverse=T)
  summary_dat$mean_deficit[summary_dat$Vegetation == "Grass non-WUI"] <- mean(values(defr), na.rm=T)
  summary_dat$sd_deficit  [summary_dat$Vegetation == "Grass non-WUI"] <- sd(values(defr), na.rm=T)
  
  rm(grass_mask, vv, all_def)
  invisible(gc())
  #---------------------------------------------------------------------------#
  
  
  
  #---------------------------------------------------------------------------#
  # Shrub
  #---------------------------------------------------------------------------#
  shrub_mask <- terra::ifel(veg_crop == "Shrub", 1, NA)
  
  #----- Make a total fire deficit map, including long deficit ---------------#
  long_crop <- terra::crop(long_shrub_def, region_shape)
  long_crop <- terra::mask(long_crop     , region_shape)
  
  short_crop <- terra::crop(short_def , region_shape)
  short_crop <- terra::mask(short_crop, region_shape)

  #----- I want to only substitute eco-level deficit where there isn't short -#
  long_crop <- terra::mask(long_crop, short_crop, inverse = T)
  all_def <- terra::subst(long_crop, NA, 0) + terra::subst(short_crop, NA, 0)
  rm(short_crop, long_crop)
  invisible(gc())
  #----- Mask back to FRI, no values if FRI was NA ---------------------------#
  all_def <- terra::mask(all_def, fri_crop)
  all_def <- terra::mask(all_def, shrub_mask)
  
  if (make_plots) {
    shrub_all_def <- all_def
  }
  #---------------------------------------------------------------------------#
  
  
  #----- Extract table stats -------------------------------------------------#
  
  # Overall area
  vv <- values(shrub_mask)
  shrub_cells <- sum(vv, na.rm=T)
  
  rm(vv)
  invisible(gc())
  
  #----- Area in surplus and deficit -----------------------------------------#
  surplus <- terra::ifel(all_def > 1.2, 1, 0)
  vv <- values(surplus)
  surplus_cells <- sum(vv, na.rm=T)
  rm(surplus, vv)
  
  deficit <- terra::ifel(all_def < -1.2, 1, 0)
  vv <- values(deficit)
  deficit_cells <- sum(vv, na.rm=T)
  rm(deficit, vv)
  invisible(gc())
  
  #----- Just right ----------------------------------------------------------#
  goldilocks <- terra::classify(all_def, mat)
  goldilocks_cells <- sum(values(goldilocks), na.rm=T)
  rm(goldilocks)
  invisible(gc())
  
  #----- Area for which we have no data --------------------------------------#
  no_dat <- !is.na(shrub_mask) & is.na(all_def)
  no_data_cells <- sum(values(no_dat))
  rm(no_dat)
  invisible(gc())
  
  #----- Mean and sd of deficit ----------------------------------------------#
  vv <- values(all_def)
  mean_deficit <- mean(vv, na.rm=T)
  sd_deficit   <- sd  (vv, na.rm=T)
  
  summary_dat <- rbind(summary_dat, data.frame(
    Vegetation         = "Shrubland",
    Numcells           = shrub_cells,
    Numsurpluscells    = surplus_cells,
    Numdeficitcells    = deficit_cells,
    Numgoldilockscells = goldilocks_cells,
    Nodatacells        = no_data_cells,
    mean_deficit       = mean_deficit,
    sd_deficit         = sd_deficit
  ))
  
  #---------------------------------------------------------------------------#
  # Shrub WUI 
  #---------------------------------------------------------------------------#
  summary_dat <- rbind(summary_dat, 
                       data.frame(Vegetation = c("Shrub WUI", "Shrub non-WUI"),
                                  Numcells = 0, Numsurpluscells = 0, 
                                  Numdeficitcells = 0, Numgoldilockscells = NA,
                                  Nodatacells = NA,
                                  mean_deficit = 0, sd_deficit = 0))
  
  #----- Overall Shrub WUI area ---------------------------------------------#
  datum <- mask(shrub_mask, wui_crop)
  numcells <- sum(!is.na(values(datum)), na.rm=T)
  summary_dat$Numcells[summary_dat$Vegetation == "Shrub WUI"] <- numcells
  datum <- mask(shrub_mask, wui_crop, inverse=T)
  numcells <- sum(!is.na(values(datum)), na.rm=T)
  summary_dat$Numcells[summary_dat$Vegetation == "Shrub non-WUI"] <- numcells
  rm(datum)
  invisible(gc())
  
  #----- Area in surplus and deficit -----------------------------------------#
  surplus <- terra::ifel(all_def > 1.2, 1, NA)
  surp <- terra::mask(surplus, wui_crop)
  numcells <- sum(values(surp), na.rm=T)
  summary_dat$Numsurpluscells[summary_dat$Vegetation == "Shrub WUI"] <- numcells
  surp <- terra::mask(surplus, wui_crop, inverse = T)
  numcells <- sum(values(surp), na.rm=T)
  summary_dat$Numsurpluscells[summary_dat$Vegetation == "Shrub non-WUI"] <- numcells

  rm(surplus, surp)
  invisible(gc())
  
  deficit <- terra::ifel(all_def < -1.2, 1, NA)
  defr <- terra::mask(deficit, wui_crop)
  numcells <- sum(values(defr), na.rm=T)
  summary_dat$Numdeficitcells[summary_dat$Vegetation == "Shrub WUI"] <- numcells
  defr <- terra::mask(deficit, wui_crop, inverse = T)
  numcells <- sum(values(defr), na.rm=T)
  summary_dat$Numdeficitcells[summary_dat$Vegetation == "Shrub non-WUI"] <- numcells
  rm(deficit, defr)
  invisible(gc())
  
  goldilocks <- terra::classify(all_def, mat)
  goldf <- terra::mask(goldilocks, wui_crop)
  summary_dat$Numgoldilockscells[summary_dat$Vegetation == "Shrub WUI"] <- sum(values(goldf), na.rm=T)
  goldf <- terra::mask(goldilocks, wui_crop, inverse = T)
  summary_dat$Numgoldilockscells[summary_dat$Vegetation == "Shrub non-WUI"] <- sum(values(goldf), na.rm=T)
  rm(goldilocks, goldf)
  invisible(gc())
  
  #----- Mean and sd of deficit ----------------------------------------------#
  defr <- terra::mask(all_def, wui_crop)
  summary_dat$mean_deficit[summary_dat$Vegetation == "Shrub WUI"] <- mean(values(defr), na.rm=T)
  summary_dat$sd_deficit  [summary_dat$Vegetation == "Shrub WUI"] <- sd(values(defr), na.rm=T)
  
  defr <- terra::mask(all_def, wui_crop, inverse=T)
  summary_dat$mean_deficit[summary_dat$Vegetation == "Shrub non-WUI"] <- mean(values(defr), na.rm=T)
  summary_dat$sd_deficit  [summary_dat$Vegetation == "Shrub non-WUI"] <- sd(values(defr), na.rm=T)
  
  rm(shrub_mask, vv, all_def, defr)
  invisible(gc())
  #---------------------------------------------------------------------------#
  
  #----- Clean up and get ready to display -----------------------------------#
  summary_dat$Total_area <- 
    ((summary_dat$Numcells * prod(res(fri_crop)))/10000) * ha_to_acres
  
  summary_dat$Surplus_area <- 
    ((summary_dat$Numsurpluscells * prod(res(fri_crop)))/10000) * ha_to_acres
  
  summary_dat$Deficit_area <- 
    ((summary_dat$Numdeficitcells * prod(res(fri_crop)))/10000) * ha_to_acres
  
  summary_dat$Goldilocks_area <- 
    ((summary_dat$Numgoldilockscells * prod(res(fri_crop)))/10000) * ha_to_acres
  
  summary_dat$Nodata_area <- 
    ((summary_dat$Nodatacells * prod(res(fri_crop)))/10000) * ha_to_acres
  
  
  summary_dat$Region <- region_name
  rm(fri_crop, veg_crop, wui_crop)
  invisible(gc())
  
  
  
  #---------------------------------------------------------------------------#
  # Plotting an overall deficit map
  #---------------------------------------------------------------------------#
  if (make_plots) {
    all_def <- sum(forest_all_def, grass_all_def, shrub_all_def, na.rm=T)
    
    three_cat_mat <- matrix(c(-100000000,         -3, 1,
                              -3        ,       -1.2, 2,
                              -1.2      ,        1.2, 3,
                               1.2      ,          3, 4,
                               3        , 1000000000, 5), byrow=T, nrow=5)
    
    five_cat_mat <- matrix(c(-100000000,         -3, 1,
                             -3        ,         -2, 2,
                             -2        ,       -1.2, 3,
                             -1.2      ,        1.2, 4,
                              1.2      ,          2, 5,
                              2        ,          3, 6,
                              3        , 1000000000, 7), byrow=T, nrow=7)
    all_def3 <- terra::classify(all_def, three_cat_mat)
    levels(all_def3) <- data.frame(ID = 1:5,
                                   category = c("Substantial deficit",
                                                "Moderate deficit",
                                                "Balanced",
                                                "Moderate surplus",
                                                "Substantial surplus"))
    
    all_def5 <- terra::classify(all_def, five_cat_mat)
    levels(all_def5) <- data.frame(ID = 1:7,
                                   category = c("Substantial deficit",
                                                "Moderate deficit",
                                                "Mild deficit",
                                                "Balanced",
                                                "Mild surplus",
                                                "Moderate surplus",
                                                "Substantial surplus"))
    
    region_name <- gsub("/", "-", region_name, fixed=T)
    
    fire_def_plot <- ggplot() + 
          geom_spatraster(data = all_def3) +
          geom_spatvector(data = region_shape, fill=NA) +
          theme_bw() +
          blank_axes +
          ggtitle("Forests fire deficit/surplus") +
          cat3_scale
    
    pp <- veg_map_plot + fire_def_plot + 
      plot_annotation(title = region_name)
    ggsave(paste0("figs/", region_name, "3cat.png"), plot= pp, 
           width = 8, units="in")
    
    #fire_def_plot <- ggplot() + 
    #  geom_spatraster(data = all_def5) +
    #  geom_spatvector(data = region_shape, fill=NA) +
    #  theme_bw() +
    #  blank_axes +
    #  ggtitle("Forests fire deficit/surplus") +
    #  cat5_scale
    
    
    #ggsave(paste0("figs/", region_name, "5cat.png"), plot= print(veg_map_plot + fire_def_plot), 
    #       width = 8, units="in")
    
    rm(all_def, forest_all_def, grass_all_def, shrub_all_def, all_def3, all_def5)
    invisible(gc())
  }
  

  return(summary_dat)
}