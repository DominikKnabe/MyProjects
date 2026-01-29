# Code to create a classification based on positive instances 
# Parts of the code were adapted from https://zenodo.org/records/15189349
rm(list = ls(all.names = TRUE))
gc()

library(terra)
library(raster)
library(dplyr)
library(ggplot2)
library(sf)
library(here)
library(stars)
library(magick)
library(sf)
library(tibble)
library(patchwork)
library(purrr)
library("biscale")
library(grid)
library(scales)

colors <- c("#007800", "#FFDC00", "#FF0000")
custom_palette <- colorRampPalette(c('#006837', '#1a9850', '#66bd63', '#a6d96a', '#d9ef8b', '#fee08b', '#fdae61', '#f46d43', '#d73027', '#a50026'))

# Load objects
    Suitability_graz <- st_read("D:/Masterarbeit/Jupyter/Data/Probabilities/sf_predicted_graz.gpkg",  layer = "layer_name")
    Meteo_pred <- read.csv("D:/Masterarbeit/Jupyter/Data/Meteo/Predicted/Meteo_pred.csv")
    Mos_graz <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Mosquito/Mosquitos_val1_31287.shp")
    
#Cleaning and date manipulation     
    Mos_graz <- Mos_graz %>% rename(date = creation_d)
    Mos_graz$date <- as.character(Mos_graz$date)
    Meteo_pred <- Meteo_pred %>% dplyr::select(date, count, Year, prediction, prediction_exc)

### Plot frequency of index values and presence observations###
    #Breaks for histograms 
    breaks <- seq(0, 1, by = 0.025)
    
    # Calculate presence observations per bin
          presence_counts_graz <- Suitability_graz %>%
            # Cut into bins based on prob_rs
            mutate(bin = cut(prob_pps, breaks = breaks, include.lowest = TRUE)) %>%
            # Group by bin
            group_by(bin) %>%
            # Count the number of presences
            summarise(presence_count = sum(Presence == "1", na.rm = TRUE))
          
            #Add midpoint value for plotting
            # For graz, last bin is missing
            presence_counts_graz <- presence_counts_graz %>%
              mutate(midpoint = 0.0125 + (row_number() - 1) * 0.025)
            
    # Calculate scale factor between the axes
            max_freq_graz <- max(hist(Suitability_graz$prob_rs, breaks = breaks, plot = FALSE)$counts)
            max_presence_graz <- max(presence_counts_graz$presence_count)
           
            #Determine overall maximum for y-limit
            max_presence_overall <- max(max_presence_graz)
            scale_factor_graz <- max_freq_graz / max_presence_graz 
    
      ### Plotting ###
      ## Graz
      a <- ggplot(Suitability_graz, aes(x = prob_rs)) +
        # Histogram for frequencies
        geom_histogram(breaks = breaks, fill = "steelblue", color = "white", alpha = 0.7) +
        # Line for presence counts scaled to primary axis
        geom_line(
          data = presence_counts_graz,
          aes(
            x = midpoint,
            y = presence_count * scale_factor_graz
          ),
          color = "red",
          size = 1
        ) +
        # Points for presence counts scaled
        geom_point(
          data = presence_counts_graz,
          aes(
            x = midpoint,
            y = presence_count * scale_factor_graz
          ),
          color = "red",
          size = 2
        ) +
        # Add secondary axis
        scale_y_continuous(
          name = "Frequency",
          sec.axis = sec_axis(~./scale_factor_graz, name = "Presences")
        ) +
        labs(
          x = "Habitat suitability [%]",
          title = "Frequency of habitat suitability indices and number of observed presences within in each bin"
        ) +
        theme_minimal() +
        theme(
          axis.title.y.right = element_text(angle = 90) 
        )
       plot(a)
       ggsave("D:/Masterarbeit/Figures/Distribution_HS_Index_with_presences.png",
                plot = a,
                height = 4, width = 6, units = "in", dpi = 300)
    
##### Calculate classification thresholds based on percentiles #####
    #### For Habitat suitability #####
    #Identify presence pixels
    obs_graz <- Suitability_graz[Suitability_graz$Presence == "1", ]
    
    # New df with probabilites of presence cells. sort descending and calculate the percentage of observations
      ind_val <- data.frame(Probability = obs_graz$prob_rs) %>%
        arrange(desc(Probability)) %>%  #sort descending
        dplyr::mutate(
          observations = (row_number() / n()) * 100 #percentage of observations
          )
    
      #define threshholds
      p50 <- ind_val %>% filter(abs(observations - 50) == min(abs(observations - 50))) %>% slice(1) %>% pull(Probability)
      p95 <- ind_val %>% filter(abs(observations - 95) == min(abs(observations - 95))) %>% slice(1) %>% pull(Probability)
      p0 <- ind_val %>%  filter(abs(observations) == min(abs(observations))) %>%  slice(1) %>%  pull(Probability)   

      ### For Meteo ###
      #Identify presence pixels
      obs_meteo <- Meteo_pred %>%
        semi_join(Mos_graz, by = "date")
     
      #New df with probabilites of presence cells. sort descending and calculate the percentage of observations
      ind_val_meteo <- data.frame(Predicted = obs_meteo$prediction_exc) %>%
        arrange(desc(Predicted)) %>%
        dplyr::mutate(
          observations = (row_number() / n()) * 100
        )
      
      #define threshholds
      p50_meteo <- ind_val_meteo %>% filter(abs(observations - 50) == min(abs(observations - 50))) %>% slice(1) %>% pull(Predicted)
      p95_meteo <- ind_val_meteo %>% filter(abs(observations - 95) == min(abs(observations - 95))) %>% slice(1) %>% pull(Predicted)
      p0_meteo <- ind_val_meteo %>%  filter(abs(observations) == min(abs(observations))) %>%  slice(1) %>%  pull(Predicted)
  
### Plotting cumulative distribution plots ###
###### For Habiat suitability ###
   theme_set(
     theme_test() +
      theme(
        legend.position.inside = c(0.8, 0.8),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), size = 8),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 8),
        axis.text = element_text(size = 8)
        )
    )

    prob_reclass <- ggplot(ind_val, aes(x = observations, y = Probability)) +
      annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = p95, alpha = 0.6, fill = "#007800") +
      annotate("rect", xmin = 0, xmax = 100, ymin = p95, ymax = p50, alpha = 0.6, fill = "#FFDC00") +
      annotate("rect", xmin = 0, xmax = 100, ymin = p50, ymax = p0, alpha = 0.6, fill = "#FF0000") +
      geom_line(linewidth = 0.6) +
      geom_vline(xintercept = c(50, 95), linetype = "dotted", alpha = 0.6, color = "white", linewidth = 0.6) +
      geom_hline(yintercept = c(p95, p50), color = "white", alpha = 0.6, linewidth = 0.6) +
      xlab("Proportion of presence cells [%]") + ylab("Habitat suitability") +
      theme(
        plot.title = element_text(face = "bold",  hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 14),  
        axis.text.y = element_text(size = 14))
    prob_reclass
    
    # Save plot
    ggsave("D:/Masterarbeit/Figures/Classification_thresholds.png",
           plot = prob_reclass)
  
    #Define classifiaction matrix for bivariate mapping
    reclass_matrix <- matrix(c(0, p95, 1, p95, p50, 2, p50, Inf, 3), ncol = 3, byrow = TRUE)
 
  ### Now the same for temporal model ###
    pred_reclass <- ggplot(ind_val_meteo, aes(x = observations, y = Predicted)) +
      annotate("rect", xmin = 0, xmax = 100, ymin = 0, ymax = p95_meteo, alpha = 0.6, fill = "#007800") +
      annotate("rect", xmin = 0, xmax = 100, ymin = p95_meteo, ymax = p50_meteo, alpha = 0.6, fill = "#FFDC00") +
      annotate("rect", xmin = 0, xmax = 100, ymin = p50_meteo, ymax = p0_meteo, alpha = 0.6, fill = "#FF0000") +
      geom_line(linewidth = 0.6) +
      geom_vline(xintercept = c(50, 95), linetype = "dotted", alpha = 0.6, color = "white", linewidth = 0.6) +
      geom_hline(yintercept = c(p95_meteo, p50_meteo), color = "white", alpha = 0.6, linewidth = 0.6) +
      xlab("Proportion of presence days [%]") + ylab("Predicted count") +
      theme(
        plot.title = element_text(face = "bold",  hjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
    axis.text.x = element_text(size = 12),  
    axis.text.y = element_text(size = 12))
    pred_reclass
  
  # Save plot
  ggsave("D:/Masterarbeit/Figures/Classification_thresholds_meteo.png",
         plot = pred_reclass)
  
  ### For all days including days without observations (e.g., winter) ###
  ind_val_meteo_2 <- data.frame(Predicted = Meteo_pred$prediction_exc) %>%
    arrange(Predicted) %>%
    dplyr::mutate(
      observations = (row_number() / n()) * 100
    )
  
  # Calculate proportion of all days within each class
  prob_low_risk <- mean(Meteo_pred$prediction_exc < p95_meteo, na.rm = TRUE) * 100
  prob_medium_risk <- mean(Meteo_pred$prediction_exc >= p95_meteo & Meteo_pred$prediction_exc < p50_meteo, na.rm = TRUE) * 100 + prob_low_risk
  prob_high_risk <- mean(Meteo_pred$prediction_exc >= p50_meteo & Meteo_pred$prediction_exc <= p0_meteo, na.rm = TRUE) * 100 + prob_medium_risk  
  
  #Dtermine y-lim
  y_max <- max(ind_val_meteo_2$Predicted, na.rm = TRUE)
  
  #Plotting
  pred_reclass_2 <- ggplot(ind_val_meteo_2, aes(x = observations, y = Predicted)) +
    annotate("rect", xmin = 0, xmax = prob_low_risk, ymin = 0, ymax = y_max, alpha = 0.6, fill = "#007800") +
    annotate("rect", xmin = prob_low_risk, xmax = prob_medium_risk, ymin = 0, ymax = y_max, alpha = 0.6, fill = "#FFDC00") +
    annotate("rect", xmin = prob_medium_risk, xmax = 100, ymin = 0, ymax = y_max, alpha = 0.6, fill = "#FF0000") +
    geom_line(linewidth = 0.6) +
    geom_vline(xintercept = c(prob_low_risk, prob_medium_risk), linetype = "dotted", alpha = 0.6, color = "white", linewidth = 0.6) +
    xlab("Proportion of days [%]") + ylab("Predcited count") +
    ggtitle("Cumulative distribution of days between 2022 and 2024") +
    theme(
      plot.title = element_text(face = "bold",  hjust = 0.5),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12))
  
  pred_reclass_2 
  
  # Save plot
  ggsave("D:/Masterarbeit/Figures/Distribution_risk_all_days.png",
         plot = pred_reclass_2,
         height = 4, width = 6, units = "in", dpi = 300)
  
# Apply classification threshholds to Meteo_pred and Suitability
  Meteo_pred <- Meteo_pred %>% mutate(
    classification = case_when(
      prediction_exc < p95_meteo ~0, # low risk
      prediction_exc >= p95_meteo & prediction_exc < p50_meteo ~ 1, # Medium risk
      prediction_exc >= p50_meteo ~ 2  # High risk
    )
  )
  
  Suitability_graz <- Suitability_graz %>% mutate(
    classification = case_when(
      prob_rs < p95  ~ 1, # Low risk
      prob_rs >= p95 & prob_rs < p50 ~ 2, # Medium risk
      prob_rs >= p50  ~ 3  # High risk
    )
  )
  
  ### Plot seasonality of weather risk ###
  Meteo_pred$date <- as.Date(Meteo_pred$date)
  temp_pred_reclass <- ggplot(Meteo_pred, aes(x = date, y = prediction_exc, fill = factor(classification))) +
    geom_col(color = NA) +  
    labs(
      x = "Date",
      y = "Predicted count"
    ) +
    theme(
      legend.position = "right",
      panel.grid.major = element_line(color = "gray80"),
      panel.grid.minor = element_line(color = "gray90"),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      legend.text = element_text(size = 10),
      legend.title = element_text(face = "bold")
    ) +
    scale_fill_manual(
      name = "Hazard level",
      values = c("#007800", "#FFDC00", "#FF0000"),
      breaks = c("0", "1", "2"),
      labels = c("low", "medium", "high")
    ) +
    scale_x_date(
      breaks = date_breaks(width = "6 months"),
      labels = date_format("%b %Y")
    )
  
  ggsave("D:/Masterarbeit/Figures/Temporal_pred_classified.png",
         plot = temp_pred_reclass)
  
###### Apply classification to suitability map #####
    map <- rast("D:/Masterarbeit/Figures/Probabilities_pred_graz_rs.tif")
    
    #reclassify raster using classification matrix
    classified <- classify(map, reclass_matrix)
    classified[classified == 0] <- NA
    plot(classified)
    # Save classified raster
    writeRaster(classified, "D:/Masterarbeit/Figures/Habitat_classified_graz.tif", overwrite=TRUE)
    
    # Now as PNG plot
    df_plot <- as.data.frame(classified, xy = TRUE)
    df_plot$Class <- factor(df_plot$Probabilities_pred, labels = c("Low", "Medium", "High"))
    
    ggplot(df_plot, aes(x = x, y = y, fill = Class)) +
      geom_raster() +
      scale_fill_manual(values = colors, na.translate = FALSE) +
      coord_equal() +
      theme_minimal() +
      theme(legend.position = "right",
            plot.title = element_text(face = "bold", hjust = 0.5)) +
      labs(title = paste("Habitat suitability of Aedes albopictus in Graz"),
           fill = "Suitability_graz",
           x = "Longitude",
           y = "Latitude")
    
    ggsave("D:/Masterarbeit/Figures/Habitat_classified_gg.png",
           plot = last_plot(), width = 10, height = 8, dpi = 300)

    #### Bivariate mapping: Suitability vs Uncertainty ####
    
    #Check distribution of uncertainty
    png("D:/Masterarbeit/Figures/Prediction_Uncertainty_Histogram.png", width = 800, height = 600, res = 150)
    
    hist(Suitability_graz$sd_pred, breaks = 20, xlab = expression(sigma))
    dev.off()
    
      # Create three classes of uncertainties using quantiles
      Suitability_graz <- Suitability_graz %>% mutate(
        uncertainty_lvl = cut(sd_pred, breaks = classInt::classIntervals(
          var = sd_pred, n = 3, style = "quantile"
        )$brks, include.lowest = TRUE, dig.lab = 3, labels = c("1", "2", "3"))
      )

    #Create factor for combined classification
    Suitability_graz <- Suitability_graz %>% mutate(
      facet_cat = as.factor(paste0(classification, " - ", uncertainty_lvl))
    )
    
    color_palette <- c(
      "1 - 1" = "#e8e8e8",   
      "1 - 2" = "#ace4e4",
      "1 - 3" = "#5ac8c8",
      "2 - 1" = "#dfb0d6",
      "2 - 2" = "#a5add3",
      "2 - 3" = "#5698b9",
      "3 - 1" = "#be64ac",
      "3 - 2" = "#8c62aa",
      "3 - 3" = "#3b4994"     
    )
    #Plotting 
    ggplot(Suitability_graz, aes(x = X, y = Y, color = facet_cat)) +
      geom_point() +
      scale_color_manual(values = color_palette) +
      coord_equal() +
      theme_minimal() +
      theme(
        strip.text = element_text(face = "bold"),
        legend.position = "none"
      ) +
      labs(
        x = "Longitude",
        y = "Latitude"
      )
    
    ggsave("D:/Masterarbeit/Figures/Facet_plot.png",
           plot = last_plot(), dpi = 600)
    
    
#########################################
#### Combine spatial and temporal #######
#########################################
    
    #Extract habitat suitability values and transform to sf
    suitability_prob <- Suitability_graz[, c("X", "Y", "prob_rs")]
    suitability_prob <- suitability_prob %>%
      st_as_sf(coords = c("X", "Y"), crs = 31287)
    
    #Normalize suitability values
    suitability_prob$prob_norm <- suitability_prob$prob_rs/sum(suitability_prob$prob_rs)
    
    #Select 2023 for time series
    Meteo_pred_GIF <- Meteo_pred[Meteo_pred$creation_y==2023,]
    
    #Sort according to date
    Meteo_pred_GIF <- Meteo_pred_GIF %>%
      arrange(date)
    
    rownames(Meteo_pred_GIF) <- NULL
    
    # Select only mosquito season
    Meteo_pred_GIF <- Meteo_pred_GIF[121:334,]
    Meteo_pred_GIF <- na.omit(Meteo_pred_GIF)
    
    #Transform normalized df to raster
    target_resolution <- 100
    extent_data <- st_bbox(suitability_prob)  # extent of Raster
    
    coords <- st_coordinates(suitability_prob)
    probs <- suitability_prob$prob_norm
    
    # Create empty raster
    raster <- raster(extent(extent_data), 
                       res = c(target_resolution, target_resolution),
                       crs = st_crs(suitability_prob)$proj4string)
    
    # Transform df to spatial point data
    spdf <- SpatialPointsDataFrame(coords, data = data.frame(probs), proj4string = CRS(st_crs(suitability_prob)$proj4string))
    
    # Rasterize the spdf using the probability field
    raster_norm <- rasterize(spdf, raster, field = "probs", na.rm = TRUE)
    #save raster
    png_filename <- "D:/Masterarbeit/Figures/suitability_raster_normalized.png"
    png(png_filename, width = 800, height = 600)
    plot(raster_norm, main = "suitability_raster_normalized", col = custom_palette(100), zlim = c(0, max(suitability_prob$prob_norm)))
    dev.off()
    
    # Create one abundance raster per day
      #Initilize stack
      abundance_stack <- stack()
      
      for (i in seq_along(Meteo_pred_GIF$prediction_exc)) {
        #distribute daily predicted mosquito counts using the normalized spatial habitat suitability scores
        abundance_raster <- raster_norm * Meteo_pred_GIF$prediction_exc[i] 
        abundance_stack <- addLayer(abundance_stack, abundance_raster)
        }
      
      # Calculate highest possible Value for upper zlim
      max_value <- max(values(abundance_stack), na.rm = TRUE)
      #transform values to consisten scale (0 to 1)
      abundance_stack_norm <- calc(abundance_stack, fun = function(x) x / max_value)
      #New max value for zlim
      max_value_2 <- max(values(abundance_stack_norm), na.rm = TRUE)
      
    #Create GIF
    #Temporarily write images from raster stack
    images <- list()
    for (i in 1:nlayers(abundance_stack_norm)) {
      plot_file <- paste0("D:/Masterarbeit/Figures/Temp/plot_day_", i, ".png")
      png(plot_file, width = 800, height = 600)
      plot(abundance_stack_norm[[i]], main = as.character(Meteo_pred_GIF$date[i]), 
           col = custom_palette(100), 
           zlim = c(0, max_value_2),
           legend.args = list(text = "Hazard", side = 3, line = 1.2, cex = 0.95, adj = 0.5, font = 2))
      dev.off()
      images[[i]] <- plot_file
    }
    
    #Read images
    img_list <- lapply(images, image_read)
    
    #Combine to GIF
    gif <- image_animate(image_join(img_list), fps = 25)
    image_write(gif, "D:/Masterarbeit/Figures/space-time_animation.gif")
    
    print(gif)

    


