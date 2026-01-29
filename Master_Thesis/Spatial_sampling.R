rm(list = ls(all.names = TRUE))
gc()
library(sf) #simple features
library(dplyr) #For mutating dfs
library(tmap) #Plotting geometries
library(mgcv) #For GAM modeling 
library(pROC) #For AUROC validation
library(ggplot2)  #Plotting histograms
library(gridExtra) #For subbplotting
library(corrplot) #Plotting Correlation Matrix
library(Hmisc) #Calculating correlation matrix
library(prioritizr) #Adjancency matrix
library(lubridate) #Date transformation
library(tidyr) #For unnesting lists containing columns

#Load data frame clipped to area of distribution for sampling of pseudo absences
    sf_train_graz <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_train_graz.gpkg",  layer = "layer_name")
    sf_train_graz$LU_Majority_Class <- as.factor(sf_train_graz$LU_Majority_Class)
    sf_train_graz$Presence <- as.factor(sf_train_graz$Presence)

    
#Load Mosquito Sightings and create day of the year attribute
    Mos <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Mosquito/Mosquitos_val1_31287.shp")
    Mos$doy <- yday(Mos$creation_d)

## define color ramp # HEX codes from https://colorbrewer2.org/
    custom_palette <- colorRampPalette(c("#f0f0f0", "#e0ecf4", "#9ebcda", "#8856a7", "#3f007d"))
    custom_palette2 <- colorRampPalette(c("black", "red"))

##########################################
###### Sampling of pseudo absences #######
##########################################

# Seed for reproducibility of random sampling
    set.seed(123)  

#### Remove absences adjacent to presences from sampling process ####
    # Create an adjacency matrix for the points
    am_points <- proximity_matrix(sf_train_graz, distance = 145)  #145 m to include all eight surrounding pixels
    
    # Identify rows of presence points
    presence_indices <- which(sf_train_graz$Presence == 1)
    
    #Create empty df
    remove_indices_0 <- numeric()
    
    for (i in presence_indices) {
      # Get the row indices of absences that are adjacent to this presence point
      adjacent_absences <- which(am_points[i, ] == 1)
      
      # Add those indices to the removal list
      remove_indices_0 <- unique(c(remove_indices_0, adjacent_absences))
    }
    
    # Now, filter out the absence points from sf_0 based on the removal indices
    sf_train_graz_filtered <- sf_train_graz[-remove_indices_0, ]

## Separate presence and absence
    sf_1 <- sf_train_graz[sf_train_graz$Presence == 1, ] 
    sf_0 <- sf_train_graz_filtered[sf_train_graz_filtered$Presence == 0, ]

### PPS and RS sampling of absences
    #Sample size absences
    num_samples <- sum(sf_train_graz$Presence == 1)*1.1 #https://opengeohub.github.io/spatial-prediction-eml/spatiotemporal-machine-learning-for-species-distribution-modeling.html
    
    # PPS-Sampling
    sf_0_pps <- sf_0[sample(1:nrow(sf_0), size = num_samples, 
                            replace = FALSE, 
                            prob = sf_0$Sampling_Effort), ]
    # Random sampling
    sf_0_rs <- sf_0[sample(1:nrow(sf_0), size = num_samples, replace = FALSE), ]
    
    # Combine the two
    sf_resamp_pps = rbind(sf_1, sf_0_pps)
    sf_resamp_rs = rbind(sf_1, sf_0_rs)
    
    sf_resamp_pps <- st_sf(sf_resamp_pps)
    sf_resamp_rs <- st_sf(sf_resamp_rs)
    

### Create map with all Presence/Absence and Sampling Effort
    #Transform Presence to factor
    sf_resamp_pps$Presence <- factor(sf_resamp_pps$Presence, levels = c("0", "1"))  # "0" for Absence, "1" for Presence
    sf_resamp_rs$Presence <- factor(sf_resamp_rs$Presence, levels = c("0", "1"))

    # Create the map with pps-sampled Presence/Absence and Sampling Effort
    map_sf_resamp_pps <- tm_shape(sf_train_graz) +
      tm_squares(
        fill = "Sampling_Effort", 
        fill.scale = tm_scale_continuous(values = custom_palette(10)), 
        size = 1, 
        fill_alpha = 1, 
        col_alpha = 0,
        fill.legend = tm_legend(title = "Sampling Effort", reverse = TRUE)
      ) +
      tm_title("PPS Sampling") + 
      tm_shape(sf_resamp_pps) +
      tm_dots(
        fill = "Presence", 
        fill.scale = tm_scale_categorical(
          values = custom_palette2(2), 
          labels = c("Absence", "Presence") 
        ),  
        size = 0.4, 
        shape = 20, 
        fill_alpha = 1, 
        col_alpha = 1,
        fill.legend = tm_legend(title = "Status")
      ) +
      tm_layout(
        outer.margins = c(0.0, 0.02, 0.02, 0.02),
        inner.margins = c(0, 0, 0, 0)
      )
    
    # Create the map with rs-sampled Presence/Absence and Sampling Effort
    map_sf_resamp_rs <- tm_shape(sf_train_graz) +
      tm_squares(
        fill = "Sampling_Effort",
        fill.scale = tm_scale_continuous(values = custom_palette(10)),
        size = 1,
        fill_alpha = 1,
        col_alpha = 0,
        fill.legend = tm_legend(title = "Sampling Effort", reverse = TRUE)
      ) +
      tm_title("Random Sampling") +
      tm_shape(sf_resamp_rs) +
      tm_dots(
        fill = "Presence",
        fill.scale = tm_scale_categorical(
          values = custom_palette2(2),
          labels = c("Absence", "Presence")
        ),
        size = 0.4,
        shape = 20,
        fill_alpha = 1,
        col_alpha = 1,
        fill.legend = tm_legend(title = "Status")
      )
    
      # Combine both maps
      tmap_arrange(map_sf_resamp_pps, map_sf_resamp_rs, ncol=2)

    
                   
### CD plots - Exploratory data analysis (EDA) ###
    
    # Numerical predictors to analyze  
    vars_eda <- c("LC_Green_Space","LC_Roof_50","LC_Roof_51",
                  "LC_Roof_52", "LC_Sealed_Ground", "LC_Sealed_Ground","LC_Vegetation_High",
                  "LC_Vegetation_Low","LC_Water", "LU_Allot_Grave",
                  "Population_Density","Rain_Inlet","TWI")

    # Drop geometry to work with data as a regular dataframe
    eda_pps <- st_drop_geometry(sf_resamp_pps)
    eda_rs <- st_drop_geometry(sf_resamp_rs)
    
    #Check correlation between predictors
    eda_cor <- eda_pps %>% dplyr::select(-c("Presence", "LU_Majority_Class", "cluster"))
    cormat <- rcorr(as.matrix(eda_cor))
    corrplot(cormat$r)

    # Compute threshold for binary response
    threshold <- sum(eda_pps$Presence == 1) / nrow(eda_pps)


## Create single plots for each variable

# corresponding adjust values, and bins for histograms
adjust_values <- 5
bins_values <- 50

# Initialize an empty list to store plots
plots_list <- list()

# Vector to store AUROC values
auroc_values <- numeric(length(vars_eda))

# Loop through each variable and create a separate plot
for (i in seq_along(vars_eda)) {
  var_name <- vars_eda[i]

  # Fit univariate GAM model (numerical)
  fit <- gam(as.formula(paste("Presence ~ s(", var_name, ", k = 3)")), data = eda_pps, family = binomial)
  probs <- predict(fit, type = "response", newdata = eda_pps)
  
  # Compute AUC (Univariate Analysis)
  auroc <- roc(response = eda_pps$Presence, predictor = probs)
  auroc_value <- round(auc(auroc), 3)
  # Store AUROC value
  auroc_values[i] <- auroc_value
  
  # Create the Conditional Density Plot with Histogram
  p <- ggplot(eda_pps, aes_string(x = var_name)) +
    # Density plot (Conditional Density)
    geom_density(aes(y = ..count../max(..count..), fill = Presence), position = "fill", adjust = adjust_values) +
    scale_fill_manual(values = c("0" = "grey", "1" = "#e34a33")) +
    # Histogram on secondary y-axis with custom bin size
    geom_histogram(aes(y = ..density../max(..density..)), bins = bins_values, fill = "grey", color = "black", alpha = 0.3) +
    # Set primary y-axis between 0 and 1
    scale_y_continuous(limits = c(0, 1), sec.axis = sec_axis(~ . * max(density(eda_pps[[var_name]], na.rm = TRUE)$y), name = "Density (Histogram)")) +
    # Labels and Title
    labs(x = var_name, y = "Conditional Density", title = paste("GAM-based AUROC:", auroc_value)) +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    # Axis limits and threshold line
    xlim(quantile(eda_pps[[var_name]], 0.01, na.rm = TRUE), quantile(eda_pps[[var_name]], 0.99, na.rm = TRUE)) +
    geom_hline(yintercept = threshold, lwd = 2, linetype = "dashed", color = "black")
  # Store the plot in the list
  plots_list[[var_name]] <- p
}

# Print AUROC values
names(auroc_values) <- vars_eda
print(auroc_values) 

grid.arrange(grobs = plots_list, nrow = 3, ncol = 6)

# Histograms of Land Use classes
#Sort in ascending order
eda_pps$LU_Majority_Class <- factor(eda_pps$LU_Majority_Class, 
                                  levels = sort(as.numeric(levels(eda_pps$LU_Majority_Class))))
# Distribution of LU Classes across Presence and Absence
plot1 <- ggplot(eda_pps, aes(x = LU_Majority_Class)) +
  geom_bar(fill = "skyblue") +
  labs(title = "LU Classes across Presences/Absences", x = "Land Use Category", y = "Count") +
  theme_minimal()

# Distribution of LU Classes only across Presence
plot2 <- ggplot(subset(eda_pps, Presence == 1), aes(x = LU_Majority_Class)) +
  geom_bar(fill = "lightcoral") +
  labs(title = "LU Classes across Presences", x = "Land Use Category", y = "Count") +
  theme_minimal()


grid.arrange(plot1, plot2, ncol = 1)

#Distribution of LU Classes across pseudo absences and presences
plot3 <- ggplot(subset(sf_resamp_pps, Presence == 0), aes(x = LU_Majority_Class)) +
  geom_bar(fill = "skyblue") +
  labs(title = "LU Classes across Pseudo absences (M1)", x = "Land Use Category", y = "Count") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 130))

# Distribution of LU Classes only across Presence
plot4 <- ggplot(subset(sf_resamp_rs, Presence == 0), aes(x = LU_Majority_Class)) +
  geom_bar(fill = "skyblue") +
  labs(title = "LU Classes across Pseudo absences (M2)", x = "Land Use Category", y = "Count") +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 130))


grid.arrange(plot2, plot3, plot4,  ncol = 1)

##################
# Count land use across presences and absences of the two approaches
plot_data <- bind_rows(
  eda_pps %>% filter(Presence == 1) %>% count(LU_Majority_Class) %>% mutate(Dataset = "Presence"),
  sf_resamp_pps %>% filter(Presence == 0) %>% count(LU_Majority_Class) %>% mutate(Dataset = "Absence M1"),
  sf_resamp_rs %>% filter(Presence == 0) %>% count(LU_Majority_Class) %>% mutate(Dataset = "Absence M2")
)

# Make sure all classes are present in each dataset
plot_data <- plot_data %>%
  complete(LU_Majority_Class, Dataset, fill = list(n=0))

# Plotting
ggplot(plot_data, aes(x=LU_Majority_Class, y=n, fill=Dataset)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.7) +
  theme_minimal() +
  labs(title="Land use class distribution by presences and pseudo absences",
       x="Land use category",
       y="Count") +
  scale_fill_manual(name = NULL, values=c("Presence"="black", "Absence M1"="skyblue", "Absence M2"="orange")) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(plot = last_plot(), "D:/Masterarbeit/Figures/Distribution_pseudoabsences.png")

#Save resampled data frames
    st_write(sf_resamp_pps, "D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_resamp_pps_no_doy.gpkg",  layer = "layer_name", append = FALSE)
    st_write(sf_resamp_rs, "D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_resamp_rs_no_doy.gpkg",  layer = "layer_name", append = FALSE)
    