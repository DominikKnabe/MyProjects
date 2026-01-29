rm(list = ls(all.names = TRUE))
gc()
library(sf) #simple features
library(dplyr) #For mutating dfs
library(terra) #Raster analysis
library(readxl) #Read excel files
library(ggplot2) # For plotting

#Load Convex Hull (area of distribution (from QGIS) as Sampling Area) and boundaries
convex_hull_graz <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Area_Of_Equilibrium/Convex_Hull_31287.shp")
city_bound_graz  <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Stadtgrenze/Stadtgrenze_31287.shp")
districts_graz <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Stadtgrenze/Bezirke.geojson")

#Clean up districts graz
districts_graz <- districts_graz[1:17, ] %>% dplyr::select(name, geometry)

#Order df according to numbering of districts and create attribute district_id
districts_graz <- districts_graz[c(1, 3, 4, 5, 6, 7, 2, 8:17), ]
districts_graz <- districts_graz %>% mutate(
  district_id = row_number()
)
#Aggregate districts to clusters for spatial CV
districts_graz <- districts_graz %>% mutate(
  cluster = case_when(
    district_id %in% c(1, 2, 3, 6) ~ 1,
    district_id %in% c(4, 5) ~ 2,
    district_id %in% c(7, 8, 9) ~ 3,
    district_id %in% c(14, 15, 16, 17) ~ 4,
    district_id %in% c(10, 11, 12, 13) ~ 5,
  )
)

#Reproject layer
districts_graz <- st_transform(districts_graz, crs = 31287)

#Plot clusters
ggplot(data = districts_graz) +
  geom_sf(aes(fill = factor(cluster)), color = "black") +  
  scale_fill_viridis_d(name = "Area") + 
  theme_minimal() +  
  labs(title = "Areas for spatial partioning of data",
       x = "Longitude",
       y = "Latitude") +
  theme(legend.position = "bottom")

ggsave(plot = last_plot(), "D:/Masterarbeit/Figures/Clusters_CV.png")

# Load Raster, that was created in QGIS, and rename bands
    # For Graz
    modeling_raster_graz<-rast("D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Modeling_Raster_100m.tif")
    new_band_names <- c("LC_Agriculture", "LC_Construction", "LC_Green_Space","LC_OpenSoil","LC_RoofSealing",
                        "LC_Roof_50", "LC_Roof_51", "LC_Roof_52","LC_Sealed_Ground","LC_Vegetation_High","LC_Vegetation_Low","LC_Water", "LU_Class50", "LU_Class51", "LU_Class52","LU_Majority_Class", "LU_Majority_Use",
                        "LU_Allot_Grave","Population_Density","Presence","Rain_Inlet","TWI","Sampling_Effort") 
    names(modeling_raster_graz) <- new_band_names
    writeRaster(modeling_raster_graz, 
                filename = "D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/modeling_raster_graz_pp.tif", 
                overwrite = TRUE)
    
    plot(modeling_raster_graz)
    
    
#Raster to dataframe
    df_pred_graz <- as.data.frame(modeling_raster_graz, xy = TRUE, na.rm = TRUE) #resulting point corrdinates represent pixel centroids
    
### Regrouping of Land Use categories ###
    #Graz: define subsets to be grouped 
    LU_10 <- c(14, 3, 10, 18, 19, 5, 29, 20)
    LU_11 <- c(21, 9) 
    LU_20 <- c(28)
    LU_30 <- c(8, 25)
    LU_40 <- c(2, 4, 26)
    LU_50 <- c(12, 23) 
    LU_51 <- c(6, 11)  
    LU_52 <- c(16, 17, 13)  
    LU_60 <- c(1, 7, 15, 22, 24, 27)
    
    #Reclassify land use classes
    df_pred_graz <- df_pred_graz %>%
      mutate(LU_Majority_Class = case_when(
        LU_Majority_Use %in% LU_10 ~ 10,
        LU_Majority_Use %in% LU_11 ~ 11,
        LU_Majority_Use %in% LU_20 ~ 20,
        LU_Majority_Use %in% LU_30 ~ 30,
        LU_Majority_Use %in% LU_40 ~ 40,
        LU_Majority_Use %in% LU_50 ~ 50,
        LU_Majority_Use %in% LU_51 ~ 51,
        LU_Majority_Use %in% LU_52 ~ 52,
        LU_Majority_Use %in% LU_60 ~ 60,
      ))

    #Proportion of each land use class
    prop.table(table(df_pred_graz$LU_Majority_Class))

### Regrouping of Land Cover Attributes ###
df_pred_graz <- df_pred_graz %>%
  mutate(LC_Soil = LC_OpenSoil + LC_Construction + LC_Agriculture) %>%   
    dplyr::select(-LC_OpenSoil, -LC_Construction, -LC_Agriculture)

#Transform presence variable and LU_Majority_Use to factor
df_pred_graz$Presence <- as.factor(df_pred_graz$Presence)
df_pred_graz$LU_Majority_Class <- as.factor(df_pred_graz$LU_Majority_Class)

#Remove unwanted attributes
df_pred_graz$LU_Majority_Use <- NULL
df_pred_graz$LU_Class51 <- NULL
df_pred_graz$LU_Class52 <- NULL
df_pred_graz$LU_Class50 <- NULL
df_pred_graz$LC_RoofSealing <- NULL

#Clip data to city boundaries
sf_pred_graz <- st_as_sf(df_pred_graz, coords = c("x", "y"), crs = 31287)
sf_pred_graz <- st_intersection(sf_pred_graz, city_bound_graz)
sf_pred_graz <- sf_pred_graz %>% dplyr::select(-GEMEINDE, -area)

# Remove water pixels
sf_pred_graz <- sf_pred_graz[sf_pred_graz$LU_Majority_Class != 30, ]

#Add spatial cluster to sf_pred_graz. For spatial cross-validation, therefore only for Graz
sf_pred_graz <- st_join(sf_pred_graz, districts_graz[, c("geometry", "cluster")])

#Divide by 100 so LC and LU variables represent the percentage of pixel being covered
    sf_pred_graz$LC_Green_Space <- sf_pred_graz$LC_Green_Space/100 
    sf_pred_graz$LC_Roof_50 <- sf_pred_graz$LC_Roof_50/100 
    sf_pred_graz$LC_Roof_51 <- sf_pred_graz$LC_Roof_51/100 
    sf_pred_graz$LC_Roof_52 <- sf_pred_graz$LC_Roof_52/100 
    sf_pred_graz$LC_Sealed_Ground <- sf_pred_graz$LC_Sealed_Ground/100 
    sf_pred_graz$LC_Vegetation_High <- sf_pred_graz$LC_Vegetation_High/100 
    sf_pred_graz$LC_Vegetation_Low <- sf_pred_graz$LC_Vegetation_Low/100 
    sf_pred_graz$LC_Water <- sf_pred_graz$LC_Water/100 
    sf_pred_graz$LU_Allot_Grave <- sf_pred_graz$LU_Allot_Grave/100 
    sf_pred_graz$LC_Soil <- sf_pred_graz$LC_Soil/100
   
#Create training set by clipping with area of mosquito spread 
    sf_train_graz <- st_intersection(sf_pred_graz, convex_hull_graz)
    sf_train_graz <- sf_train_graz %>% dplyr::select(-id, -area, -perimeter)

# Plot data frames
    plot(sf_pred_graz)
    plot(sf_train_graz)
   
#Save dfs
    st_write(sf_train_graz, "D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_train_graz.gpkg",  layer = "layer_name", append = FALSE)
    st_write(sf_pred_graz, "D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_pred_graz.gpkg",  layer = "layer_name", append = FALSE)
    
#Statsitics SE
    mean(sf_pred_graz$Sampling_Effort, na.rm = TRUE)
    sd(sf_pred_graz$Sampling_Effort, na.rm = TRUE)
    