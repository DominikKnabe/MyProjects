rm(list = ls(all.names = TRUE))
gc()
library(sperrorest) #Data partitioning for cross-validation
library(magick) # image processing and editing
library(tmap) #thematic mapping
library(dplyr) #For mutating df
library(sf) #Simple features
library(raster) #For raster data
library(ggplot2) #For plotting
library(vip) #visualizing variable importance
library(terra) #Raster analysis
library(tibble) #for data wrangling
library(mgcv) #For model fitting
library(viridisLite) #color palettes
library(pROC) # AUROC calculation
library(mgcViz) # Visualization of partial effects and interaction effects
library(patchwork) #For plotting

## define color ramp # HEX codes from https://colorbrewer2.org/
custom_palette <- colorRampPalette(c('#006837', '#1a9850', '#66bd63', '#a6d96a', '#d9ef8b', '#fee08b', '#fdae61', '#f46d43', '#d73027', '#a50026'))
custom_palette2 <- colorRampPalette(c("black", "red"))

##Read data and manipulations

    # PPS and RS are resampled data from sf_train (cropped to current species distribution) to fit the model
    sf_train_pps <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_resamp_pps_no_doy.gpkg",  layer = "layer_name")
    sf_train_rs <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_resamp_rs_no_doy.gpkg",  layer = "layer_name")
    
    sf_train_graz <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_train_graz.gpkg",  layer = "layer_name")
   
    # Data frames that cover the whole study areas (not cropped to species current distribution)
    sf_pred_graz <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/100m/Postprocessed/sf_pred_graz.gpkg",  layer = "layer_name")
    
    ##Create factors
      sf_train_pps$Presence <- as.factor(sf_train_pps$Presence)
      sf_train_pps$LU_Majority_Class <- as.factor(sf_train_pps$LU_Majority_Class)
      sf_train_pps$cluster <- as.factor(sf_train_pps$cluster)
      
      sf_train_rs$Presence <- as.factor(sf_train_rs$Presence)
      sf_train_rs$LU_Majority_Class <- as.factor(sf_train_rs$LU_Majority_Class)
      
      sf_pred_graz$Presence <- as.factor(sf_pred_graz$Presence)
      sf_pred_graz$LU_Majority_Class <- as.factor(sf_pred_graz$LU_Majority_Class)
      
#########################
### Fitting GAM Model ###
#########################
  
  ## define formulas including interaction terms
  #pps
  fo_all_pps <- Presence ~  
        factor(LU_Majority_Class) + 
        s(Population_Density, k = 3) + 
        s(LU_Allot_Grave, k = 3) +
        s(Rain_Inlet, k = 3) + 
        s(LC_Soil, k = 3) +
        s(TWI, k = 3) + 
        s(LC_Roof_51, k = 3) +
        s(LC_Vegetation_Low, k = 3) +
        s(LC_Vegetation_High, k = 3) +     
        s(LC_Water, k=3) +
        s(LC_Sealed_Ground, k=3) +
        s(LC_Green_Space, k = 3) +
        s(LC_Roof_52, k = 3) +
        s(LC_Roof_50, k = 3) +
        ti(LC_Roof_51, LC_Vegetation_Low, k=c(4,4)) +
        ti(LC_Roof_51, LC_Vegetation_High, k=c(4,4)) +
        ti(LC_Roof_51, LC_Green_Space, k=c(4,4)) + 
        ti(LC_Roof_52, LC_Vegetation_High, k=c(4,4)) +
        ti(LC_Roof_52, LC_Vegetation_Low, k=c(4,4)) + 
        ti(LC_Roof_52, LC_Green_Space, k=c(4,4))
     
      
    #rs
    fo_all_rs <- Presence ~  
      factor(LU_Majority_Class) + 
      s(Population_Density, k = 3) + 
      s(LU_Allot_Grave, k = 3) +
      s(Rain_Inlet, k = 3) + 
      s(LC_Soil, k = 3) +
      s(TWI, k = 3) + 
      s(LC_Roof_51, k = 3) +
      s(LC_Vegetation_Low, k = 3) +
      s(LC_Vegetation_High, k = 3) +     
      s(LC_Water, k=3) +
      s(LC_Sealed_Ground, k=3) +
      s(LC_Green_Space, k = 3) +
      s(LC_Roof_52, k = 3) +
      s(LC_Roof_50, k = 3) +
      s(Sampling_Effort, k=3) +
      ti(LC_Roof_51, LC_Vegetation_Low, k=c(4,4)) +
      ti(LC_Roof_51, LC_Vegetation_High, k=c(4,4)) +
      ti(LC_Roof_51, LC_Green_Space, k=c(4,4)) + 
      ti(LC_Roof_52, LC_Vegetation_High, k=c(4,4)) +
      ti(LC_Roof_52, LC_Vegetation_Low, k=c(4,4)) + 
      ti(LC_Roof_52, LC_Green_Space, k=c(4,4)) 
      
    ##fit models
    #pps
    myfit_all_pps = mgcv::gam(fo_all_pps, data=sf_train_pps, family=binomial, method = "REML", select = TRUE)
    summary(myfit_all_pps)
    
    #rs
    myfit_all_rs = mgcv::gam(fo_all_rs, data=sf_train_rs, family=binomial, method = "REML", select = TRUE)
    summary(myfit_all_rs)
 
    # Viszualize interaction effects  
        #PPS
        vis.gam(myfit_all_pps, plot.type="contour", view=c("LC_Roof_51","LC_Vegetation_Low"), too.far = 0.05)
        vis.gam(myfit_all_pps, plot.type="contour", view=c("LC_Roof_51","LC_Vegetation_High"), too.far = 0.05)
        vis.gam(myfit_all_pps, plot.type="contour", view=c("LC_Roof_52","LC_Green_Space"), too.far = 0.05)
        vis.gam(myfit_all_pps, plot.type="contour", view=c("LC_Roof_52","LC_Vegetation_High"), too.far = 0.05)
        vis.gam(myfit_all_pps, plot.type="contour", view=c("LC_Roof_52","LC_Vegetation_Low"), too.far = 0.05)
        
        #RS
        par(mfrow=c(2,2))
        vis.gam(myfit_all_rs, plot.type="contour", view=c("LC_Roof_51","LC_Green_Space"), too.far = 0.2)
        vis.gam(myfit_all_rs, plot.type="contour", view=c("LC_Roof_51","LC_Vegetation_Low"), too.far = 0.2)
        vis.gam(myfit_all_rs, plot.type="contour", view=c("LC_Roof_51","LC_Vegetation_High"), too.far = 0.2)
        vis.gam(myfit_all_rs, plot.type="contour", view=c("LC_Roof_52","LC_Green_Space"), too.far = 0.2)
        vis.gam(myfit_all_rs, plot.type="contour", view=c("LC_Roof_52","LC_Vegetation_High"), too.far = 0.2)
        vis.gam(myfit_all_rs, plot.type="contour", view=c("LC_Roof_52","LC_Vegetation_Low"), too.far = 0.2)
    
    par(mfrow=c(3,4))
    
    #Create visGAM objects for partial effect plots
    vis_gam_pps <- getViz(myfit_all_pps)
    vis_gam_rs <- getViz(myfit_all_rs)
    
    print(plot(vis_gam_pps, allTerms = T), pages = 1) 
    print(plot(vis_gam_rs, allTerms = T), pages = 1) 
    
    #Partial effects RS
    P_LU_Allot <- plot( sm(vis_gam_rs, 2) )
    P_LU_Allot <- P_LU_Allot + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
      l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
    
    P_TWI <- plot( sm(vis_gam_rs, 5) )
    P_TWI <- P_TWI + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
      l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
    
    P_LC_Roof_51 <- plot( sm(vis_gam_rs, 6) )
    P_LC_Roof_51 <- P_LC_Roof_51 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
      l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
    
    P_LC_Roof_52 <- plot( sm(vis_gam_rs, 12) )
    P_LC_Roof_52 <- P_LC_Roof_52 + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
      l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()

    P_Sampling_Effort <- plot( sm(vis_gam_rs, 14) )
    P_Sampling_Effort <- P_Sampling_Effort + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
      l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
    
    P_LC_Vegetation_High <- plot( sm(vis_gam_rs, 8) )
    P_LC_Vegetation_High <- P_LC_Vegetation_High + l_fitLine(colour = "red") + l_rug(mapping = aes(x=x, y=y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) + 
      l_points(shape = 19, size = 1, alpha = 0.1) + theme_classic()
    
    gridPrint(P_LU_Allot, P_LC_Roof_51, P_LC_Roof_52, P_Sampling_Effort, P_LC_Vegetation_High, P_TWI, ncol=3)

###################################
##### Graz AUROC calculations ######
#################################### 
    
    #Predict habitat suitability values for sf_train (complete df for graz cropped to area of species distribution): for model comparison
    sf_train_graz$prob_pps <- predict(myfit_all_pps, type = "response", newdata = sf_train_graz)
    sf_train_graz$prob_rs <- predict(myfit_all_rs, exclude = c("s(Sampling_Effort)"), type = "response", newdata = sf_train_graz)
    
    #Prediction on resampled dfs for fitting AUROC calculation
    sf_train_pps$prob_pps <- predict(myfit_all_pps, type = "response", newdata = sf_train_pps)
    sf_train_rs$prob_rs <- predict(myfit_all_rs, exclude = c("s(Sampling_Effort)"), type = "response", newdata = sf_train_rs)
    
    #Calculate fitting ROC
    roc_graz_pps <- roc(sf_train_graz$Presence, sf_train_graz$prob_pps)
    roc_graz_rs <- roc(sf_train_graz$Presence, sf_train_graz$prob_rs)
    fitting_roc_pps <- roc(sf_train_pps$Presence, sf_train_pps$prob_pps)
    fitting_roc_rs <- roc(sf_train_rs$Presence, sf_train_rs$prob_rs)
    
    #Calculate fitting AUROC
    auroc_graz_pps <- auc(roc_graz_pps)
    auroc_graz_rs <- auc(roc_graz_rs)
    fitting_auroc_pps <- auc(fitting_roc_pps)
    fitting_auroc_rs <- auc(fitting_roc_rs)
    
    print(auroc_graz_pps)
    print(auroc_graz_rs)
    print(fitting_auroc_pps)
    print(fitting_auroc_rs)
   
###################################    
### Prediction of Probabilities ###
##################################
    #Graz
    sf_pred_graz$prob_pps <- predict(myfit_all_pps, type = "response", newdata = sf_pred_graz)
    sf_pred_graz$prob_rs <- predict(myfit_all_rs, type = "response", exclude = c("s(Sampling_Effort)"), newdata = sf_pred_graz)
   
    #Export predicitons
    st_write(sf_pred_graz, "D:/Masterarbeit/Jupyter/Data/Probabilities/sf_predicted_graz.gpkg",  layer = "layer_name", append = FALSE)
    
    #AUROC after prediction
    roc_graz_pred_pps <- roc(sf_pred_graz$Presence, sf_pred_graz$prob_pps)
    roc_graz_pred_rs <- roc(sf_pred_graz$Presence, sf_pred_graz$prob_rs)
    
    auroc_graz_pred_pps <- auc(roc_graz_pred_pps)
    auroc_graz_pred_rs <- auc(roc_graz_pred_rs)
   
    print(auroc_graz_pred_pps)
    print(auroc_graz_pred_rs)
    
######################################
### Plot predictions ################# 
######################################
    
    #Pixel size
    target_resolution <- 100
    extent_data_graz <- st_bbox(sf_pred_graz)  # extent of Raster
    #Color palette
    custom_palette <- viridisLite::turbo(100)
    
    coords_graz <- st_coordinates(sf_pred_graz)
    
    probs_graz_pps <- sf_pred_graz$prob_pps
    probs_graz_rs <- sf_pred_graz$prob_rs
    
    raster_i_graz <- raster(extent(extent_data_graz), 
                       res = c(target_resolution, target_resolution),
                       crs = CRS("+init=epsg:31287"))
    
    # Transform df to spatial point data
    spdf_graz_pps <- SpatialPointsDataFrame(coords_graz, data = data.frame(probs_graz_pps), proj4string = CRS(st_crs(sf_pred_graz)$proj4string))
    spdf_graz_rs <- SpatialPointsDataFrame(coords_graz, data = data.frame(probs_graz_rs), proj4string = CRS(st_crs(sf_pred_graz)$proj4string))
    
    #Rasterization 
    raster_i_graz_pps <- rasterize(spdf_graz_pps, raster_i_graz, field = "probs_graz_pps", na.rm = TRUE)
    raster_i_graz_rs <- rasterize(spdf_graz_rs, raster_i_graz, field = "probs_graz_rs", na.rm = TRUE)
   
    #Saving
    writeRaster(raster_i_graz_pps, filename = "D:/Masterarbeit/Figures/Probabilities_pred_graz_pps.tif", overwrite = TRUE)
    writeRaster(raster_i_graz_rs, filename = "D:/Masterarbeit/Figures/Probabilities_pred_graz_rs.tif", overwrite = TRUE)
    
    #Plotting 
        #PPS
        png("D:/Masterarbeit/Figures/Habitat_suitability_Graz_pps.png", width = 800, height = 600)
        plot(raster_i_graz_pps, col = custom_palette, zlim = c(0, 1),
             main = "Habitat suitability in Graz (PPS)",
             xlab = "Longitude", ylab = "Latitude")
        
        #RS
        png("D:/Masterarbeit/Figures/Habitat_suitability_Graz_rs.png", width = 800, height = 600)
        plot(raster_i_graz_rs, col = custom_palette, zlim = c(0, 1),
             main = "Habitat suitability in Graz (RS)",
             xlab = "Longitude", ylab = "Latitude")
        
        dev.off()
    
########################################
### FEATURE IMPORTANCE based on auroc###
########################################
    
    #For PPS
        ds <- as_tibble(sf_train_pps)
        mynsim = 30 # number of simulations
        
        set.seed(666)
        
        # Select the relevant columns
        cn_s <- colnames(myfit_all_pps$model)
        cn_s[cn_s == "factor(LU_Majority_Class)"] = "LU_Majority_Class"  # Correcting column name
        train_sel_s <- ds %>%
          dplyr::select(any_of(cn_s)) %>%
          dplyr::select(-Presence)
        
        # Define the target variable
        target_s <- ds$Presence
        
        # Calculate variable importance
        result_pps <- vi_permute(
          nsim = mynsim,
          object = myfit_all_pps, 
          train = train_sel_s, 
          target = target_s, 
          event_level = "second",            # Because predicted instance is 1 and not 0
          metric = "roc_auc",                # auroc based...
          pred_wrapper = predict.gam
        ) %>%
          dplyr::arrange(-Importance) %>%
          dplyr::mutate(process = "Mosquito Presence")
        ###
        ##
        #
        
        result_pps
        
        # Ensure the Variable column is ordered by Importance in descending order
        result_pps <- result_pps %>%
          arrange(desc(Importance)) %>%  # Sort data frame descending 
          mutate(Variable = factor(Variable, levels = rev(Variable)))  # Reverse factor levels
        result_pps
        
        #Plot feature importance
        Importance_Mosquito_pps <- ggplot(result_pps, aes(x = Importance, y = Variable)) +
          geom_point(color = "#747474", size = 3) +  
          geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev), 
                         height = 0.2, color = "#747474", linewidth = 0.5) +  # Increase line width
          theme_minimal(base_size = 16) +  # Increase font size
          labs(title = "Feature Importance M1", x = "Importance", y = "Variable") +
          theme(axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                plot.title = element_text(size = 18, face = "bold"),
                panel.border = element_rect(color = "black", fill = NA, size = 1))
        Importance_Mosquito_pps
        
        ggsave(plot=Importance_Mosquito_pps, "D:/Masterarbeit/Figures/Feature_Importance_PPS.png")
        
    #For RS
        ds <- as_tibble(sf_train_rs)
        mynsim = 30 # number of simulations 
        
        set.seed(666)
        
        # Select the relevant columns
        cn_s <- colnames(myfit_all_rs$model)
        cn_s[cn_s == "factor(LU_Majority_Class)"] = "LU_Majority_Class"  # Correcting column name
        train_sel_s <- ds %>%
          dplyr::select(any_of(cn_s)) %>%
          dplyr::select(-Presence)
        
        # Define the target variable
        target_s <- ds$Presence
        
        # Calculate variable importance
        result_rs <- vi_permute(
          nsim = mynsim,
          object = myfit_all_rs, 
          train = train_sel_s, 
          target = target_s, 
          event_level = "second",            # Because my predicted instance is 1 and not 0
          metric = "roc_auc",                # auroc based...
          pred_wrapper = predict.gam
        ) %>%
          dplyr::arrange(-Importance) %>%
          dplyr::mutate(process = "Mosquito Presence")
        ###
        ##
        #
        
        result_rs
        
        # Ensure the Variable column is ordered by Importance in descending order
        result_rs <- result_rs %>%
          arrange(desc(Importance)) %>%  # Sort data frame
          mutate(Variable = factor(Variable, levels = rev(Variable)))  # Reverse factor levels
        result_rs
        
        # Plotting 
        Importance_Mosquito_rs <- ggplot(result_rs, aes(x = Importance, y = Variable)) +
          geom_point(color = "#747474", size = 3) +  
          geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev), 
                         height = 0.2, color = "#747474", linewidth = 0.5) +  # Increase line width
          theme_minimal(base_size = 16) +  # Increase font size
          labs(title = "Feature Importance M2", x = "Importance", y = "Variable") +
          theme(axis.title.x = element_text(size = 15),
                axis.title.y = element_text(size = 15),
                plot.title = element_text(size = 18, face = "bold"),
                panel.border = element_rect(color = "black", fill = NA, size = 1))
        Importance_Mosquito_rs

        ggsave(plot=Importance_Mosquito_rs, "D:/Masterarbeit/Figures/Feature_Importance_RS.png")
        
################################################
### Spatial cross validation using sperrorest (Brenning, 2012)###
################################################
    #Transform geometry to coords
    coords_pps <- st_coordinates(sf_train_pps)
    coords_rs <- st_coordinates(sf_train_rs)
    
    # Bind coordinates back as separate columns
    sf_train_pps <- sf_train_pps %>%
      st_drop_geometry() %>%
      cbind(coords_pps)
    
    sf_train_rs <- sf_train_rs %>%
      st_drop_geometry() %>%
      cbind(coords_rs)

    nfolds = 5 #Five fold since there are 5 clusters
    nreps = 1 #one repetition is enough since partitioning is factor based
    
    #Factor partioning of data into train and test
    parti_pps <- partition_factor_cv(sf_train_pps, coords = c("X", "Y"), nfold = nfolds, repetition = nreps, seed1 = 123, fac = "cluster")
    parti_rs <- partition_factor_cv(sf_train_rs, coords = c("X", "Y"), nfold = nfolds, repetition = nreps, seed1 = 123, fac = "cluster")

    #Initzialize data frame for results
    results_cv_mos_pps <- data.frame(Repetition = integer(), Fold = integer(), AUROC = numeric(), stringsAsFactors = FALSE)
    results_cv_mos_rs <- data.frame(Repetition = integer(), Fold = integer(), AUROC = numeric(), stringsAsFactors = FALSE)
    
    #Loop through folds and reps to fit model, predict and calculate AUROC
    #pps
    for (j in 1:nreps) {
      partiloop <- parti_pps[[j]]  # Access the jth element
      for (i in 1:nfolds) {
        first <- partiloop[[i]][[2]]  # Access the second element of the ith fold
        test <- sf_train_pps[first,]; ntesti <- nrow(test); train <- sf_train_pps[-first,]  # Exclude the train set
        # Fit and calculate AUROC
        myfit <- mgcv::gam(fo_all_pps, data = train, family = binomial, method = "REML", select = TRUE)
        test$prob <- predict.gam(myfit, type = "response", newdata = test)
        myroc <- pROC::roc(response = test$Presence, predictor = test$prob, auc = TRUE)
        auroc <- round(myroc$auc, 5)
        print(auroc)
        # Store results
        results_cv_mos_pps <- rbind(results_cv_mos_pps, data.frame(Repetition = j, Fold = i, AUROC = auroc))}}
    
    #rs
    for (j in 1:nreps) {
      partiloop <- parti_rs[[j]]  # Access the jth element
      for (i in 1:nfolds) {
        first <- partiloop[[i]][[2]]  # Access the second element of the ith fold
        test <- sf_train_rs[first,]; ntesti <- nrow(test); train <- sf_train_rs[-first,]  # Exclude the test set
        # Fit and calculate AUROC
        myfit <- mgcv::gam(fo_all_rs, data = train, family = binomial, method = "REML", select = TRUE)
        test$prob <- predict.gam(myfit, type = "response", exclude = "s(Sampling_Effort)", newdata = test)
        myroc <- pROC::roc(response = test$Presence, predictor = test$prob, auc = TRUE)
        auroc <- round(myroc$auc, 5)
        print(auroc)
        # Store results
        results_cv_mos_rs <- rbind(results_cv_mos_rs, data.frame(Repetition = j, Fold = i, AUROC = auroc))}}

    #average results
    average_results_pps <- aggregate(AUROC ~ Repetition, data = results_cv_mos_pps, FUN = mean)
    average_results_rs <- aggregate(AUROC ~ Repetition, data = results_cv_mos_rs, FUN = mean)
    print(average_results_pps)
    print(average_results_rs)

    ######################################################
    ##### Random CV for uncertainty calculation ##########
    ######################################################
    
    coords_2 <- st_coordinates(sf_pred_graz)
    
    sf_pred_graz <- sf_pred_graz %>%
      st_drop_geometry() %>%
      cbind(coords_2)
    
    nfolds = 5 #folds
    nreps = 10 #repetitions
    
    # Initzialize list to store predicitons of each fold
    predictions_list <- vector("list", length = nfolds)
    
    #Random partioning of data into train and test
    parti <- partition_cv(sf_train_rs, nfold = nfolds, repetition = nreps, seed1 = 123)
    
    #initzialize empty df
    results_cv_mos <- data.frame(Repetition = integer(), Fold = integer(), AUROC = numeric(), stringsAsFactors = FALSE)
    
    #Loop through folds and reps to fit model and predict habitat suitability scores
    for (j in 1:nreps) {
      partiloop <- parti[[j]]  # Access the jth element
      for (i in 1:nfolds) {
        first <- partiloop[[i]][[2]]  # Access the second element of the ith fold
        train <- sf_train_rs[-first,]  # Exclude the test set
        pred <- sf_pred_graz
        # Fitting
        myfit <- mgcv::gam(fo_all_rs, data = train, family = binomial, method = "REML", select = TRUE)
        #Prediciton
        pred$prob <- predict.gam(myfit, type = "response", exclude = c("s(Sampling_Effort)"), newdata = pred)
                #Store predictions of each fold
        predictions_list[[i]] <- rbind(predictions_list[[i]], pred[, c("X", "Y", "prob")])
      }
      #Print repetition number
      print(j)
    }
    #Combine lists
    all_predictions <- do.call(rbind, predictions_list)
    
    # Calculate standard deviation of predictions by location
    uncertainty <- all_predictions %>%
      group_by(X, Y) %>%
      summarise(sd_pred = sd(prob, na.rm = TRUE))
    
    #Plot uncertainties
    ggplot(uncertainty, aes(x = X, y = Y, fill = sd_pred)) +
      geom_tile() +
      scale_fill_viridis_c(option = "viridis", name = "Prediction\nStd Dev") +
      coord_fixed() +
      labs(title = "Spatial Uncertainty in Predictions",
           x = "Longitude",
           y = "Latitude") +
      theme_minimal()
    
    ggsave(plot = last_plot(), "D:/Masterarbeit/Figures/Model_uncertainties.png")
    #Merge uncertainties and sf_pred_graz
    sf_pred_graz <- sf_pred_graz %>%
      left_join(uncertainty %>% dplyr::select(X, Y, sd_pred), by = c("X", "Y"))
    
    st_write(sf_pred_graz, "D:/Masterarbeit/Jupyter/Data/Probabilities/sf_predicted_graz.gpkg",  layer = "layer_name", append = FALSE)
    
   