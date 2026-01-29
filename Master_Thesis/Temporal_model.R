rm(list = ls(all.names = TRUE))
gc()
library(dplyr) #For mutating dfs
library(mgcv) #For model fitting
library(MASS) #For negative binomial distribution
library(sperrorest) #Partitioning in cross validation
library(vip) #Feature Importance
library(ggplot2)  #Plotting histograms
library(mgcViz) #visualization of partial effects

#load postprocessed data
df_train <- read.csv("D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/Meteo/df_train_meteo_MA.csv")
#Rename columns
df_train <- df_train %>%
  dplyr::rename(Sampling_Effort = SE) %>%
  dplyr::rename(year = creation_y)

df_train$year <- as.factor(df_train$year)

########################
######### GAMM ##########
########################

# 1. Set: Daily values and 28-days time-lagged averages/accumulations
fo <- count ~
  s(t_max, k = 3) +
  s(rh_mean, k = 3) +
  s(ws_mean, k = 3) +
  s(ws_cv, k = 3) +
  s(p, k = 3) +
  s(t_mean_mean_28, k = 3) + 
  s(rh_mean_mean_28, k = 3) +
  s(ws_mean_mean_28, k = 3) + 
  s(p_acc_28, k = 3) +
  s(doy, k = 12, bs = "cc") +  s(year, bs = "re") +  s(Sampling_Effort, k = 3)

# 2. Set: Daily values and 28-days time-lagged variations
fo_cv <- count ~
  s(t_max, k = 3) +
  s(rh_mean, k = 3) +
  s(ws_mean, k = 3) +
  s(ws_cv, k = 3) +
  s(p, k = 3) +
  s(t_mean_cv_28, k = 3) + 
  s(ws_mean_cv_28, k = 3) + 
  s(rh_mean_cv_28, k = 3) +
  s(p_cv_28, k = 3) +
  s(doy, k = 12, bs = "cc") +   s(year, bs = "re") +  s(Sampling_Effort, k = 3)

# # 1. Set: Daily values, 28-days time-lagged averages/accumulations and 28-days time-lagged variations
fo_all <- count ~ 
  s(t_max, k = 3) + 
  s(rh_mean, k = 3) + 
  s(ws_mean, k = 3) + 
  s(ws_cv, k = 3) +
  s(p, k = 3) + 
  s(t_mean_mean_28, k = 3) + 
  s(t_mean_cv_28, k = 3) +
  s(p_acc_28, k = 3) + 
  s(p_cv_28, k = 3) +
  s(rh_mean_mean_28, k = 3) + 
  s(rh_mean_cv_28, k = 3) + 
  s(ws_mean_mean_28, k = 3) + 
  s(ws_mean_cv_28, k = 3) + 
  s(doy, k = 12, bs = "cc") +  s(year, bs = "re") +   s(Sampling_Effort, k = 3)

#### Model Fitting using three sets of predictors ####
myfit = mgcv::gam(fo, data=df_train, family= nb(), select = TRUE)
summary(myfit)

plot(myfit, page=1, scale=-1)

myfit_cv = mgcv::gam(fo_cv, data=df_train, family= nb(), select = TRUE)
summary(myfit_cv)

plot(myfit_cv, page=1, scale=0)

myfit_all = mgcv::gam(fo_all, data=df_train, family= nb(), select = TRUE)
summary(myfit_all)

plot(myfit_all, page=1, scale=0)

#R-squared
summary(myfit)$r.sq        
summary(myfit_cv)$r.sq     
summary(myfit_all)$r.sq 

##########################
#### Partial effects #####
##########################
    
#cretae visGAM object
    vis_gam <- getViz(myfit)
    print(plot(vis_gam, allTerms = T), pages = 1) 
    
    #PLot partial effects
    P_t_mean_28 <- plot(sm(vis_gam, 6))
    P_t_mean_28 <- P_t_mean_28 +
      l_fitLine(colour = "red") +
      l_rug(mapping = aes(x = x, y = y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) +
      l_points(shape = 19, size = 1, alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_classic()
    
    P_doy <- plot(sm(vis_gam, 10))
    P_doy <- P_doy +
      l_fitLine(colour = "red") +
      l_rug(mapping = aes(x = x, y = y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) +
      l_points(shape = 19, size = 1, alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_classic()
    
    P_t_max <- plot(sm(vis_gam, 1))
    P_t_max <- P_t_max +
      l_fitLine(colour = "red") +
      l_rug(mapping = aes(x = x, y = y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) +
      l_points(shape = 19, size = 1, alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_classic()
    
    P_rh_mean <- plot(sm(vis_gam, 2))
    P_rh_mean <- P_rh_mean +
      l_fitLine(colour = "red") +
      l_rug(mapping = aes(x = x, y = y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) +
      l_points(shape = 19, size = 1, alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_classic()
    
    P_SE <- plot(sm(vis_gam, 12))
    P_SE <- P_SE +
      l_fitLine(colour = "red") +
      l_rug(mapping = aes(x = x, y = y), alpha = 0.8) +
      l_ciLine(mul = 5, colour = "blue", linetype = 2) +
      l_points(shape = 19, size = 1, alpha = 0.1) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_classic()
    
    gridPrint(P_t_mean_28, P_doy, P_t_max, P_rh_mean, P_SE, ncol = 3)

##########################
#### Cross validation ####
##########################
    
nfolds = 5 #folds
nreps = 5 #repetitions

# Random partioninig
parti_factor <- partition_cv(df_train, nfold = nfolds, repetition = nreps, seed1 = 123)

# Initzialize empty df for results
results_cv_meteo <- data.frame(
  Repetition = integer(),
  Fold = integer(),
  R2 = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)
results_cv_meteo_cv <- data.frame(
  Repetition = integer(),
  Fold = integer(),
  R2 = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)
results_cv_meteo_all <- data.frame(
  Repetition = integer(),
  Fold = integer(),
  R2 = numeric(),
  RMSE = numeric(),
  stringsAsFactors = FALSE
)

# CV loop for fo
for (j in 1:nreps) {
  partiloop <- parti_factor[[j]]  # Access the j-th element
  for (i in 1:nfolds) {
    first <- partiloop[[i]]$test  # Access the test set of the i-th fold
    test <- df_train[first,]; ntesti <- nrow(test);
    train <- df_train[-first,]  # Exclude the test set
    # Fit and predict
    fit_cv <- mgcv::gam(fo, data = train, family = nb(), select = TRUE)
    test$pred <- predict.gam(fit_cv, type = "response", newdata = test)
    
    #RMSE
    mse <- mean((test$count - test$pred)^2, na.rm = TRUE) 
    rmse <- sqrt(mse)
    
    #Calculate R2
    ss_res <- sum((test$count - test$pred)^2, na.rm = TRUE)
    ss_tot <- sum((test$count - mean(test$count, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- 1 - ss_res/ss_tot
    
    # Store results
    results_cv_meteo <- rbind(results_cv_meteo, data.frame(Repetition = j, Fold = i, R2 = r2, RMSE = rmse))
  }
}

# CV loop  for fo_cv    
for (j in 1:nreps) {
  partiloop <- parti_factor[[j]]  # Access the j-th element
  for (i in 1:nfolds) {
    first <- partiloop[[i]]$test  # Access the test set of the i-th fold
    test <- df_train[first,]; ntesti <- nrow(test);
    train <- df_train[-first,]  # Exclude the test set
    # Fit and predict
    fit_cv <- mgcv::gam(fo_cv, data = train, family = nb(), select = TRUE)
    test$pred <- predict.gam(fit_cv, type = "response", newdata = test)
    
    #RMSE
    mse <- mean((test$count - test$pred)^2, na.rm = TRUE) 
    rmse <- sqrt(mse)
    
    #Calculate R2
    ss_res <- sum((test$count - test$pred)^2, na.rm = TRUE)
    ss_tot <- sum((test$count - mean(test$count, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- 1 - ss_res/ss_tot
    
    # Store results
    results_cv_meteo_cv <- rbind(results_cv_meteo_cv, data.frame(Repetition = j, Fold = i, R2 = r2, RMSE = rmse))
  }
}

# CV loop  for fo_all    
for (j in 1:nreps) {
  partiloop <- parti_factor[[j]]  # Access the j-th element
  for (i in 1:nfolds) {
    first <- partiloop[[i]]$test  # Access the test set of the i-th fold
    test <- df_train[first,]; ntesti <- nrow(test);
    train <- df_train[-first,]  # Exclude the test set
    # Fit and predict
    fit_cv <- mgcv::gam(fo_all, data = train, family = nb(), select = TRUE)
    test$pred <- predict.gam(fit_cv, type = "response", newdata = test)
    
    #RMSE
    mse <- mean((test$count - test$pred)^2, na.rm = TRUE) 
    rmse <- sqrt(mse)
    
    #Calculate R2
    ss_res <- sum((test$count - test$pred)^2, na.rm = TRUE)
    ss_tot <- sum((test$count - mean(test$count, na.rm = TRUE))^2, na.rm = TRUE)
    r2 <- 1 - ss_res/ss_tot
    
    # Store results
    results_cv_meteo_all <- rbind(results_cv_meteo_all, data.frame(Repetition = j, Fold = i, R2 = r2, RMSE = rmse))
  }
}

# Print results
average_results_r2 <- aggregate(R2 ~ Repetition, data = results_cv_meteo, FUN = mean)
average_results_r2_cv <- aggregate(R2 ~ Repetition, data = results_cv_meteo_cv, FUN = mean)
average_results_r2_all <- aggregate(R2 ~ Repetition, data = results_cv_meteo_all, FUN = mean)
print(mean(average_results_r2$R2))
print(mean(average_results_r2_cv$R2))
print(mean(average_results_r2_all$R2))

average_results_rmse <- aggregate(RMSE ~ Repetition, data = results_cv_meteo, FUN = mean) 
average_results_rmse_cv <- aggregate(RMSE ~ Repetition, data = results_cv_meteo_cv, FUN = mean) 
average_results_rmse_all <- aggregate(RMSE ~ Repetition, data = results_cv_meteo_all, FUN = mean) 
print(mean(average_results_rmse$RMSE)) 
print(mean(average_results_rmse_cv$RMSE)) 
print(mean(average_results_rmse_all$RMSE))


############################################
################# Prediction ###############
############################################
#Predictors to be zeroed out during prediction
my_exclude <- c("s(year)" , "s(Sampling_Effort)")  

#predict
df_train$prediction <- predict(myfit, type = "response", newdata = df_train)
df_train$prediction_exc <- predict(myfit, type = "response", exclude = my_exclude, newdata = df_train)

#difference between prediction and true count
print(sum(df_train$count)-sum(df_train$prediction, na.rm = TRUE))

#Plot without my_exclude
plot <- ggplot(df_train, aes(x = as.Date(date))) + 
  geom_line(aes(y = count, color = "Observed Count"), linewidth = 0.7) +  
  geom_line(aes(y = prediction, color = "Predicted Count"), linewidth = 0.9, linetype = "dashed") +
  labs(
    x = "Date",
    y = "Count"
  ) +
  scale_color_manual(values = c("Observed Count" = "blue", "Predicted Count" = "red")) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_blank(),            
    legend.text = element_text(size = 12),    
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )
    print(plot)

    ggsave(plot = last_plot(), "D:/Masterarbeit/Figures/Meteo_pred.png")
#Plot with my_exclude
    plot_exc <- ggplot(df_train, aes(x = as.Date(date))) + 
      geom_line(aes(y = count, color = "Observed Count"), linewidth = 0.7) +  
      geom_line(aes(y = prediction_exc, color = "Predicted Count"), linewidth = 0.9, linetype = "dashed") +
      labs(
        x = "Date",
        y = "Count"
      ) +
      scale_color_manual(values = c("Observed Count" = "blue", "Predicted Count" = "red")) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_blank(),            
        legend.text = element_text(size = 12),    
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    
    print(plot_exc)
    ggsave(plot = last_plot(), "D:/Masterarbeit/Figures/Meteo_pred_exc.png")
#Root Mean Square Error
    mse_gam <- mean((df_train$count - df_train$prediction)^2, na.rm = TRUE)
    rmse_gam <- sqrt(mse_gam)
    rmse_gam

write.csv(df_train, "D:/Masterarbeit/Jupyter/Data/Meteo/Predicted/Meteo_pred.csv")

############################################
########## Feature Importance ##############
############################################
#transfrom to tibble
ds <- as_tibble(df_train)
mynsim = 30  
#Seed for reproducability
set.seed(666)

# Select the relevant columns
cn_s <- colnames(myfit$model)
train_sel_s <- ds %>%
  dplyr::select(any_of(cn_s)) %>%
  dplyr::select(-count) 

# Define the target variable
target_s <- ds$count

#define pred_wrapper since we want to exclude year and Sampling_Effort
pfun_prob <- function(object, newdata) { 
  # prediction wrapper  
  predict.gam(object, newdata = newdata, type = "response", exclude = my_exclude, select = TRUE)}
 
# Calculate variable importance
result_s <- vi_permute(
  nsim = mynsim,
  object = myfit, 
  train = train_sel_s, 
  target = target_s, 
  metric = "RMSE",            
  pred_wrapper = pfun_prob
) %>%
  dplyr::arrange(-Importance) %>%
  dplyr::mutate(process = "Mosquito Presence")
###

# Ensure the Variable column is ordered by Importance in descending order
result_s <- result_s %>%
  arrange(desc(Importance)) %>%  # Sort data frame
  mutate(Variable = factor(Variable, levels = rev(Variable)))  # Reverse factor levels
result_s

# Sort the data by Importance in descending order
Importance_Mosquito <- ggplot(result_s, aes(x = Importance, y = Variable)) +
  geom_point(color = "#747474", size = 3) +  
  geom_errorbarh(aes(xmin = Importance - StDev, xmax = Importance + StDev), 
                 height = 0.2, color = "#747474", linewidth = 0.5) +  
  theme_minimal(base_size = 16) +  
  labs(title = "Feature Importance in temporal model", x = "Importance", y = "Variable") +
  theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 18, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1))
Importance_Mosquito

#When we see negative values for the feature importance, it can happen that the predictions on the shuffled 
#data are more accurate than the real data. This occurs when the feature doesn’t matter, 
#but random chance causes the predictions on the shuffled data to be more accurate. (https://someshfengde.medium.com/machine-learning-explainability-permutation-importance-7a9a69bf5943)




