rm(list = ls(all.names = TRUE))
gc()
library(sf) #simple features
library(dplyr) #For mutating dfs
library(ggplot2)  #Plotting histograms
library(gridExtra) #For subplotting
library(lubridate) #Date transformations
library(tidyr) #data wrangling
library(raster) #Claculating CV Score
library(zoo) #For rollapply()
library(corrplot) #Plotting Correlation Matrix
library(Hmisc) #Calculating correlation matrix
library(MASS) #For negative binomial distribution
library(mgcv) # For gam()
library(lme4) # For glmm
library(patchwork) #stacking subplots

#Load Weather data
#(t = Lufttemp [°C], rh = rel. Feuchte [%], ffx = max. Windgeschwindigkeit [m/s], ws_mean = mittere Windspeed [m/s],
                    #p = NS 24h [mm],  so_h = Sonnenscheindauer [h], cglo_j = Globalstrahlung [J/cm2])
  #daily
  Meteo_d <- read.csv("D:/Masterarbeit/Jupyter/Data/Meteo/Messstationen Tagesdaten v2 Datensatz_20210101_20241231.csv")  
  #hourly (wind only)
  Meteo_h <- read.csv("D:/Masterarbeit/Jupyter/Data/Meteo/Messstationen Stundendaten v2 Datensatz_20210101T0000_20241231T0000.csv")  

#Load validated Mosquito Sightings from Mosquito Alert for Graz (after preprocessing in MosAl_explo_prepro.ipynb)
Mos <- st_read("D:/Masterarbeit/Jupyter/Data/QGIS/Mosquito/Mosquitos_val1_31287.shp")

#Load daily Sampling_Effort (after MosAl_explo_prepro.ipynb)
Sampling_Effort <- read.csv("D:/Masterarbeit/Jupyter/Data/Sampling Effort/Postprocessed/SE_pp.csv")
Sampling_Effort$date <- as.Date(Sampling_Effort$date)
Sampling_Effort$X <- NULL


#Manipulate Meteo data
  #Remove station 16401 as there is no data
  Meteo_d <- Meteo_d[Meteo_d$station != 16401,]
  #Remove _flag attributes
  Meteo_d <- Meteo_d %>% dplyr::select(-ends_with("_flag"))
  #Transform time to date
  Meteo_d <- Meteo_d %>%
    mutate(date = as.Date(sub("T.*","",time))) %>%
    dplyr::select(-time)
  #Remove 2021
  Meteo_d$year <- year(Meteo_d$date)
  Meteo_d <- Meteo_d %>% filter(Meteo_d$year != 2021) %>% dplyr::select(-"year")
  #Drop unwanted columns
  Meteo_d <- Meteo_d %>% 
     dplyr::select(-"ffx", -"cglo_j")
  #Rename columns
  col_names <- c("station","t_max", "t_min", "t_mean", "rh_mean", "ws_mean", "p", "so_h", "date")
  colnames(Meteo_d) <- col_names 
  
  #Replace negative rain values (No precip) with 0
  Meteo_d <- Meteo_d %>%
    mutate(
      p = ifelse(p < 0, 0, p)
    )
    #Calcualte daily means from both stations 
  Meteo_d <- Meteo_d %>% 
              group_by(date)%>%
                summarise(
                  t_max = mean(t_max, na.rm = TRUE),
                  t_min = mean(t_min, na.rm = TRUE),
                  t_mean = mean(t_mean, na.rm = TRUE),
                  rh_mean = mean(rh_mean, na.rm = TRUE),
                  ws_mean = mean(ws_mean, na.rm = TRUE),
                  p = mean(p, na.rm = TRUE), 
                  so_h = mean(so_h, na.rm = TRUE),
                )
  
#Manipulation of Mosquito data
  Mos <- Mos %>%
    mutate(creation_d = ymd(creation_d))
  
  #Clean up df
    Mos <- Mos %>% dplyr::select("creation_d","creation_y", "creation_m")
  #Remove 2021 due to sparsity of data
    Mos_22_24 <- Mos %>%
      filter(creation_y != 2021)
  
  # Prepare daily counts of mosquitoes and create month and year attribute
    #Create sequence of days
      date_range <- seq(min(Meteo_d$date), max(Meteo_d$date), by = "day")
    #Df from sequence
    full_dates <- data.frame(creation_d = date_range)
   
 
    daily_data <- Mos_22_24 %>%
      count(creation_d) %>% #aggregate daily counts
      right_join(full_dates, by = "creation_d") %>% # Right join to keep all days
      mutate(
        creation_m = month(creation_d),   # create month and year attribute
        creation_y = year(creation_d))%>%
      rename(count = n, date = creation_d)%>% #rename columns
      replace_na(list(count = 0)) #replace NAN with 0
  
    #drop geometry column
    daily_data <- st_drop_geometry(daily_data)
    daily_data <- daily_data %>% filter(daily_data$creation_y != 2021)

    # transfrom factor for second y-axis
    transform_factor <- max(daily_data$count)/max(Sampling_Effort$SE)
    
    #Plot SE and daily counts for 2022-2024
    ggplot(daily_data, aes(x = date, y = count)) +
      geom_col(fill = "orange") +
      geom_line(data = Sampling_Effort, aes(x = date, y = SE * transform_factor), color = "blue", alpha = 0.4) +
      labs(title = "Daily Mosquito Sightings between 2022 and 2024", x = "Year", y = "Daily observations") +
      theme_minimal() +
      scale_y_continuous(
        sec.axis = sec_axis(~ . / transform_factor, name = "Sampling_Effort")
         ) +
      scale_x_date(date_breaks = "1 year", date_labels = "%Y", date_minor_breaks = "1 month")  
    
    #Save plot
    ggsave(plot = last_plot(), filename = "D:/Masterarbeit/Figures/Daily_sightings_22_24.png")
    
#Calculate mean and variation of hourly windspeed
    #Hourly means between stations
    daily_ws_mean <- Meteo_h %>%
      group_by(time) %>% # Group by hour
      summarise(
        ws_mean = mean(ff, na.rm = TRUE)
      ) %>%
      mutate(
        date = date(time)
      )
    
    #CV Score per day                                   
    daily_ws_cv <- daily_ws_mean %>%
      group_by(date) %>%
      summarise(
        ws_cv = raster::cv(ws_mean)
      )
   #Join daily df with ws_CV 
   Meteo_d <- left_join(Meteo_d, daily_ws_cv, by = "date")
   
   # Claculate Meteo parameters across 14 and 28 days in advance
    # Compute rolling mean, standard deviation and coefficient of variation for each parameter
   Meteo_d <- Meteo_d %>%
     arrange(date) %>%
     mutate(
       # 14-days time-lagged
       t_max_mean_14 = rollapply(t_max, width = 14, FUN = mean, fill = NA, align = "right"),
       t_max_sd_14 = rollapply(t_max, width = 14, FUN = sd, fill = NA, align = "right"),
       t_max_cv_14 = ifelse(t_max_mean_14 != 0, (t_max_sd_14 / t_max_mean_14) * 100, NA),
       
       t_min_mean_14 = rollapply(t_min, width = 14, FUN = mean, fill = NA, align = "right"),
       t_min_sd_14 = rollapply(t_min, width = 14, FUN = sd, fill = NA, align = "right"),
       t_min_cv_14 = ifelse(t_min_mean_14 != 0, (t_min_sd_14 / t_min_mean_14) * 100, NA),
       
       t_mean_mean_14 = rollapply(t_mean, width = 14, FUN = mean, fill = NA, align = "right"),
       t_mean_sd_14 = rollapply(t_mean, width = 14, FUN = sd, fill = NA, align = "right"),
       t_mean_cv_14 = ifelse(t_mean_mean_14 != 0, (t_mean_sd_14 / t_mean_mean_14) * 100, NA),
       
       rh_mean_mean_14 = rollapply(rh_mean, width = 14, FUN = mean, fill = NA, align = "right"),
       rh_mean_sd_14 = rollapply(rh_mean, width = 14, FUN = sd, fill = NA, align = "right"),
       rh_mean_cv_14 = ifelse(rh_mean_mean_14 != 0, (rh_mean_sd_14 / rh_mean_mean_14) * 100, NA),
       
       ws_mean_mean_14 = rollapply(ws_mean, width = 14, FUN = mean, fill = NA, align = "right"),
       ws_mean_sd_14 = rollapply(ws_mean, width = 14, FUN = sd, fill = NA, align = "right"),
       ws_mean_cv_14 = ifelse(ws_mean_mean_14 != 0, (ws_mean_sd_14 / ws_mean_mean_14) * 100, NA),
       
       p_mean_14 = rollapply(p, width = 14, FUN = mean, fill = NA, align = "right"),
       p_sd_14 = rollapply(p, width = 14, FUN = sd, fill = NA, align = "right"),
       p_cv_14 = ifelse(p_mean_14 != 0, (p_sd_14 / p_mean_14) * 100, NA),
       p_acc_14 = rollapply(p, width = 14, FUN = sum, fill = NA, align = "right"),
       
       so_h_mean_14 = rollapply(so_h, width = 14, FUN = mean, fill = NA, align = "right"),
       so_h_sd_14 = rollapply(so_h, width = 14, FUN = sd, fill = NA, align = "right"),
       so_h_cv_14 = ifelse(so_h_mean_14 != 0, (so_h_sd_14 / so_h_mean_14) * 100, NA),
       )
   
   # 28-days time-lagged
   Meteo_d <- Meteo_d %>%
     arrange(date) %>%
     mutate(
       t_max_mean_28 = rollapply(t_max, width = 28, FUN = mean, fill = NA, align = "right"),
       t_max_sd_28 = rollapply(t_max, width = 28, FUN = sd, fill = NA, align = "right"),
       t_max_cv_28 = ifelse(t_max_mean_28 != 0, (t_max_sd_28 / t_max_mean_28) * 100, NA),
       
       t_min_mean_28 = rollapply(t_min, width = 28, FUN = mean, fill = NA, align = "right"),
       t_min_sd_28 = rollapply(t_min, width = 28, FUN = sd, fill = NA, align = "right"),
       t_min_cv_28 = ifelse(t_min_mean_28 != 0, (t_min_sd_28 / t_min_mean_28) * 100, NA),
       
       t_mean_mean_28 = rollapply(t_mean, width = 28, FUN = mean, fill = NA, align = "right"),
       t_mean_sd_28 = rollapply(t_mean, width = 28, FUN = sd, fill = NA, align = "right"),
       t_mean_cv_28 = ifelse(t_mean_mean_28 != 0, (t_mean_sd_28 / t_mean_mean_28) * 100, NA),
       
       rh_mean_mean_28 = rollapply(rh_mean, width = 28, FUN = mean, fill = NA, align = "right"),
       rh_mean_sd_28 = rollapply(rh_mean, width = 28, FUN = sd, fill = NA, align = "right"),
       rh_mean_cv_28 = ifelse(rh_mean_mean_28 != 0, (rh_mean_sd_28 / rh_mean_mean_28) * 100, NA),
       
       ws_mean_mean_28 = rollapply(ws_mean, width = 28, FUN = mean, fill = NA, align = "right"),
       ws_mean_sd_28 = rollapply(ws_mean, width = 28, FUN = sd, fill = NA, align = "right"),
       ws_mean_cv_28 = ifelse(ws_mean_mean_28 != 0, (ws_mean_sd_28 / ws_mean_mean_28) * 100, NA),
       
       p_mean_28 = rollapply(p, width = 28, FUN = mean, fill = NA, align = "right"),
       p_sd_28 = rollapply(p, width = 28, FUN = sd, fill = NA, align = "right"),
       p_cv_28 = ifelse(p_mean_28 != 0, (p_sd_28 / p_mean_28) * 100, NA),
       p_acc_28 = rollapply(p, width = 28, FUN = sum, fill = NA, align = "right"),
       
       so_h_mean_28 = rollapply(so_h, width = 28, FUN = mean, fill = NA, align = "right"),
       so_h_sd_28 = rollapply(so_h, width = 28, FUN = sd, fill = NA, align = "right"),
       so_h_cv_28 = ifelse(so_h_mean_28 != 0, (so_h_sd_28 / so_h_mean_28) * 100, NA),
     )
 #Clean df
   Meteo_d <- Meteo_d %>%
     dplyr::select(-contains("_sd"), -"p_mean_28", -"p_mean_14")
 
 #Join Mos_Count with Meteo_d and sampling effort
   df_train <- left_join(daily_data, Meteo_d, by = "date")
   df_train <- left_join(df_train, Sampling_Effort, by = "date")
   #Create day of the year attribute
   df_train$doy <- yday(df_train$date)
   

#Save Meteo_data
   write.csv(Meteo_d, "D:/Masterarbeit/Jupyter/Data/Meteo/Postprocessed/Meteo_d.csv",  row.names = FALSE)
#Save df_train_meteo_MA  
   write.csv(df_train, "D:/Masterarbeit/Jupyter/Data/QGIS/Sampling_data/Meteo/df_train_meteo_MA.csv",
             row.names = FALSE)
 
###########################
### Correlation matrix ####
###########################
    
   vars_eda <- c(
     "SE",
     "t_mean",
     "t_max",
     "t_min",
     "rh_mean",
     "ws_mean",
     "ws_cv",
     "p",
     "so_h",
     "so_h_mean_28",
     "so_h_mean_14",
     "t_mean_mean_28",
     "t_mean_mean_14",
     "t_max_mean_28",
     "t_max_mean_14",
     "t_min_mean_28",
     "t_min_mean_14",
     "p_acc_28",
     "p_acc_14",
     "ws_mean_mean_28",
     "ws_mean_mean_14",
     "rh_mean_mean_28",
     "rh_mean_mean_14"
   )
   
    ## Poisson or negative binomial distribution? 
    count_mean <- mean(df_train$count, na.rm = TRUE)
    count_variance <- var(df_train$count, na.rm = TRUE)
    
    print(count_mean)
    print(count_variance)
    
    ## Variance exceeds mean with factor 7 --> Negative binomial family 
    
    ###########################
    ### Correlation matrix ####
    ###########################
    
    vars_eda <- c(
      "SE",
      "t_mean",
      "t_max",
      "t_min",
      "rh_mean",
      "ws_mean",
      "ws_cv",
      "p",
      "so_h",
      "so_h_mean_28",
      "so_h_mean_14",
      "t_mean_mean_28",
      "t_mean_mean_14",
      "t_max_mean_28",
      "t_max_mean_14",
      "t_min_mean_28",
      "t_min_mean_14",
      "p_acc_28",
      "p_acc_14",
      "ws_mean_mean_28",
      "ws_mean_mean_14",
      "rh_mean_mean_28",
      "rh_mean_mean_14"
    )
    #Correlation of weather parameters only during active season 
    eda_cor <- df_train %>%
      dplyr::filter(doy >= 121 & doy <= 304) %>%
      dplyr::select(count, all_of(vars_eda))
    
    #Pearson-correlations of each attribute in vars_eda with count
    cor_values <- sapply(eda_cor, function(x) cor(x, eda_cor$count, use = "complete.obs")) #ignores missing values 
    #Create df for readability
    cor_df <- data.frame(Predictor = names(cor_values), Correlation = cor_values)
    
    # Sorting in descending order
    sorted_cor_df <- cor_df %>%
      arrange(desc(Correlation))
    
    print(sorted_cor_df)
    
    ### Correlation Matrix ###
    
    vars_mat <- c(
      "count",
      "SE",
      "t_max",
      "rh_mean",
      "ws_mean",
      "ws_cv",
      "p",
      "so_h",
      "so_h_mean_28",
      "t_mean_mean_28",
      "p_acc_28",
      "ws_mean_mean_28",
      "rh_mean_mean_28"
    )
    
    #Select data for matrix
    eda_mat <- eda_cor[, vars_mat]
    
    #Create matrix
    cormat <- rcorr(as.matrix(eda_mat))
    corrplot(cormat$r)
    
    #Plotting
    png("D:/Masterarbeit/Figures/CorMat_Meteo.png")
    corrplot(cormat$r)
    dev.off()
    
   