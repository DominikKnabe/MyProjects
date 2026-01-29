# Predicting the spatio-temporal presence of *Aedes albopictus*

This repository was used to predict the spatial habitat suitability and the daily activity of the Asian tiger mosquito in Graz, Austria using citizen science data from Mosquito Alert and can be adapted for similar use cases in other study areas. 
Spatial habitat suitability was predicted using a binomial generalized additive model (GAM) based on presence/pseudo-absence data while daily mosquito counts were predicted using a negative binomial generalized additive mixed model (GAMM). To predict spatial habitat suitability two approaches were tested, to account for biased observation data due to heterogeneous sampling effort:  

**(1)** Probability-proportional-to-size (PPS) sampling of pseudo-absences based on sampling effort  
**(2)** Employ numerical sampling effort as predictor during model fitting (zeroed out during prediction)  

Predicted daily mosquito counts were spatially distributed using the predicted habitat suitability scores to visualize the presence of *Aedes albopictus* in space and time.  
To provide an overview of the associated research, the abstract of the thesis is provided below. 

Mosquito Alert observation and sampling effort data is accessible from: https://labs.mosquitoalert.com/metadata_public_portal/README.html.  
Hourly weather data is accessible from: https://data.hub.geosphere.at/dataset/klima-v2-1h  
Daily weather data is accessible from: https://data.hub.geosphere.at/dataset/klima-v2-1d  

Further predictor data (land use, land cover, etc.) was provided by the Magistrat Graz, Stadtvermessungsamt and is not publicly available. 

## File descriptions
**Time_series_Aedes_albopictus_presence.gif:** Visualized time series of the predicted spatio-temporal presence in Graz, Austria  
**MosAl_explo_prepro.ipynb:** Data preparation and exploration of Mosquito Alert data (presences and sampling effort)  
**Spatial_Data_Prep.R:** Preparation of spatial predictor data  
**Spatial_sampling.R:** Sampling of pseudo-absences (two approaches) and eda  
**Spatial_model.R:** Model fitting, predictions and validation (random and spatial cross-validation) of spatial habitat suitability model  
**Temporal_data_prep.R:** Data preparation of temporal (time-lagged) predictors, aggregation of daily mosquito counts, correlation matrix  
**Temporal_model.R:** Temporal model fitting (three different sets of predictors were tested), predictions and random cross-validation  
**Visualization_classification.r:** Risk-classification of spatial and temporal predictions, bivariate mapping, time-series  
**Workflow.png:** Flow chart to describe the individual steps of the spatio-temporal analysis

## Abstract 
The Asian tiger mosquito (Aedes albopictus) is one of the most invasive species in the world and a vector for numerous arboviruses. Since its arrival in Austria in 2012, populations have been growing, raising public health concerns. In this work, the spatio-temporal presence of Aedes albopictus was modeled for Graz, Austria, using citizen science data from Mosquito Alert. A framework of generalized additive (mixed) models (GAM/GAMM) was employed to predict spatial habitat suitability and the daily number of mosquito counts using a set of static (land cover, land use, etc.) and dynamic (weather parameters) predictors. Two approaches were tested to account for observer bias in Mosquito Alert data due to heterogenous sampling effort. The spatial model identified residential areas, allotment gardens and cemeteries as areas at risk, while industrial and rural areas showed low suitability. The presence of green spaces in residential areas was found to have a positive effect on habitat suitability. The temporal model highlighted the seasonal activity of Aedes albopictus, which was primarily driven by temperature with peak mosquito season in late summer. Precipitation was not identified as relevant driver for temporal mosquito activity. The results showed that citizen science data can be used exclusively, when bias-corrected, to yield biologically plausible results to gain insight into invasive mosquito distributions. The results provide a basis for enhance mosquito surveillance and targeted strategies for vector management in Austrian cities.
