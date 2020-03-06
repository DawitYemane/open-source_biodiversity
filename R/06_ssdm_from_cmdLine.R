# -------------- This will apply ensemble SDM (Species distribution Model) for each species and the three functional groups --------
# ------------- load libraries to be used for data reading, processing, and preparation ----------
libs <-c('tidyverse','sf','raster','SSDM','parallel');lapply(libs, library,character.only=T)

# ------------------ load the species occurrence/encounter and environmental data --------------------

env_rast <-load_var(path = proc_env_Dir,format='.grd',tmp = TRUE,Norm = FALSE,categorical = 'bot_sediment')

rast_vals <- values(env_rast)%>%as_tibble()
rast_vif_al<-car::vif(lm(bot_tempmean~bot_carbonphytomean+bot_chlomean+bot_depth+bot_dissoxmean+
                           bot_ironmean+bot_lightbotmean+bot_nitratemean+bot_phosphatemean+bot_ppmean+
                           bot_salinitymean+bot_sediment+bot_silicatemean+bot_TRI+sur_carbonphytomean+
                           sur_chlomean+sur_dissoxmean+sur_ironmean+sur_nitratemean+sur_phosphatemean+
                           sur_ppmean+sur_salinitymean+sur_silicatemean+sur_tempmean,data=rast_vals))
rast_vif_sel<-car::vif(lm(bot_depth~bot_dissoxmean+
                            bot_ironmean+bot_lightbotmean+bot_ppmean+
                            bot_salinitymean+bot_TRI+
                            sur_chlomean+
                            sur_salinitymean+sur_silicatemean+sur_tempmean,data=rast_vals))

excl_ind <- which(!names(env_rast)%in%c(names(rast_vif_sel),'bot_depth'))
env_rast_sel <-dropLayer(env_rast,excl_ind)

occ_env_fish_all <-load_occ(path=proc_obs_Dir,Xcol = 'long',Ycol = 'lat',file='processed_occ_data_fish_all.csv',
                            Spcol = 'scientificName',Env = env_rast_sel,GeoRes = FALSE,sep=',')
occ_env_benthos_all <-load_occ(path=proc_obs_Dir,Xcol = 'long',Ycol = 'lat',file='processed_occ_data_benthos_all.csv',
                               Spcol = 'scientificName',Env = env_rast_sel,GeoRes = FALSE,sep=',')
occ_env_zoo_all <-load_occ(path=proc_obs_Dir,Xcol = 'long',Ycol = 'lat',file='processed_occ_data_zooplankton_all.csv',
                           Spcol = 'scientificName',Env = env_rast_sel,GeoRes = FALSE,sep=',')


# -------------- To speedup the process parallel computation was applied -------------------------

#  ----  Utility function for parallel SDM -------------------------------


# ----------- Prepare data for parallel computation -----------------------


# ----------   Parallel ensemble SDM: benthos -----------------------------


# -------------- Parallel ensemble SDM: fish  -----------------------------


# ------------------ Parallel ensemble SDM: zooplankton -------------------









