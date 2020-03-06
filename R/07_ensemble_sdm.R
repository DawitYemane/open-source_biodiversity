# ---- Processs  and prepare the outputs from ensemble distriution maps ---------------------------
# ------ load libraries to be used for processing of distribution maps --------------------
libs <-c('tidyverse','sf','raster','SSDM','parallel','here')
lapply(libs, library,character.only=T)

# R -e 'library(tidyverse);library(sf);library(raster);library(SSDM);library(parallel);source('R/07_prepare_biodiversity_maps.R')  
# ----- load environmental and biological data ----------------------------
env_rast <-load_var(path = 'data/tidy/env/sa/',format='.grd',tmp = TRUE,Norm = FALSE,categorical = 'bot_sediment')


# ------------------------check for multicollinearity and exclude highly correlat ed variables -------------------
rast_vals <- raster::values(env_rast)%>%as_tibble()
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


# ------------- retain only least correlated variables --------------------
excl_ind <- which(!names(env_rast)%in%c(names(rast_vif_sel),'bot_depth'))
env_rast_sel <-dropLayer(env_rast,excl_ind)


# -------- write  selected/retained environmental variables to file -------
raster::writeRaster(x = env_rast_sel,filename = here('data/sel_env/sa/.grd'),
                    suffix=names(env_rast_sel), bylayer=TRUE,overwrite=TRUE)



# -------- load occurrence  data for all three groups ---------------------
occ_env_fish_all <-load_occ(path='data/tidy/bio/',Xcol = 'long',Ycol = 'lat',file='processed_occ_data_fish_all.csv',
                            Spcol = 'scientificName',Env = env_rast_sel,GeoRes = FALSE,sep=',')
occ_env_benthos_all <-load_occ(path='data/tidy/bio/',Xcol = 'long',Ycol = 'lat',file='processed_occ_data_benthos_all.csv',
                               Spcol = 'scientificName',Env = env_rast_sel,GeoRes = FALSE,sep=',')
occ_env_zoo_all <-load_occ(path='data/tidy/bio/',Xcol = 'long',Ycol = 'lat',file='processed_occ_data_zooplankton_all.csv',
                           Spcol = 'scientificName',Env = env_rast_sel,GeoRes = FALSE,sep=',')


# -----   names of SDM  that can be considered ----------------------------
mods<-c('CTA', 'SVM','GLM','GAM','MARS','GBM','RF','ANN','MAXENT')

# -----  utility functions to be used  for parallel fitting of SDMs (multiple species at once) ---------------
par_ens_spp_data <- function(species,occ_data){
  library(SSDM)
  name_sdm <- species
  Spoccurrences <- subset(occ_data, occ_data[which(names(occ_data) == 'scientificName')] == species)
  
  #tmppath <-modelEns
  
  cat("Ensemble modelling :", name_sdm, "\n\n")
  

# --- try-catch block just in case ----------------------------------------
  enm <- try(ensemble_modelling(algorithms = mods[-c(2,4,6,9)], Spoccurrences, 
                                Env = env_rast_sel, Xcol='long', Ycol='lat', rep = 3, name = name_sdm, 
                                save = TRUE, path = temp_modelEns,trees=500,ensemble.metric = 'AUC',
                                ensemble.thresh = 0.7,
                                tmp = FALSE, verbose = TRUE),silent = TRUE)
  
  if (inherits(enm, "try-error")) {
    
    enm <- NULL
  }
  else {
    if (!is.null(enm)) {
      tSp<-ifelse(stringr::str_detect(species,'[()]'),stringr::str_replace_all(species,pattern = '[()]',''),species)
      tSp2 <-ifelse(stringr::str_detect(tSp,' '),stringr::str_replace_all(tSp,' ','_'),tSp)
      assign(tSp2,enm)
      save(list=tSp2,file = paste(temp_modelEns,species,'.RData',sep=''))
      
    }
    
  }
  rm(enm)
  gc()
}


#  wrapper funcion to run  the above function in parallel -----------------
run_par_ens_data <-function(cl,species,the_data){
  
  return(parLapplyLB(cl,species,function(x) {
    par_ens_spp_data(x,occ_data=the_data)}))
  
}




# ----- data  preparation for parallel computation ------------------------

# ------- list of species -------------------------------------------------
spp_benthos_all <-unique(as.character(occ_env_benthos_all$scientificName))
spp_fish_all <-unique(as.character(occ_env_fish_all$scientificName))
spp_zoo_all <-unique(as.character(occ_env_zoo_all$scientificName))


#  ----  data frame for indexing species for parallel computation ---------
sp_tbl_benthos_all <-tibble(species=spp_benthos_all,ind=round(seq(1,15,length.out = length(spp_benthos_all))))
sp_tbl_fish_all <-tibble(species=spp_fish_all,ind=round(seq(1,20,length.out = length(spp_fish_all))))
sp_tbl_zoo_all <-tibble(species=spp_zoo_all,ind=round(seq(1,5,length.out = length(spp_zoo_all))))


#  --- custom function to  create directories  to store ensemble SDM ------------------
chk_dir_create <-function(tPath){if(dir.exists(tPath)){tPath}else{dir.create(tPath);tPath}}


#  ---  check  create  paths to store results -----------------------------
path_benthos_all <- here::here('analysis/ensemble_out_sa/all_species_ens/ens_benthos_all/')%>%chk_dir_create(.)
path_fish_all <- here::here('analysis/ensemble_out_sa/all_species_ens/ens_fish_all/')%>%chk_dir_create(.)
path_zoo_all <- here::here('analysis/ensemble_out_sa/all_species_ens/ens_zoo_all/')%>%chk_dir_create(.)


# ------ check list of species for which an ensemble sdm was fitted and written to a file ------------------
list_mod_benthos <-lapply(list.files(path_benthos_all)[str_detect(list.files(path_benthos_all),'.RData')],
                          function(x) str_replace_all(string = x,pattern = '.RData',''))%>%unlist()
list_mod_zoo <-lapply(list.files(path_zoo_all)[str_detect(list.files(path_zoo_all),'.RData')],
                      function(x) str_replace_all(string = x,pattern = '.RData',''))%>%unlist()
list_mod_fish <-lapply(list.files(path_fish_all)[str_detect(list.files(path_fish_all),'.RData')],
                       function(x) str_replace_all(string = x,pattern = '.RData',''))%>%unlist()

# ----------- Parallel computation all three groups -----------------------
# -------   Benthos -------------------------------------------------------
if(!(all(list_mod_benthos%in%sp_tbl_benthos_all$species)&!is.null(list_mod_benthos))){
  #### ensemble modelling benthos: all
  for(i in unique(sp_tbl_benthos_all$ind)){
    species2 <-sp_tbl_benthos_all$species[sp_tbl_benthos_all$ind==i]
    temp_modelEns <-path_benthos_all
    #the_data = occ_env_al_8690
    obs_to_export_temp <-c('occ_env_benthos_all','species2','mods','env_rast_sel','temp_modelEns','par_ens_spp_data')
    cl <- makeCluster(6, outfile = "")
    # clusterEvalQ(cl,{library(SSDM)})
    clusterExport(cl, varlist = obs_to_export_temp,envir = environment())
    
    run_par_ens_data(cl,species2,the_data = occ_env_benthos_all)
    gc()
    stopCluster(cl)
    Sys.sleep(2)
  }
  
}

Sys.sleep(10)

# --------- Zooplankton ---------------------------------------------------
if(!(all(list_mod_zoo%in%sp_tbl_zoo_all$species)&!is.null(list_mod_zoo))){
  for(i in unique(sp_tbl_zoo_all$ind)){
    species2 <-sp_tbl_zoo_all$species[sp_tbl_zoo_all$ind==i]
    temp_modelEns <-path_zoo_all
    #the_data = occ_env_al_8690
    obs_to_export_temp <-c('occ_env_zoo_all','species2','mods','env_rast_sel','temp_modelEns','par_ens_spp_data')
    cl <- makeCluster(6, outfile = "")
    # clusterEvalQ(cl,{library(SSDM)})
    clusterExport(cl, varlist = obs_to_export_temp,envir = environment())
    
    run_par_ens_data(cl,species2,the_data = occ_env_zoo_all)
    gc()
    stopCluster(cl)
    Sys.sleep(2)
  }
  
}

Sys.sleep(10)


# ------------ Fishes -----------------------------------------------------
if(!(all(list_mod_fish%in%sp_tbl_fish_all$species)&!is.null(list_mod_fish))){
  for(i in unique(sp_tbl_fish_all$ind)){
    species2 <-sp_tbl_fish_all$species[sp_tbl_fish_all$ind==i]
    temp_modelEns <-path_fish_all
    #the_data = occ_env_al_8690
    obs_to_export_temp <-c('occ_env_fish_all','species2','mods','env_rast_sel','temp_modelEns','par_ens_spp_data')
    cl <- makeCluster(6, outfile = "")
    # clusterEvalQ(cl,{library(SSDM)})
    clusterExport(cl, varlist = obs_to_export_temp,envir = environment())
    
    run_par_ens_data(cl,species2,the_data = occ_env_fish_all)
    gc()
    stopCluster(cl)
    Sys.sleep(2)
  }
  
}
