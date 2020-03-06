#  -----  To be  used to prepare biodiversity maps ------------------------

#  -------------- load  required libraries --------------------------------
libs <- c('SSDM', 'tidyverse','raster','sf','rworldmap', 'rworldxtra', 'rgdal');lapply(libs,library, character.only=TRUE)

#   ------ prepare files and data for biodiversity construction -----------

# ----- taxonomic details  of all fish groups -----------------------------
taxon_info <-read_csv('data/tidy/bio/processed_occ_data_fish_all_withTaxonInfo.csv')

#  ----- subset to three major fish classes -------------------------------
elasmos <- taxon_info%>%filter(class%in%c('Elasmobranchii'))
Cephalos <- taxon_info%>%filter(class%in%c('Cephalopoda'))
Actinopetr <- taxon_info%>%filter(class%in%c('Actinopterygii'))

#  ---- locations of the three major groups -------------------------------
path_benthos_all='analysis/ensemble_out_sa/all_species_ens/ens_benthos_all/'
path_fish_all='analysis/ensemble_out_sa/all_species_ens/ens_fish_all/'
path_zoo_all='analysis/ensemble_out_sa/all_species_ens/ens_zoo_all/'


# --------- list of sdm objects for benthos -------------------------------
stk_bnt_al<-list.files(path_benthos_all)[str_detect(list.files(path_benthos_all),'.RData')]


# ----- list of sdm objects for fishes ------------------------------------
stk_fish_al<-list.files(path_fish_all)[str_detect(list.files(path_fish_all),'.RData')]
stk_elasmo <- stk_fish_al[str_replace_all(stk_fish_al,'.RData','')%in%elasmos$scientificName]
stk_cephalo <- stk_fish_al[str_replace_all(stk_fish_al,'.RData','')%in%Cephalos$scientificName]
stk_actino <- stk_fish_al[str_replace_all(stk_fish_al,'.RData','')%in%Actinopetr$scientificName]

# ------- list of sdm objects for zooplanktons ----------------------------
stk_zoo_al<-list.files(path_zoo_all)[str_detect(list.files(path_zoo_all),'.RData')]

# ---------- load all sdm objects -----------------------------------------
stk_bnt_lst_al <-sapply(paste(path_benthos_all,stk_bnt_al,sep = ''), 
                        function(x) mget(load(x)),simplify = TRUE)

stk_fish_lst_al <-sapply(paste(path_fish_all,stk_fish_al,sep = ''), 
                         function(x) mget(load(x)),simplify = TRUE)

stk_elasmo_lst_al <-sapply(paste(path_fish_all,stk_elasmo,sep = ''), 
                           function(x) mget(load(x)),simplify = TRUE)
stk_cephalo_lst_al <-sapply(paste(path_fish_all,stk_cephalo,sep = ''), 
                            function(x) mget(load(x)),simplify = TRUE)
stk_actino_lst_al <-sapply(paste(path_fish_all,stk_actino,sep = ''), 
                           function(x) mget(load(x)),simplify = TRUE)

stk_zoo_lst_al <-sapply(paste(path_zoo_all,stk_zoo_al,sep = ''), 
                        function(x) mget(load(x)),simplify = TRUE)

# --- source file to be used for  the construction of biodiversity maps --------------------
source('R/custom_biodiveristy.R')


# ----------  stacking: bssdm ---------------------------------------------
bnt_biod_al_bssdm <- stacked_biodiversity(sdm_list = stk_bnt_lst_al,method = 'bSSDM')

fish_biod_al_bssdm <- stacked_biodiversity(sdm_list = stk_fish_lst_al,method = 'bSSDM')

zoo_biod_al_bssdm <- stacked_biodiversity(sdm_list = stk_zoo_lst_al,method = 'bSSDM')


# ------------ stacking: pssdm --------------------------------------------
bnt_biod_al_pssdm <- stacked_biodiversity(sdm_list = stk_bnt_lst_al,method = 'pSSDM')

fish_biod_al_pssdm <- stacked_biodiversity(sdm_list = stk_fish_lst_al,method = 'pSSDM')

zoo_biod_al_pssdm <- stacked_biodiversity(sdm_list = stk_zoo_lst_al,method = 'pSSDM')


#  -------- stacking: maximumLikelIhood -----------------------------------
bnt_biod_al_ml <- stacked_biodiversity(sdm_list = stk_bnt_lst_al,method = 'MaximumLikelihood')

fish_biod_al_ml <- stacked_biodiversity(sdm_list = stk_fish_lst_al,method = 'MaximumLikelihood')
celphao_biod_al_ml <- stacked_biodiversity(sdm_list = stk_cephalo_lst_al,method = 'MaximumLikelihood')
elasmo_biod_al_ml <- stacked_biodiversity(sdm_list = stk_elasmo_lst_al,method = 'MaximumLikelihood')
actino_biod_al_ml <- stacked_biodiversity(sdm_list = stk_actino_lst_al,method = 'MaximumLikelihood')


zoo_biod_al_ml <- stacked_biodiversity(sdm_list = stk_zoo_lst_al,method = 'MaximumLikelihood')

# ------ check plots of biodiversity maps from the three methods of stacking -------

# ---- plot benthos -------------------------------------------------------
bnt_biod_al_bssdm$diversity.map%>%plot()

bnt_biod_al_pssdm$diversity.map%>%plot()

bnt_biod_al_ml$diversity.map%>%plot()


# ---- plot fishes --------------------------------------------------------
fish_biod_al_bssdm$diversity.map%>%plot()

fish_biod_al_pssdm$diversity.map%>%plot()

fish_biod_al_ml$diversity.map%>%plot()

celphao_biod_al_ml$diversity.map%>%plot()
elasmo_biod_al_ml$diversity.map%>%plot()
actino_biod_al_ml$diversity.map%>%plot()

#  ------ plot zooplankton ------------------------------------------------
zoo_biod_al_bssdm$diversity.map%>%plot()

zoo_biod_al_pssdm$diversity.map%>%plot()

zoo_biod_al_ml$diversity.map%>%plot()

# -------- prepare base maps and eez shape  file --------------------------
world <- getMap(resolution = "high")
### convert the spatialPolygonsDataFrame to an sf object
world_sf <- st_as_sf(world)
sa_shape <- world_sf%>%dplyr::filter(NAME=='South Africa')%>%st_crop(xmin=15,xmax=35,ymin=-37,ymax=-20)%>%
  group_by(NAME)%>%summarise()

# --- eez shape file ------------------------------------------------------
shp_eezWAfr = readOGR(dsn = 'spatial/shape/',layer = 'eez_200NM_shape_wAfr')
eez_wAfr_sf <- st_as_sf(shp_eezWAfr)
eez_sAfr_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa')
eez_sAfr_sp <- as_Spatial(eez_sAfr_sf)

eez_sAfr_sf2 <-eez_wAfr_sf%>%filter(Territory1=='South Africa')%>%group_by(Territory1)%>%summarise()

eez_bclme_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa'|Territory1=='Democratic Republic of the Congo'|
                                       Territory1=='Namibia'|Territory1=='Angola')%>%mutate(LME='bclme')%>%
  st_crop(xmin=5,xmax=20,ymin=-40,ymax=-3)%>%group_by(LME)%>%summarise()


# --- prepare for plotting the biodiversity maps from different me --------
all_biodiversity_ml <-cbind(as_tibble(mask(fish_biod_al_ml$diversity.map,eez_sAfr_sp)%>%rasterToPoints())%>%dplyr::rename(fish_al=layer),
                            mask(bnt_biod_al_ml$diversity.map,eez_sAfr_sp)%>%rasterToPoints()%>%as_tibble()%>%dplyr::select(layer)%>%dplyr::rename(benthos_al=layer),
                            mask(zoo_biod_al_ml$diversity.map,eez_sAfr_sp)%>%rasterToPoints()%>%as_tibble()%>%dplyr::select(layer)%>%dplyr::rename(zooplankton_al=layer))
class_biodiversity_ml <-cbind(as_tibble(mask(celphao_biod_al_ml$diversity.map,eez_sAfr_sp)%>%rasterToPoints())%>%dplyr::rename(Cephalopoda=layer),
                              mask(elasmo_biod_al_ml$diversity.map,eez_sAfr_sp)%>%rasterToPoints()%>%as_tibble()%>%dplyr::select(layer)%>%dplyr::rename(Elasmobranchii=layer),
                              mask(actino_biod_al_ml$diversity.map,eez_sAfr_sp)%>%rasterToPoints()%>%as_tibble()%>%dplyr::select(layer)%>%dplyr::rename(Actinopterygii=layer))

all_biodiversity_bssdm <-cbind(as_tibble(mask(fish_biod_al_bssdm$diversity.map,eez_sAfr_sp)%>%rasterToPoints())%>%dplyr::rename(fish_al=layer),
                               mask(bnt_biod_al_bssdm$diversity.map,eez_sAfr_sp)%>%rasterToPoints()%>%as_tibble()%>%dplyr::select(layer)%>%dplyr::rename(benthos_al=layer),
                               mask(zoo_biod_al_bssdm$diversity.map,eez_sAfr_sp)%>%rasterToPoints()%>%as_tibble()%>%dplyr::select(layer)%>%dplyr::rename(zooplankton_al=layer))

all_biodiversity_pssdm <-cbind(as_tibble(mask(fish_biod_al_pssdm$diversity.map,eez_sAfr_sp)%>%rasterToPoints())%>%dplyr::rename(fish_al=layer),
                               mask(bnt_biod_al_pssdm$diversity.map,eez_sAfr_sp)%>%rasterToPoints()%>%as_tibble()%>%dplyr::select(layer)%>%dplyr::rename(benthos_al=layer),
                               mask(zoo_biod_al_pssdm$diversity.map,eez_sAfr_sp)%>%rasterToPoints()%>%as_tibble()%>%dplyr::select(layer)%>%dplyr::rename(zooplankton_al=layer))


#   ------ further processing of the biodiversity data --------------------
all_biodiversity_long_ml <-all_biodiversity_ml%>%group_by(x,y)%>%
  pivot_longer(names_to = 'group_vars',values_to = 'richness',cols = fish_al:zooplankton_al)%>%
  mutate(func_grps = str_replace_all(group_vars,pattern = "_[a-z:0-9]+",'')%>%str_to_title())%>%
  ungroup()

class_biodiversity_long_ml <-class_biodiversity_ml%>%group_by(x,y)%>%
  pivot_longer(names_to = 'group_vars',values_to =  'richness',cols=Cephalopoda:Actinopterygii)%>%
  mutate(func_grps = str_replace_all(group_vars,pattern = "_[a-z:0-9]+",'')%>%str_to_title())%>%
  ungroup()


all_biodiversity_long_bssdm <-all_biodiversity_bssdm%>%group_by(x,y)%>%
  pivot_longer(names_to = 'group_vars',values_to = 'richness',cols=fish_al:zooplankton_al)%>%
  mutate(func_grps = str_replace_all(group_vars,pattern = "_[a-z:0-9]+",'')%>%str_to_title())%>%
  ungroup()

all_biodiversity_long_pssdm <-all_biodiversity_pssdm%>%group_by(x,y)%>%
  pivot_longer(names_to =  'group_vars',values_to = 'richness',cols = fish_al:zooplankton_al)%>%
  mutate(func_grps = str_replace_all(group_vars,pattern = "_[a-z:0-9]+",'')%>%str_to_title())%>%
  ungroup()


#  ---- biodiversity data into raster -------------------------------------
all_biodiversity_ml_zoo <-all_biodiversity_long_ml%>%filter(func_grps=='Zooplankton')%>%
  dplyr::select(x,y,richness)%>%dplyr::rename(zooplankton=richness)%>%as.matrix()
all_biodiversity_ml_benthos <-all_biodiversity_long_ml%>%filter(func_grps=='Benthos')%>%
  dplyr::select(x,y,richness)%>%dplyr::rename(benthos=richness)%>%as.matrix()
all_biodiversity_ml_fish <-all_biodiversity_long_ml%>%filter(func_grps=='Fish')%>%
  dplyr::select(x,y,richness)%>%dplyr::rename(fish=richness)%>%as.matrix()

fish_rast <-rasterFromXYZ(all_biodiversity_ml_fish,crs = "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
benthos_rast <-rasterFromXYZ(all_biodiversity_ml_benthos,crs = "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
zoo_rast <-rasterFromXYZ(all_biodiversity_ml_zoo,crs = "+init=epsg:4326 +proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")


# ---- diversity raster stack ---------------------------------------------
all_div_rast <-stack(fish_rast,benthos_rast,zoo_rast)
system.time(writeRaster(x = all_div_rast,filename = 'analysis/biodiversity/processed_5nm_SA_Diversity.grd',overwrite=TRUE))

#### write processed richness data for three functional groups ##############################
save(list = c('all_biodiversity_long_ml','all_biodiversity_long_bssdm','all_biodiversity_long_pssdm',
              'class_biodiversity_long_ml'),
     file='analysis/biodiversity/processed_biodiversity_all_methods.RData')


# -------- plots of biodiversity maps -------------------------------------
all_biodiversity_long_ml%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill='gray50')+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()


# ---- plot based on the maixmum likeliohood method -----------------------
fish_al_ml<-all_biodiversity_long_ml%>%filter(func_grps=='Fish')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill='gray50')+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()
benthos_al_ml<-all_biodiversity_long_ml%>%filter(func_grps=='Benthos')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill='gray50')+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()
zoo_al_ml<-all_biodiversity_long_ml%>%filter(func_grps=='Zooplankton')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill='gray50')+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()

gridExtra::grid.arrange(fish_al_ml,benthos_al_ml,zoo_al_ml,nrow=1)

ceph_al_ml<-class_biodiversity_long_ml%>%filter(func_grps=='Cephalopoda')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill='gray50')+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()
elasmo_al_ml<-class_biodiversity_long_ml%>%filter(func_grps=='Elasmobranchii')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill='gray50')+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()
actino_al_ml<-class_biodiversity_long_ml%>%filter(func_grps=='Actinopterygii')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill='gray50')+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()

gridExtra::grid.arrange(ceph_al_ml,elasmo_al_ml,actino_al_ml,nrow=1)

# ---- plot based on bssdm (stacking binaries method) ---------------------
fish_al_bssdm<-all_biodiversity_long_bssdm%>%filter(func_grps=='Fish')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()
benthos_al_bssdm<-all_biodiversity_long_bssdm%>%filter(func_grps=='Benthos')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()
zoo_al_bssdm<-all_biodiversity_long_bssdm%>%filter(func_grps=='Zooplankton')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()

gridExtra::grid.arrange(fish_al_bssdm,benthos_al_bssdm,zoo_al_bssdm,nrow=1)


# ---- plot based on the pssdm (stacking prpbability methods) -------------
fish_al_pssdm<-all_biodiversity_long_pssdm%>%filter(func_grps=='Fish')%>%ggplot()+
  geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()
benthos_al_pssdm<-all_biodiversity_long_pssdm%>%filter(func_grps=='Benthos')%>%ggplot()+geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()
zoo_al_pssdm<-all_biodiversity_long_pssdm%>%filter(func_grps=='Zooplankton')%>%ggplot()+
  geom_tile(aes(x,y,fill=richness))+scale_fill_gradientn(colours = fields::tim.colors(100))+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~func_grps)+theme_bw()

gridExtra::grid.arrange(fish_al_pssdm,benthos_al_pssdm,zoo_al_pssdm,nrow=1)







