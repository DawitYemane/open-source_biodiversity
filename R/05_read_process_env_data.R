#  --------------- Used to process environmental raster layers ------------
# -------- libraries to be used to process environmental data -------------
libs <- c('tidyverse','sf','raster','SSDM','ncdf4','rgdal','rworldmap', 'rworldxtra'); lapply(libs, library, character.only=TRUE)

rasterOptions(maxmemory = 1e08,tmpdir = 'temp_raster/',tmptime = 7,progress = 'text')
# ----- raster file processing related -----------------------------------
layer_infos_b=read_csv('data/raw/env/selected_bottomLayer_bioOracle.csv')
layer_infos_s=read_csv('data/raw/env/selected_surfaceLayer_bioOracle.csv')


#  ------ prepare environmental layers for processing ---------------------
# ---- procss raster layer name -------------------------------------------
N_forR_b= layer_infos_b$layer_code[str_detect(layer_infos_b$layer_code,'min|max')]
N_forMean_b = layer_infos_b$layer_code[!str_detect(layer_infos_b$layer_code,'min|max')]

N_forR_s= layer_infos_s$layer_code[str_detect(layer_infos_s$layer_code,'min|max')]
N_forMean_s = layer_infos_s$layer_code[!str_detect(layer_infos_s$layer_code,'min|max')]

fulName_mean_b <- paste('data/raw/env/bottom/',paste(N_forMean_b,'_lonlat.tif',sep = ''),sep='')%>%as.list()
fulName_mean_s <- paste('data/raw/env/surface/',paste(N_forMean_s,'_lonlat.tif',sep = ''),sep='')%>%as.list()

# ------ load all the environmental layers  (bottom and surface) ----------
system.time(alBot_layer <- stack(fulName_mean_b))
system.time(alSur_layer <- stack(fulName_mean_s))

# ---- high resolution bathymetry and sediment  ---------------------------
system.time(tRast_m<-raster('data/raw/env/GEBCO_2014_2D_-30.0_-42.0_60.0_45.0.nc'))
system.time(sed_dat <- raster('data/raw/env/seabed_lithology_v1.nc'))

####### common projection to use
wgs84_projCRS = CRS("+init=epsg:4326")

# --------- load  EEZ shape files -----------------------------------------
shp_eezWAfr = readOGR(dsn = 'spatial/shape/',layer = 'eez_200NM_shape_wAfr')

# ------ get global base maps  --------------------------------------------
world <- getMap(resolution = "high")
### convert the spatialPolygonsDataFrame to an sf object
world_sf <- st_as_sf(world)


# ------ South Africa basemap ---------------------------------------------
sa_shape_al <- world_sf%>%dplyr::filter(NAME=='South Africa')%>%st_crop(xmin=15,xmax=35,ymin=-37,ymax=-22)%>%
  group_by(NAME)%>%summarise()

# --------- BCLME basemap -------------------------------------------------
bclme_al <- world_sf%>%dplyr::filter(NAME=='South Africa'|NAME=='Namibia'|NAME=='Angola')%>%
  st_crop(xmin=11.5,xmax=38,ymin=-37,ymax=-5)%>%
  group_by(NAME)%>%summarise()

# -----  Process eez shape files ------------------------------------------
##### convert to sf object 
eez_wAfr_sf <- st_as_sf(shp_eezWAfr)
eez_sAfr_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa')
eez_bclme_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa'|Territory1=='Namibia'|Territory1=='Angola'|
                                       Territory1=='Democratic Republic of the Congo')%>%
  mutate(LME = 'BCLME')
## merge the polygons of the eez of the three countries into a single polygon ################
eez_bclme_sf_crp <-eez_bclme_sf%>%st_crop(xmin = 8.2,ymin = -38.175,xmax = 20,ymax = -5.029)
eez_all_bclme<-eez_bclme_sf_crp%>%group_by(LME)%>%summarise()
eez_all_bclme%>%ggplot()+geom_sf()+ 
  geom_sf(data = bclme_al, fill = NA)+theme_bw()
eez_all_bclme_sp <-as_Spatial(eez_all_bclme)

eez_sAfr_sf%>%ggplot()+geom_sf()+ 
  geom_sf(data = sa_shape_al, fill = NA)+theme_bw()
eez_sAfr_sp <- as_Spatial(eez_sAfr_sf)

# --------------------------------- Subsetting environments ---------------
# ---------   subset environmental layers to the EEZ ----------------------

# -------------- crop and aggregate to coarser resolution: bathymetry --------
system.time(tRast_crop <-crop(tRast_m,extent(2,40,-40,-5)))
system.time(tRast_crop_agg <-raster::aggregate(tRast_crop,fact=10,fun=mean))
system.time(ocean_mask_eez <- mask(tRast_crop_agg,eez_all_bclme_sp))
system.time(ocean_mask_eezSA <- raster::crop(tRast_crop_agg,extent(10,40,-40,-26.5))%>%mask(x = .,mask = eez_sAfr_sp))
ocean_mask_eez[ocean_mask_eez>0]=NA
ocean_mask_eezSA[ocean_mask_eezSA>0]=NA

# ----- crop, resample/project sediment to get it to the same resolution as the rest of the data ---------------
system.time(sedRast_crop <-crop(sed_dat,extent(2,40,-40,-5)))
system.time(sedRast_mask_eez <- mask(sedRast_crop,eez_all_bclme_sp))
system.time(sedRast_mask_eezSA <- raster::crop(sed_dat,extent(10,40,-40,-26.5))%>%mask(.,eez_sAfr_sp))

system.time(sedRast_mask_eez_hRes <- projectRaster(from = sedRast_mask_eez,to =ocean_mask_eez,method = 'ngb'))
system.time(sedRast_mask_eezSA_hRes <- projectRaster(from = sedRast_mask_eezSA,to = ocean_mask_eezSA,method = 'ngb'))

# ---- mask both the fine and coarse resolution bathymetry data ---------------------------
system.time(ocean_mask_eez_hRes <-mask(tRast_crop,eez_all_bclme_sp))
system.time(ocean_mask_eezSA_hRes <- raster::crop(tRast_crop,extent(10,40,-40,-26.5))%>%mask(x = .,mask = eez_sAfr_sp))

# ----- crop both  environmental layers to the EEZ: bottom ----------------
system.time(alBot_lCrop <- crop(alBot_layer,extent(2,40,-40,-5)))
system.time(alBot_lMask <- mask(alBot_lCrop,eez_all_bclme_sp))
system.time(alBot_lMask_sa <- crop(alBot_layer,extent(10,40,-40,-26.5))%>%mask(.,eez_sAfr_sp))

# ----- crop both  environmental layers to the EEZ: surface ----------------
system.time(alSur_lCrop <- crop(alSur_layer,extent(2,40,-40,-5)))
system.time(alSur_lMask <- mask(alSur_lCrop,eez_all_bclme_sp))
system.time(alSur_lMask_sa <- crop(alSur_layer,extent(10,40,-40,-26.5))%>%mask(.,eez_sAfr_sp))

# --------  crop and convert to fine scale resolution based on the high res bathymetry data --------------------
system.time(alBot_lMask_hRes <- projectRaster(from = alBot_lMask,to =ocean_mask_eez_hRes,method = 'ngb'))
system.time(alBot_lMask_sa_hRes <- projectRaster(from = alBot_lMask_sa,to = ocean_mask_eezSA_hRes,method = 'ngb'))
# -----------  compute topograpic indices ---------------------------------
topInd_tri_al = terrain(tRast_crop,opts='TRI',unit = 'degrees')%>%raster::aggregate(.,fact=10,fun=mean)
topInd_slp_al = terrain(tRast_crop,opts='slope',unit = 'degrees')%>%raster::aggregate(.,fact=10,fun=mean)
topInd_fldir_al = terrain(tRast_crop,opts='flowdir',unit = 'degrees')%>%raster::aggregate(.,fact=10,fun=mean)
topInd_asp_al = terrain(tRast_crop,opts='aspect',unit = 'degrees')%>%raster::aggregate(.,fact=10,fun=mean)


#  ------------ crop topographic indices to bclme eez ---------------------
topInd_tri_bclme <- mask(topInd_tri_al,mask = eez_all_bclme_sp)
topInd_slp_bclme <- mask(topInd_slp_al,mask = eez_all_bclme_sp)
topInd_asp_bclme <- mask(topInd_asp_al,mask = eez_all_bclme_sp)
topInd_fldir_bclme <- mask(topInd_fldir_al,mask = eez_all_bclme_sp)


# ------------- crop topographic indices to South African EEZ -------------
topInd_tri_sa <- raster::crop(topInd_tri_al,extent(10,40,-40,-26.5))%>%mask(.,mask = eez_sAfr_sp)
topInd_slp_sa <- raster::crop(topInd_slp_al,extent(10,40,-40,-26.5))%>%mask(.,mask = eez_sAfr_sp)
topInd_asp_sa <- raster::crop(topInd_asp_al,extent(10,40,-40,-26.5))%>%mask(.,mask = eez_sAfr_sp)
topInd_fldir_sa <- raster::crop(topInd_fldir_al,extent(10,40,-40,-26.5))%>%mask(.,mask = eez_sAfr_sp)


# ------- create raster stacks  based on topographic indices: bclme and sa --------------
bathy_drvd_phys_bclme <-stack(topInd_asp_bclme,topInd_fldir_bclme,topInd_slp_bclme,topInd_tri_bclme,ocean_mask_eez)
bathy_drvd_phys_sa <-stack(topInd_asp_sa,topInd_fldir_sa,topInd_slp_sa,topInd_tri_sa,ocean_mask_eezSA)
names(bathy_drvd_phys_bclme)<- c('aspect','flowdir','slope','tri','depth')
names(bathy_drvd_phys_sa)<- c('aspect','flowdir','slope','tri','depth')

# ------------------- Combine all raster stacks ---------------------------
names(ocean_mask_eez) = 'depth'
names(ocean_mask_eezSA) = 'depth'
ocean_mask_eez_pos <-calc(ocean_mask_eez,fun = function(x) {x*-1})
ocean_mask_eezSA_pos <-calc(ocean_mask_eezSA,fun = function(x) {x*-1})

sedRast_mask_eez_hRes_f <-as.factor(sedRast_mask_eez_hRes)
sedRast_mask_eezSA_hRes_f <-as.factor(sedRast_mask_eezSA_hRes)


# ------------ unstack the raster stack file and add new rasters:  bclme --------
unst_botLayer = unstack(alBot_lMask) ## first unstack the raster layer 
unst_botLayer[[12]]=sedRast_mask_eez_hRes_f
unst_botLayer[[13]]=topInd_tri_bclme ## add the additional layers
unst_botLayer[[14]]=ocean_mask_eez_pos

index_stack = stack(unst_botLayer)

# ---------------- unstack the raster stack file and add new rasters: sa --------
unst_botLayerSA <- unstack(alBot_lMask_sa)
unst_botLayerSA[[12]] = sedRast_mask_eezSA_hRes_f
unst_botLayerSA[[13]] = topInd_tri_sa
unst_botLayerSA[[14]] = ocean_mask_eezSA_pos
index_stackSA <-stack(unst_botLayerSA)


#  ------------- rename variables in the raster stack ---------------------
tName_b=names(index_stack)
tName_s=names(alSur_lMask)

alName_b=str_split(tName_b[-c(12:14)],'_')%>%lapply(., function(x) x[2])%>%unlist()
alName_s = str_split(tName_s,'_')%>%lapply(.,function(x) x[2])%>%unlist()

alName_b=c(alName_b,'sediment','TRI','depth')

names(index_stack)=alName_b
names(index_stackSA) = alName_b
names(alSur_lMask)=alName_s
names(alSur_lMask_sa)=alName_s

# ------------- Write processed environmental layers to file -------------
writeRaster(index_stackSA,filename = here::here('data/tidy/env/sa/bot.grd'),format='raster',bylayer=TRUE,
            overwrite=TRUE,suffix=names(index_stackSA))
writeRaster(index_stackSA,filename = here::here('data/tidy/env/sa/bot.tif'),format='GTiff',bylayer=TRUE,overwrite=TRUE,
            suffix=names(index_stackSA))
writeRaster(alSur_lMask_sa,filename = here::here('data/tidy/env/sa/sur.grd'),format='raster',bylayer=TRUE,
            overwrite=TRUE,suffix=names(alSur_lMask_sa))
writeRaster(alSur_lMask_sa,filename = here::here('data/tidy/env/sa/bot.tif'),format='GTiff',bylayer=TRUE,overwrite=TRUE,
            suffix=names(alSur_lMask_sa))


# ------------- write bathymetry derived variables  -----------------------
writeRaster(bathy_drvd_phys_sa,filename = here::here('bathy_derived_indices/sa/bathy.grd'),format='raster',bylayer=TRUE,
            overwrite=TRUE,suffix=names(bathy_drvd_phys_sa))
writeRaster(bathy_drvd_phys_sa,filename = here::here('bathy_derived_indices/sa/bathy.tif',sep = ''),format='GTiff',bylayer=TRUE,overwrite=TRUE,
            suffix=names(bathy_drvd_phys_sa))


# -------- write processed  raster stack to a file ------------------------
system.time(writeRaster(x = index_stack,filename = 'data/tidy/env/processed_5nm_bclme_data_bottom.grd',overwrite=TRUE))
system.time(writeRaster(x = index_stackSA,filename = 'data/tidy/env/processed_5nm_SA_data_bottom.grd',overwrite=TRUE))
system.time(writeRaster(x = alSur_lMask,filename = 'data/tidy/env/processed_5nm_bclme_data_surface.grd',overwrite=TRUE))
system.time(writeRaster(x = alSur_lMask_sa,filename = 'data/tidy/env/processed_5nm_SA_data_surface.grd',overwrite=TRUE))


# write processed environmental layers to file: bclme ---------------------
index_stack_Bot_bclme_compNoSed <-dropLayer(index_stack,12)#%>%crop(.,extent(5,20,-40,-5))
alSur_lMask_bclme <-alSur_lMask#,extent(5,20,-40,-5))

writeRaster(index_stack_Bot_bclme_compNoSed,filename = here::here('data/tidy/env/bclme/bot.grd'),format='raster',bylayer=TRUE,
            overwrite=TRUE,suffix=names(index_stack_Bot_bclme_compNoSed))
writeRaster(index_stack_Bot_bclme_compNoSed,filename = here::here('data/tidy/env/bot.tif'),format='GTiff',bylayer=TRUE,overwrite=TRUE,
            suffix=names(index_stack_Bot_bclme_compNoSed))
writeRaster(alSur_lMask_bclme,filename = here::here('data/tidy/env/bclme/sur.grd'),format='raster',bylayer=TRUE,
            overwrite=TRUE,suffix=names(alSur_lMask_bclme))
writeRaster(alSur_lMask_bclme,filename = here::here('data/tidy/env/bclme/sur.tif'),format='GTiff',bylayer=TRUE,overwrite=TRUE,
            suffix=names(alSur_lMask_bclme))
writeRaster(bathy_drvd_phys_bclme,filename = here::here('bathy_derived_indices/bclme/bathy.grd'),format='raster',bylayer=TRUE,
            overwrite=TRUE,suffix=names(bathy_drvd_phys_bclme))
writeRaster(bathy_drvd_phys_bclme,filename = here::here('bathy_derived_indices/bclme/bathy.tif'),format='GTiff',bylayer=TRUE,overwrite=TRUE,
            suffix=names(bathy_drvd_phys_bclme))






