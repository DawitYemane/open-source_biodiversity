#  -----  To  be used to prepare historical distribution maps -------------
#  --- load required libraries --------------------------------------------
libs <-c('tidyverse','sf','raster','SSDM','rworldmap','here','rgdal','rworldxtra')
lapply(libs, library,character.only=T)


#  ------ load and prepare selected environmental raster ------------------
env_rast_sel <- load_var(path = 'data/sel_env/sa/',format = '.grd')
t_names <- str_replace_all(names(env_rast_sel),'X_','')
names(env_rast_sel) <- t_names


# ----  prepare and process environmental grid ----------------------------
bot_depth_rast <- env_rast_sel[[1]]
bot_depth_rast10 <-raster::aggregate(bot_depth_rast,fact=2)

rast_df <- rasterToPoints(bot_depth_rast)
cell_no <- cellFromXY(object = bot_depth_rast,xy = rast_df[c(1:2)])
rast_df <- rast_df%>%as.data.frame()%>%mutate(cell = cell_no)

## 10NM resolution
rast_df10 <- rasterToPoints(bot_depth_rast10)
cell_no10 <- cellFromXY(object = bot_depth_rast10,xy = rast_df10[,c(1:2)])
rast_df10 <- rast_df10%>%as.data.frame()%>%mutate(cell = cell_no10)

#  ----- load all historical occurrence data used in distribution modelling ------------------------------------
# ---  raw occurrence data ------------------------------------------------
occ_benthos <- read_csv('data/raw/bio/benthos_SA.csv')
occ_fish <- read_csv('data/raw/bio/fish_SA.csv')
occ_phytoZoo <- read_csv('data/raw/bio/phytoZooplankton_SA.csv')


#  ------  process  occurrence data ---------------------------------------
benthos_occ <- occ_benthos%>%filter(!is.na(scientificName))%>%
  mutate(Year_grp  = case_when(
    year<=1900 ~ '<= 1900',
    year>1900 &year <=1920 ~ '1900 - 1920',
    year>1920 &year <=1940 ~ '1920 - 1940',
    year>1940 &year <=1960 ~ '1940 - 1960',
    year>1960 &year <=1980 ~ '1960 - 1980',
    year>1980 &year <=2000 ~ '1980 - 2000',
    year>2000 ~ '>= 2000'
  ),Year_grp =reorder(factor(Year_grp),year,mean))%>%dplyr::select(long,lat,Year_grp,scientificName)
fish_occ <- occ_fish%>%filter(!is.na(scientificName))%>%
  mutate(Year_grp  = case_when(
    year<=1900 ~ '<= 1900',
    year>1900 &year <=1920 ~ '1900 - 1920',
    year>1920 &year <=1940 ~ '1920 - 1940',
    year>1940 &year <=1960 ~ '1940 - 1960',
    year>1960 &year <=1980 ~ '1960 - 1980',
    year>1980 &year <=2000 ~ '1980 - 2000',
    year>2000 ~ '>= 2000'
  ),Year_grp =reorder(factor(Year_grp),year,mean))%>%
  dplyr::select(long,lat,Year_grp,scientificName)
zoo_occ <- occ_phytoZoo%>%filter(adult=='zooplankton'&!is.na(scientificName))%>%
  mutate(Year_grp  = case_when(
    year<=1900 ~ '<= 1900',
    year>1900 &year <=1920 ~ '1900 - 1920',
    year>1920 &year <=1940 ~ '1920 - 1940',
    year>1940 &year <=1960 ~ '1940 - 1960',
    year>1960 &year <=1980 ~ '1960 - 1980',
    year>1980 &year <=2000 ~ '1980 - 2000',
    year>2000 ~ '>= 2000'
  ),Year_grp =reorder(factor(Year_grp),year,mean))%>%
  dplyr::select(long,lat,Year_grp,scientificName)


# --- prepare occurrence data at 10NMile resolution -----------------------
benthos_occ10 <-benthos_occ%>%mutate(cell = cellFromXY(object = bot_depth_rast10,as.matrix(benthos_occ[,c(1:2)])))%>%
  dplyr::select(Year_grp,scientificName,cell)%>%distinct()
fish_occ10 <-fish_occ%>%mutate(cell = cellFromXY(object = bot_depth_rast10,as.matrix(fish_occ[,c(1:2)])))%>%
  dplyr::select(Year_grp,scientificName,cell)%>%distinct()
zoo_occ10 <-zoo_occ%>%mutate(cell = cellFromXY(object = bot_depth_rast10,as.matrix(zoo_occ[,c(1:2)])))%>%
  dplyr::select(Year_grp,scientificName,cell)%>%distinct()

benthos_summ_gr10 <- benthos_occ10%>%group_by(Year_grp,cell)%>%count()
fish_summ_gr10 <- fish_occ10%>%group_by(Year_grp,cell)%>%count()
zoo_summ_gr10 <- zoo_occ10%>%group_by(Year_grp,cell)%>%count()


#  ---- join raster cells to occurence data and compute richness ----------
map_rich_benthos_grp10 <- rast_df10%>%left_join(benthos_summ_gr10)%>%filter(!is.na(bot_depth))%>%
  mutate(n_count = ifelse(is.na(n),0,n))%>%mutate(n_grp = case_when(
    n_count==0 ~ '0',
    n_count>0&n_count<=10 ~ '1-10',
    n_count>10&n_count<=20 ~ '11-20',
    n_count>20&n_count<=50 ~ "21-50",
    n_count>50 ~ '> 50'
  ),n_grp = reorder(factor(n_grp),n_count,mean))
map_rich_fish_grp10 <- rast_df10%>%left_join(fish_summ_gr10)%>%filter(!is.na(bot_depth))%>%
  mutate(n_count = ifelse(is.na(n),0,n))%>%mutate(n_grp = case_when(
    n_count==0 ~ '0',
    n_count>0&n_count<=10 ~ '1-10',
    n_count>10&n_count<=20 ~ '11-20',
    n_count>20&n_count<=50 ~ "21-50",
    n_count>50 ~ '> 50'
  ),n_grp = reorder(factor(n_grp),n_count,mean))
map_rich_zoo_grp10 <- rast_df10%>%left_join(zoo_summ_gr10)%>%filter(!is.na(bot_depth))%>%
  mutate(n_count = ifelse(is.na(n),0,n))%>%mutate(n_grp = case_when(
    n_count==0 ~ '0',
    n_count>0&n_count<=10 ~ '1-10',
    n_count>10&n_count<=20 ~ '11-20',
    n_count>20&n_count<=50 ~ "21-50",
    n_count>50 ~ '> 50'
  ),n_grp = reorder(factor(n_grp),n_count,mean))

# -------- prepare base maps and eez shape  file --------------------------
world <- getMap(resolution = "high")
### convert the spatialPolygonsDataFrame to an sf object
world_sf <- st_as_sf(world)
sa_shape <- world_sf%>%dplyr::filter(NAME=='South Africa')%>%st_crop(xmin=15,xmax=35,ymin=-37,ymax=-20)%>%
  group_by(NAME)%>%summarise()

# --- eez shape file ------------------------------------------------------
shp_eezWAfr = readOGR(dsn = 'spatial/shape/',layer = 'eez_200NM_shape_wAfr')
eez_wAfr_sf <- st_as_sf(shp_eezWAfr)
eez_sAfr_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa')%>%group_by(Territory1)%>%summarise()
eez_sAfr_sp <- as_Spatial(eez_sAfr_sf)

eez_bclme_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa'|Territory1=='Democratic Republic of the Congo'|
                                       Territory1=='Namibia'|Territory1=='Angola')%>%mutate(LME='bclme')%>%
  lwgeom::st_make_valid()%>%st_crop(xmin=5,xmax=20,ymin=-40,ymax=-3)%>%group_by(LME)%>%summarise()



# ----- plot historical richness maps -------------------------------------
col_code = c("0"='gray',"1-10"='blue', "11-20"= 'green',"21-50"='orange', "> 50"='red' )

map_rich_benthos_grp10%>%filter(!is.na(Year_grp))%>%ggplot()+geom_tile(aes(x,y,fill=n_grp))+
  scale_fill_manual(values = col_code,name='number of species')+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~Year_grp)+theme_bw()
map_rich_fish_grp10%>%filter(!is.na(Year_grp))%>%ggplot()+geom_tile(aes(x,y,fill=n_grp))+
  scale_fill_manual(values = col_code,name='number of species')+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~Year_grp)+theme_bw()
map_rich_zoo_grp10%>%filter(!is.na(Year_grp))%>%ggplot()+geom_tile(aes(x,y,fill=n_grp))+
  scale_fill_manual(values = col_code,name='number of species')+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~Year_grp)+theme_bw()


#  ------- write hisotrical distribution map related objects to a  --------
save(list = c('map_rich_benthos_grp10','map_rich_fish_grp10','map_rich_zoo_grp10',
              'col_code','sa_shape','eez_sAfr_sf','sa_shape'),file = 'analysis/historical_map_related.RData')

