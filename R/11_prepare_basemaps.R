# ------ to be used to prepare figure 2 (basemaps) ------------------------
libLst<-c('tidyverse','raster','sf','rworldmap','rgdal','purrr','ggsn','grid','ggrepel',"rworldxtra")
lapply(libLst,library,character.only=T)


# --- location  coordinates ------------------------------------------------
loc_cords <-tibble(name=c('St Helena Bay', 'Saldanha','False bay',"L` Agulhas",'Mossel Bay', 'Plattenberg Bay',
                          'Cape St Franscis','Port Elizabeth','Durban'),
                   lon=c(18.008457,17.924521,18.648847, 20.026122, 22.143147, 23.388734, 24.850262, 
                         25.611973, 31.011833),
                   lat=c(-32.744267,-33.038005,-34.192827,-34.823214,-34.191164,-34.105154,
                         -34.191058,-33.945024,-29.919980))
incl_name = c('St Helena Bay', 'Mossel Bay', 'Plattenberg Bay','Port Elizabeth','Durban')

loc_cords_sel <- loc_cords%>%dplyr::filter(name%in%incl_name)

# --- load raster data to be used  when preparing basemaps -------------------------
system.time(load('data/tidy/env/sa/bottom_raster_sa.RData'))

d_sa <-index_mask_sa_phys_bot[[8]]

d_sa_cont <- rasterToContour(d_sa,nlevels=3,levels=c(200,500,1000))%>%
  raster::crop(.,extent(c(10,35,-37,-25)))%>%st_as_sf()%>%
  mutate(n_level =  reorder(factor(as.character(level)),X = as.numeric(as.character(level))))

# ----- get eez shape files and prepare country shape files --------------
shp_eezWAfr = readOGR(dsn = 'spatial/shape',layer = 'eez_200NM_shape_wAfr')
eez_wAfr_sf <- st_as_sf(shp_eezWAfr)
eez_sAfr_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa')%>%group_by(Territory1)%>%summarise()
eez_sAfr_sp <- as_Spatial(eez_sAfr_sf)

# -----  get global  country shape file  and process ----------------------
world <- getMap(resolution = "high")
### convert the spatialPolygonsDataFrame to an sf object
world_sf <- st_as_sf(world)
al_afr_shape <- world_sf%>%dplyr::filter(REGION=='Africa')%>%lwgeom::st_make_valid()%>%
  st_crop(xmin=-21,xmax=59,ymin=-36,ymax=37)%>%
  group_by(NAME)%>%summarise()
sa_shape_al <- world_sf%>%dplyr::filter(NAME=='South Africa')%>%st_crop(xmin=15,xmax=35,ymin=-37,ymax=-22)%>%
  group_by(NAME)%>%summarise()

sa_shape_sp <- as_Spatial(sa_shape_al)

bathy_base<- ggplot()+ 
  geom_sf(data = sa_shape_al,  color = "black", fill = "grey80") +
  geom_sf(data = d_sa_cont,aes(linetype=n_level),show.legend = 'line')+
  geom_sf(data = eez_sAfr_sf, fill=NA)+
  scale_linetype_manual(name='depth',values = c("dotted", "solid", "dashed"))+
  geom_text_repel(data = loc_cords_sel,aes(x = lon,y = ,lat,label=name),nudge_y = 0.9)+
  xlab('longitude')+ylab('latitude')+ 
  theme_bw()

bathy_base_scale <- bathy_base+
  ggsn::north(data = sa_shape_al,location = 'topright',symbol = 1,scale = 0.2)+#,x = 0.6,y = 0.8,scale = .15)
  ggsn::scalebar(eez_sAfr_sf, dist = 100, dist_unit = "km", st.size = 2, 
                 transform = TRUE, model = "WGS84",location = 'bottomright',border.size = 1)

fil_col_opts<-c("gray44", "gray81", "khaki","floralwhite", "ghostwhite","gray81", "gray46", "gray33")


inset_map<-ggplot()+geom_sf(data=al_afr_shape,colour=fil_col_opts[1],fill=fil_col_opts[5])+
  geom_sf(data = sa_shape_al,fill='gray46')+
  theme_bw()+labs(x=NULL,y=NULL)+
  theme(axis.text.x =element_blank(),axis.text.y= element_blank(), axis.ticks=element_blank(),axis.title.x =element_blank(),
        axis.title.y= element_blank())

#  ---  put together all figures  -----------------------------------------
tiff(filename = 'figures/map_study_reg_with_inset_sa.tif',width = 20,height = 18,units = 'cm',compression = 'lzw',res = 320)
pdf(file = 'figures/map_study_reg_with_inset_sa.pdf',width = 20,height = 18)
grid.newpage()
v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.28, height = 0.17, x = 0.2, y = 0.8) #plot area for the inset map
print(bathy_base_scale,vp=v1) 
print(inset_map,vp=v2)
dev.off()

#  ----- write prepared basemap related objects to a file -----------------
save(list = c('al_afr_shape','sa_shape_al','d_sa_cont','eez_sAfr_sf',
              'loc_cords_sel','sa_shape_sp', 'eez_sAfr_sp'),file = 'analysis/base_map_related.RData')



