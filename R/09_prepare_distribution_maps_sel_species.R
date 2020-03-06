# ---- To be used to prepare distribution maps for selected species from the three groups---------------
# ------- load required libraries -----------------------------------------
libs <-c('tidyverse','sf','raster','SSDM','stringr','rgeos');lapply(libs, library,character.only=T)


#   ----- sets of species  to be considered for distribution maps ---------
sel_fish <-c("Sardinops sagax", "Engraulis encrasicolus", "Etrumeus whiteheadi", 
             "Merluccius capensis", "Merluccius paradoxus", "Trachurus capensis", 
             "Lampanyctodes hectoris", "Austroglossus pectoralis", "Thyrsites atun", 
             "Genypterus capensis", "Scomber japonicus", "Loligo vulgaris reynaudi", 
             "Callorhinchus capensis", "Chelidonichthys queketti", "Maurolicus muelleri", 
             "Chelidonichthys capensis")
sel_inv2 <-c("Chaceon chuni", "Echinus gilchristi", "Exodromidia spinosa", 
             "Marthasterias glacialis", "Mursia cristiata", "Toraster tuberculatus","Emarginula natalensis",
             "Nassarius eusulcatus","Goneplax rhomboides","Luidia atlantidea","Xenophora solarioides","Vexillum sculptile")

sel_zoo <- c("Calanoides carinatus", "Calanoides", "Calanus agulhensis", 
             "Centropages", "Clausocalanus", "Ctenocalanus", "Euphausia lucens", 
             "Funchalia woodwardi", "Funchalia", "Metridia", "Oithona", "Paracalanus pygmaeus", 
             "Pleuromamma", "Rhincalanus", "Sergestes", "Sergia")


#  ---- check list of sdms to be considered -------------------------------
# --- paths  for the sdm object -------------------------------------------
path_benthos_all='analysis/ensemble_out_sa/all_species_ens/ens_benthos_all/'
path_fish_all='analysis/ensemble_out_sa/all_species_ens/ens_fish_all/'
path_zoo_all='analysis/ensemble_out_sa/all_species_ens/ens_zoo_all/'

# -- sdm objects: benthos -------------------------------------------------
stk_bnt_al<-list.files(path_benthos_all)[str_detect(list.files(path_benthos_all),'.RData')]


# --------- sdm objects: fishes -------------------------------------------
stk_fish_al<-list.files(path_fish_all)[str_detect(list.files(path_fish_all),'.RData')]

# ---- sdm objects: zooplankton -------------------------------------------
stk_zoo_al<-list.files(path_zoo_all)[str_detect(list.files(path_zoo_all),'.RData')]


# ------ Process sdm object names and load sdm objects ---------------------------------------------
sel_bnt_al <- stk_bnt_al[str_replace_all(stk_bnt_al,'.RData','')%in%sel_inv2]

sel_fish_al <- stk_fish_al[str_replace_all(stk_fish_al,'.RData','')%in%sel_fish]
sel_zoo <- str_replace_all(stk_zoo_al,'.RData','')


# ----- load sdm objects --------------------------------------------------
stk_bnt_lst_al <-sapply(paste(path_benthos_all,sel_bnt_al,sep = ''), 
                        function(x) mget(load(x)),simplify = TRUE)

stk_fish_lst_al <-sapply(paste(path_fish_all,sel_fish_al,sep = ''), 
                         function(x) mget(load(x)),simplify = TRUE)

stk_zoo_lst_al <-sapply(paste(path_zoo_all,stk_zoo_al,sep = ''), 
                        function(x) mget(load(x)),simplify = TRUE)


#  ---- rename sdm objects ------------------------------------------------
names(stk_bnt_lst_al) = str_replace_all(sel_bnt_al,'.RData','')
names(stk_fish_lst_al) =str_replace_all(sel_fish_al,'.RData','')
names(stk_zoo_lst_al) = str_replace_all(stk_zoo_al,'.RData','')


# -------- prepare base maps and eez shape  file --------------------------
world <- getMap(resolution = "high")
### convert the spatialPolygonsDataFrame to an sf object
world_sf <- st_as_sf(world)
sa_shape <- world_sf%>%dplyr::filter(NAME=='South Africa')%>%lwgeom::st_make_valid()%>%
  st_crop(xmin=15,xmax=35,ymin=-37,ymax=-20)%>%
  group_by(NAME)%>%summarise()

# --- eez shape file ------------------------------------------------------
shp_eezWAfr = readOGR(dsn = 'spatial/shape/',layer = 'eez_200NM_shape_wAfr')
eez_wAfr_sf <- st_as_sf(shp_eezWAfr)
eez_sAfr_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa')%>%group_by(Territory1)%>%summarise()
eez_sAfr_sp <- as_Spatial(eez_sAfr_sf)

eez_bclme_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa'|Territory1=='Democratic Republic of the Congo'|
                                       Territory1=='Namibia'|Territory1=='Angola')%>%mutate(LME='bclme')%>%
  lwgeom::st_make_valid()%>%
  st_crop(xmin=5,xmax=20,ymin=-40,ymax=-3)%>%group_by(LME)%>%summarise()

# ------ Utility functions to be used to process sdm objects --------------

# ----------- get predicted distribution -----------------------------
get_stack_dist <- function(dist_list){
  n_layers <-length(dist_list)
  tbl_store <-tibble(species=names(dist_list),dist_vals=vector('list',length(dist_list)))
  
  for(i in 1:length(dist_list)){
    t_dat <-as_tibble(mask(dist_list[[i]]@projection,eez_sAfr_sp)%>%rasterToPoints())
    tbl_store$dist_vals[[i]]=t_dat
  }
  tbl_store
}


# ----------- get the occurrence/encounter  records -----------------------
get_presence_data <-function(dist_list){
  n_layers <-length(dist_list)
  tbl_store <-tibble(species=names(dist_list),dist_vals=vector('list',length(dist_list)))
  
  for(i in 1:length(dist_list)){
    t_dat <-as_tibble(dist_list[[i]]@data%>%dplyr::select(X,Y,Presence)%>%filter(Presence==1))
    if(nrow(t_dat)<=300){t_dat}else{set.seed(1234);t_dat <- t_dat[sample(1:nrow(t_dat),300,replace = FALSE),]}
    tbl_store$dist_vals[[i]]=t_dat
  }
  tbl_store
}


#  ------- get variable importance ----------------------------------------
get_var_importance <-function(stk_list){
  n_layers = length(stk_list)
  store_tbl <-tibble(species=names(stk_list), var_imp = vector(mode = 'list',length = length(stk_list)))
  for(i in 1:n_layers){
    imp_vars <- stk_list[[i]]@variable.importance%>%t()%>%as.data.frame()%>%
      rownames_to_column(var = 'var_name')%>%arrange(desc(Axes.evaluation))%>%
      mutate(rank= order(Axes.evaluation))
    store_tbl$var_imp[[i]] <- imp_vars
  }
  
  store_tbl
}


# -------------------- predicted distribution: ensemble --------------------

fish_dist_sel <- get_stack_dist(stk_fish_lst_al)%>%unnest(dist_vals)
benthos_dist_sel <- get_stack_dist(stk_bnt_lst_al)%>%unnest(dist_vals)
zoo_dist_sel <- get_stack_dist(stk_zoo_lst_al)%>%unnest(dist_vals)

# ------- extract ocurrence/presence data ---------------------------------
fish_pres_sel <- get_presence_data(stk_fish_lst_al)%>%unnest(dist_vals)
benthos_pres_sel <- get_presence_data(stk_bnt_lst_al)%>%unnest(dist_vals)
zoo_pres_sel <- get_presence_data(stk_zoo_lst_al)%>%unnest(dist_vals)


# ------  write predicted distribution and occcurence data ----------------
save(list = c("benthos_dist_sel", "benthos_pres_sel", "fish_dist_sel", "fish_pres_sel", 
              "zoo_dist_sel", "zoo_pres_sel"),file = 'analysis/distribution/processed_sample_distn_map.RData')


# ---- check plots of distribution map with occurrence locations ----------
fish_dist_sel%>%filter(species%in%sel_fish[1:8])%>%ggplot()+geom_tile(aes(x,y,fill=Probability))+
  scale_fill_gradientn(colors = fields::tim.colors(100),name='probability of\n occurrence')+
  geom_point(data=fish_pres_sel%>%filter(species%in%sel_fish[1:8]),aes(X,Y),shape=1,size=0.3)+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~species)+theme_bw()

fish_dist_sel%>%filter(species%in%sel_fish[9:16])%>%ggplot()+geom_tile(aes(x,y,fill=Probability))+
  scale_fill_gradientn(colors = fields::tim.colors(100),name='probability of\n occurrence')+
  geom_point(data=fish_pres_sel%>%filter(species%in%sel_fish[9:16]),aes(X,Y),shape=1,size=0.3)+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~species)+theme_bw()


benthos_dist_sel%>%filter(species%in%sel_inv2[1:6])%>%ggplot()+geom_tile(aes(x,y,fill=Probability))+
  scale_fill_gradientn(colors = fields::tim.colors(100),name='probability of\n occurrence')+
  geom_point(data=benthos_pres_sel%>%filter(species%in%sel_inv2[1:6]),aes(X,Y),shape=1,size=0.3)+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~species)+theme_bw()

benthos_dist_sel%>%filter(species%in%sel_inv2[7:12])%>%ggplot()+geom_tile(aes(x,y,fill=Probability))+
  scale_fill_gradientn(colors = fields::tim.colors(100),name='probability of\n occurrence')+
  geom_point(data=benthos_pres_sel%>%filter(species%in%sel_inv2[7:12]),aes(X,Y),shape=1,size=0.3)+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~species)+theme_bw()

zoo_dist_sel%>%filter(species%in%sel_zoo[1:8])%>%ggplot()+geom_tile(aes(x,y,fill=Probability))+
  scale_fill_gradientn(colors = fields::tim.colors(100),name='probability of\n occurrence')+
  geom_point(data=zoo_pres_sel%>%filter(species%in%sel_zoo[1:8]),aes(X,Y),shape=1,size=0.3)+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~species)+theme_bw()

zoo_dist_sel%>%filter(species%in%sel_zoo[9:16])%>%ggplot()+geom_tile(aes(x,y,fill=Probability))+
  scale_fill_gradientn(colors = fields::tim.colors(100),name='probability of\n occurrence')+
  geom_point(data=zoo_pres_sel%>%filter(species%in%sel_zoo[9:16]),aes(X,Y),shape=1,size=0.3)+
  geom_sf(data=sa_shape,fill=NA)+geom_sf(data=eez_sAfr_sf,fill=NA)+labs(y='latitude',x='longitude')+
  facet_wrap(~species)+theme_bw()


#  -------- extract index of relative importance --------------------------

#  ---- get sdm objects ---------------------------------------------------
stk_bnt_lst_vImp <-sapply(paste(path_benthos_all,stk_bnt_al,sep = ''), 
                          function(x) mget(load(x)),simplify = TRUE)

stk_fish_lst_vImp <-sapply(paste(path_fish_all,stk_fish_al,sep = ''), 
                           function(x) mget(load(x)),simplify = TRUE)

stk_zoo_lst_vImp <-sapply(paste(path_zoo_all,stk_zoo_al,sep = ''), 
                          function(x) mget(load(x)),simplify = TRUE)


# ---- rename sdm objects -------------------------------------------------
names(stk_bnt_lst_vImp) = str_replace_all(stk_bnt_al,'.RData','')
names(stk_fish_lst_vImp) =str_replace_all(stk_fish_al,'.RData','')
names(stk_zoo_lst_vImp) = str_replace_all(stk_zoo_al,'.RData','')


#  ---- extract the index of importance -----------------------------------
vImp_benthos <- get_var_importance(stk_bnt_lst_vImp)
vImp_fish <- get_var_importance(stk_fish_lst_vImp)
vImp_zoo <- get_var_importance(stk_zoo_lst_vImp)


# ----- summarise index of importance -------------------------------------
sum_vImp_fish <- vImp_fish%>%unnest(var_imp)%>%group_by(var_name)%>%summarise(mean_rank = mean(rank,na.rm=T))
sum_vImp_benthos <- vImp_benthos%>%unnest(var_imp)%>%group_by(var_name)%>%summarise(mean_rank = mean(rank,na.rm=T))
sum_vImp_zoo <- vImp_zoo%>%unnest(var_imp)%>%group_by(var_name)%>%summarise(mean_rank = mean(rank,na.rm=T))


all_var_imp  = rbind(vImp_fish%>%unnest(var_imp)%>%mutate(type='Fish'),
                     vImp_benthos%>%unnest(var_imp)%>%mutate(type='Benthos'),
                     vImp_zoo%>%unnest(var_imp)%>%mutate(type='Zooplankton'))


# --------- plot the Index of importance ----------------------------------
all_var_imp%>%ggplot()+geom_boxplot(aes(x = var_name,rank))+geom_violin(aes(var_name,rank))+
  labs(y="index of importance",x='variables')+
  facet_wrap(~type,nrow = 3)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))


all_var_imp%>%ggplot()+geom_boxplot(aes(var_name,rank),outlier.colour = NA)+
  labs(y="index of importance",x='variables')+
  facet_wrap(~type,nrow = 3)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))
all_var_imp%>%ggplot()+geom_violin(aes(var_name,Axes.evaluation),draw_quantiles = c(0.25,0.5,0.75))+
  labs(y="index of importance",x='variables')+
  facet_wrap(~type,nrow = 3)+theme_bw()+
  theme(axis.text.x = element_text(angle = 90))


# ---- write to file the index of importance ------------------------------
write_csv(all_var_imp,'analysis/distribution/variable_importance_all_group.csv')

