# ----  link the aphia based trait datt to the  occurrence data ------------

# ------ load  libraries to be used  --------------------------------------
libs <- c('tidyverse', 'sf','rgdal', 'raster');lapply(libs, library, character.only=TRUE)

#  -------- load the attributes data from  worrrms ------------------------
load('data/tidy/bio/all_attribute_table.RData')

attr_tibble_store_sel <-attr_tibble_store%>%
  mutate(sel_rslt = purrr::map(rslt,function(x) x%>%ungroup()%>%dplyr::select(AphiaID,adult)))%>%unnest(cols = sel_rslt)%>%
  dplyr::rename(aphiaID = AphiaID)%>%dplyr::select(-c(index,rslt))%>%filter(!is.na(aphiaID))


#  ------- load the occurrence data ---------------------------------------
db_conn <- dbConnect(RSQLite::SQLite(),'data/raw/bio/allDat_db.sqlite')
dbListTables(db_conn)
dbListFields(db_conn,'allData')
al_fields <-dbListFields(db_conn,'allData')

############ database query using dplyr grammer ########################
tbl_conn_ref = tbl(db_conn,'allData')

system.time(selDat_allData <- tbl_conn_ref %>% dplyr::select(decimalLongitude, decimalLatitude,
                                                             phylum, class,order, family, genus, scientificName, species,
                                                             eventDate,institutionID,ownerInstitutionCode,aphiaID)%>%as_tibble())#%>%filter(phylum=='Chordata')

#  ---- rename coordinates ------------------------------------------------
names(selDat_allData)[c(1:2)]=c("long","lat")

#   ---- process collection date ------------------------------------------
selDat_allData <-selDat_allData%>%mutate(collection_date = as_datetime(eventDate),
                                         year=year(collection_date),month=month(collection_date))

aphia_species <- selDat_allData%>%filter(!is.na(aphiaID)&!is.na(scientificName))%>%
  dplyr::select(aphiaID,species,scientificName)%>%distinct()

# --- join the attribute table containing broad functional groups ---------
aphia_spp_attr <-aphia_species%>%left_join(attr_tibble_store_sel)


# --- join species attribute data to the whole occurrence data ------------
selDat_allData<-selDat_allData%>%left_join(aphia_spp_attr%>%dplyr::select(species,scientificName,adult))


# ----------------------- subset to EEZ ----------------
#  ----- get eez shape files and subset the data --------------------------
eez_wAfr <- rgdal::readOGR(dsn = "spatial/shape/", 
                           layer = "eez_200NM_shape_wAfr", stringsAsFactors = TRUE)


# --- conver to an sf object for processing -------------------------------
eez_wAfr_sf <- st_as_sf(eez_wAfr)
eez_sAfr_sf <- eez_wAfr_sf%>%filter(Territory1=='South Africa')
eez_sAfr_sf%>%ggplot()+geom_sf() # check the eez for south africa
eez_bclmeAgl_sf <- eez_wAfr_sf%>%
  filter(Territory1=='South Africa'|Territory1=='Namibia'|Territory1=='Angola'|
           Territory1=='Democratic Republic of the Congo')
eez_bclme_sf<-eez_bclmeAgl_sf%>%st_crop(xmin = 8.2,ymin = -38.175,xmax = 20,ymax = -5.029)%>%mutate(lme='BCLME')
eez_bclme_sf%>%ggplot()+geom_sf()+theme_bw() # check eez for the bclme
eez_bclme_sfJ <-eez_bclme_sf%>%group_by(lme)%>%summarise()

#   ---- subset the occurrence data ----------------------------------------
crd_dat_al = selDat_allData[,c(1:2)]
system.time(crd_sf <- st_as_sf(crd_dat_al,coords = c('long','lat'),crs=4326))
system.time(occInd_bn<-st_intersects(eez_bclme_sfJ,crd_sf)) #occuring within BCLME EEZ
system.time(occInd_sa<-st_intersects(eez_sAfr_sf,crd_sf)) # occurring with SA EEZ


# subset  coordinates of occurrence  within EEZ ---------------------------
selDat_allDataF_bn <- selDat_allData[occInd_bn[[1]],]%>%filter(!is.na(genus))
selDat_allDataF_sa <- selDat_allData[occInd_sa[[1]],]%>%filter(!is.na(genus))

#  --- Subset the data by region: bclme/sa and by functional groups --------------------------------

sel_classes = c('Cephalopoda','Elasmobranchii','Actinopterygii','Holocephali')


# -------- BCLME ----------------------------------------------------------
sel_bclme_fish <- selDat_allDataF_bn%>%filter(adult%in%'fish'|class%in%'Cephalopoda')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()
sel_bclme_benthos <- selDat_allDataF_bn%>%filter((adult%in%'benthos'|adult%in%'macrobenthos'|
                                                    adult%in%'epibenthos')&!class%in%'Cephalopoda')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()
sel_bclme_phytoZooplankton <- selDat_allDataF_bn%>%filter(adult%in%'phytoplankton'|adult%in%'zooplankton')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()
sel_bclme_birds <- selDat_allDataF_bn%>%filter(adult%in%'birds')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()
sel_bclme_mamals <- selDat_allDataF_bn%>%filter(adult%in%'mammals')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()


# --------- South Africa --------------------------------------------------
sel_beng_fish <- selDat_allDataF_sa%>%filter(adult%in%'fish'|class%in%'Cephalopoda')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()
sel_beng_benthos <- selDat_allDataF_sa%>%filter((adult%in%'benthos'|adult%in%'macrobenthos'|
                                                   adult%in%'epibenthos')&!class%in%'Cephalopoda')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()
sel_beng_phytoZooplankton <- selDat_allDataF_sa%>%filter(adult%in%'phytoplankton'|adult%in%'zooplankton')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()
sel_beng_birds <- selDat_allDataF_sa%>%filter(adult%in%'birds')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()
sel_beng_mamals <- selDat_allDataF_sa%>%filter(adult%in%'mammals')%>%
  dplyr::select(long,lat,year,month,class,family,genus,species,scientificName,adult)%>%distinct()


#  ----- write processed occurrence data to a file ------------------------

# ---- SA --------------------------------------------------------------
write_csv(sel_beng_benthos,'data/raw/bio/benthos_SA.csv')
write_csv(sel_beng_birds,'data/raw/bio/birds_SA.csv')
write_csv(sel_beng_fish,'data/raw/bio/fish_SA.csv')
write_csv(sel_beng_mamals,'data/raw/bio/mammals_SA.csv')
write_csv(sel_beng_phytoZooplankton,'data/raw/bio/phytoZooplankton_SA.csv')


# ----- BCLME ----------------------------------------------------------------
write_csv(sel_bclme_benthos,'data/raw/bio/benthos_bclme.csv')
write_csv(sel_bclme_birds,'data/raw/bio/birds_bclme.csv')
write_csv(sel_bclme_fish,'data/raw/bio/fish_bclme.csv')
write_csv(sel_bclme_mamals,'data/raw/bio/mammals_bclme.csv')
write_csv(sel_bclme_phytoZooplankton,'data/raw/bio/phytoZooplankton_bclme.csv')



