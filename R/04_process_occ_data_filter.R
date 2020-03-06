# ------ further processing of occurrence data: filtering -----------------
#  ---- load required libraries -------------------------------------------
libs <- c('tidyverse','SSDM');lapply(libs, library,character.only=TRUE) 

# load occurrence data ----------------------------------------------------
occ_benthos <- read_csv('data/raw/bio/benthos_SA.csv')
occ_fish <- read_csv('data/raw/bio/fish_SA.csv')
occ_phytoZoo <- read_csv('data/raw/bio/phytoZooplankton_SA.csv')


# ----- process occurencee data: benthos (total numbers of record)----------------------------------
benthos_occ_yr <- occ_benthos%>%dplyr::select(year,scientificName)%>%group_by(year,scientificName)%>%
  summarise(n=n())%>%ungroup()
benthos_occ_sel8695 <-occ_benthos%>%filter(year>=1986&year<=1995)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
benthos_occ_sel9500 <-occ_benthos%>%filter(year>1995)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
benthos_occ_selAl <-occ_benthos%>%filter(year>=1986)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()


# ----- process occurencee data: phytoplankton/zooplankton (total numbers of record)-------------------------
phytozoo_occ_yr <- occ_phytoZoo%>%dplyr::select(year,scientificName,adult)%>%group_by(year,scientificName)%>%
  summarise(n=n())%>%ungroup()
phyto_occ_sel8695 <-occ_phytoZoo%>%filter(adult=='phytoplankton'&year>=1986&year<=1995)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
zoo_occ_sel8695 <-occ_phytoZoo%>%filter(adult=='zooplankton'&year>=1986&year<=1995)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
phyto_occ_sel9500 <-occ_phytoZoo%>%filter(adult=='phytoplankton'&year>1995)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
zoo_occ_sel9500 <-occ_phytoZoo%>%filter(adult=='zooplankton'&year>1995)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
zoo_occ_selAl <-occ_phytoZoo%>%filter(adult=='zooplankton'&year>=1986)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
phyto_occ_selAl <-occ_phytoZoo%>%filter(adult=='phytoplankton'&year>=1986)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()


# ----- process occurencee data: fishes (total numbers of records)------------------------------
fish_occ_yr <- occ_fish%>%dplyr::select(year,scientificName)%>%group_by(year,scientificName)%>%
  summarise(n=n())%>%ungroup()
fish_occ_sel8695 <-occ_fish%>%filter(year>=1986&year<=1995)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
fish_occ_sel9500 <-occ_fish%>%filter(year>1995)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()
fish_occ_selAl <-occ_fish%>%filter(year>=1986)%>%
  dplyr::select(long,lat,scientificName)%>%distinct()%>%group_by(scientificName)%>%count()%>%ungroup()

# ----  subset species that occur atleast 20 times (this falls later to a minimum of about 10 records ------
# -- for benthos by period ------------------------------------------------
sel_benthos_al8690 <- benthos_occ_sel8695%>%filter(n>=20)
sel_benthos_al9500 <- benthos_occ_sel9500%>%filter(n>=20)
sel_benthos_al <- benthos_occ_selAl%>%filter(n>=20)

#  ---- for zooplankton by period -----------------------------------------
sel_zoo_al8690 <- zoo_occ_sel8695%>%filter(n>=20)
sel_zoo_al9500 <- zoo_occ_sel9500%>%filter(n>=20)
sel_zoo_al <- zoo_occ_selAl%>%filter(n>=20)


# ---- for fishes by period -----------------------------------------------
sel_fish_al8690 <- fish_occ_sel8695%>%filter(n>=20) 
sel_fish_al9500 <- fish_occ_sel9500%>%filter(n>=20)
sel_fish_al <- fish_occ_selAl%>%filter(n>=20)


# ----  Subset occurrence data for the three groups: based on above computed occurrence threshold -----------
#   subsetted: fish -------------------------------------------------------
occ_fish_8690 <- occ_fish%>%filter(!is.na(year)&year>=1986&year<=1995)%>%
  filter(scientificName%in%sel_fish_al8690$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()
occ_fish_9500 <-occ_fish%>%filter(!is.na(year)&year>1995)%>%
  filter(scientificName%in%sel_fish_al9500$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()
occ_fish_Al <- occ_fish%>%filter(!is.na(year)&year>=1986)%>%
  filter(scientificName%in%sel_fish_al$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()

occ_fish_Al_taxInfo <- occ_fish%>%filter(!is.na(year)&year>=1986)%>%
  filter(scientificName%in%sel_fish_al$scientificName)%>%
  dplyr::select(class,family,genus,species,scientificName)%>%distinct()

# subsetted:  benthos -----------------------------------------------------
occ_benthos_8690 <- occ_benthos%>%filter(!is.na(year)&year>=1986&year<=1995)%>%
  filter(scientificName%in%sel_benthos_al8690$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()
occ_benthos_9500 <-occ_benthos%>%filter(!is.na(year)&year>1995)%>%
  filter(scientificName%in%sel_benthos_al9500$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()
occ_benthos_Al <- occ_benthos%>%filter(!is.na(year)&year>=1986)%>%
  filter(scientificName%in%sel_benthos_al$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()


# ----  subsetted:  zooplankton -------------------------------------------
occ_zoo_8690 <- occ_phytoZoo%>%filter(adult=='zooplankton'&!is.na(year)&year>=1986&year<=1995)%>%
  filter(scientificName%in%sel_zoo_al8690$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()
occ_zoo_9500 <-occ_phytoZoo%>%filter(adult=='zooplankton'&!is.na(year)&year>1995)%>%
  filter(scientificName%in%sel_zoo_al9500$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()
occ_zoo_Al <- occ_phytoZoo%>%filter(adult=='zooplankton'&!is.na(year)&year>=1986)%>%
  filter(scientificName%in%sel_zoo_al$scientificName)%>%dplyr::select(long,lat,scientificName)%>%distinct()


#   ----    write processed data to a file --------------------------------
write_csv(x = occ_benthos_8690,path = 'data/tidy/bio/processed_occ_data_benthos_8690.csv')
write_csv(x = occ_benthos_9500,path = 'data/tidy/bio/processed_occ_data_benthos_9500.csv')
write_csv(x = occ_benthos_Al,path = 'data/tidy/bio/processed_occ_data_benthos_all.csv')
write_csv(x = occ_fish_8690,path = 'data/tidy/bio/processed_occ_data_fish_8690.csv')
write_csv(x = occ_fish_9500,path = 'data/tidy/bio/processed_occ_data_fish_9500.csv')
write_csv(x = occ_fish_Al,path = 'data/tidy/bio/processed_occ_data_fish_all.csv')
write_csv(x = occ_fish_Al_taxInfo,path = 'data/tidy/bio/processed_occ_data_fish_all_withTaxonInfo.csv')

write_csv(x = occ_zoo_8690,path = 'data/tidy/bio/processed_occ_data_zooplankton_8690.csv')
write_csv(x = occ_zoo_9500,path = 'data/tidy/bio/processed_occ_data_zooplankton_9500.csv')
write_csv(x = occ_zoo_Al,path = 'data/tidy/bio/processed_occ_data_zooplankton_all.csv')


