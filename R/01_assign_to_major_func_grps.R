# ------ Group taxons into three  major functional groups  ----------------

#  ----- load libraries and R source files to be used  --------------------
libs <- c('tidyverse','worrms',"DBI",'parallel');lapply(libs, library, character.only=TRUE)


# --- load processed OBIS occcurrence data  -------------------------------
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



#   ----  check the APHIA ID of each taxon --------------------------------------------
unq_aphiaID <- selDat_allData%>%dplyr::select(aphiaID)%>%distinct()
unq_aphiaID <- unq_aphiaID%>%mutate(index = round(seq(1,100,length.out = nrow(unq_aphiaID)))) # for batch processing

#   ---- subset the data sets with non-missing aphiaID, species na --------
aphia_species <- selDat_allData%>%filter(!is.na(aphiaID)&!is.na(scientificName))%>%
  dplyr::select(aphiaID,species,scientificName)%>%distinct()

#  ---  extraction of aphia ID can be done sequentially or in parallel ------
# --- NOTE: this takes long -----------------------------------------------
system.time(att_tble_spp <- unq_aphiaID%>%
              group_by(aphiaID) %>%
              do(get_worms_fgrp(AphiaID = .$aphiaID)))


# ---- better version  of the above approach ------------------------------
attr_tibble_store <- tibble(index=1:100,rslt = vector('list',length = 100))

for(i in unique(unq_aphiaID$index)){
  sppID <-tibble(AphiaID=unq_aphiaID$aphiaID[unq_aphiaID$index==i])
  
  sppID%>%
    group_by(AphiaID) %>%
    do(get_worms_fgrp(AphiaID = .$AphiaID))->t_dat
  
  attr_tibble_store$rslt[[i]]=t_dat
  Sys.sleep(1)
  
  cat(i, " ... batch done\n")
  cat((length(unique(unq_aphiaID$index))-i), "... more batch to do\n\n")
}

#  ---  write extracted taxon atrributes (includng  broad  functional groups --------
save(list = 'attr_tibble_store',file = 'data/tidy/bio/all_attribute_table.RData')



