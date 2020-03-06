# ------  download, processs , and store processed OBIS data --------------

# --------- load all the required libraries  ------------------------------
libs <- c('tidyverse','robis',"DBI");lapply(libs,library,character.only=TRUE)


# -----  download the  occurrrence/encounter data for a selected region --------
#  ----   boundary of the spatial domain given as WKT  --------------------
# ------ # Note- the code for downloding data takes long to run and results in big file -----------
tSearch_obis <- occurrence(geometry = "POLYGON((-20 -40,-20 40,60 40,60 -40,-20 -40))")


#  --- first write the raw occurrence data to a file ------------------------
write_csv(tSearch_obis, path = 'data/raw/bio/raw-OBIS_data.csv')


# ---- process downloaded occurrence data ---------------------------------
# --- raw occurrence data downloaded above--------------------------------------
system.time(raw_occ <- read_csv('data/raw/bio/raw-OBIS_data.csv')) 

#  --- write the entire data into an sqlite database for ease of accesss and processing ------
Create_db <- dbConnect(RSQLite::SQLite(),"data/raw/bio/allDat_db.sqlite")
dbWriteTable(Create_db,"allData",raw_occ)
dbDisconnect(Create_db)

#   read only part of the database for processing  ----------------------

#  ---- make connection to the database and explore its content -----------
theConn = dbConnect(RSQLite::SQLite(),'data/raw/bio/allDat_db.sqlite')
dbListTables(theConn)
dbListFields(theConn,'allData')


# ----- extract only desired fieldfs from the database --------------------
theConn_ref = tbl(theConn,'allData')
selDat_allData <- theConn_ref %>% dplyr::select(decimalLongitude, decimalLatitude,depth,
                                                maximumDepthInMeters,minimumDepthInMeters,
                                                phylum, order, family, genus, scientificName, species,
                                                habitat)%>%as.data.frame()#%>%filter(phylum=='Chordata')
names(selDat_allData)[c(1:5)]=c("long","lat",'depth','max_depth',"min_depth")


#  ------ save data with selectd fields to a file  ------------------------
write_csv(selDat_allData,path = 'data/tidy/bio/selected_fields_obis.csv')
dbDisconnect(theConn)
dbDisconnect(theConn_ref)

# -----  Or the data can futher be processed as follows: first connect to the database ------------------


# ---- As the occurrence/encounter data is for any taxonomic group in a region (it has to be substted
#  to the EEZ of for the selected region --------
#  ----- get the shape file for the EEZ of the region ---------------------


#  --- subset the occcurrence/enccounter data to the EEZ ------------------


#  --------- write the occurrence/encounter data to a file ----------------







