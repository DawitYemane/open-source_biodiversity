# --------- To be used to download, process, and write to a file environmental data --------

# ---- load libraries to be used to download and process BioORACLE  envrionmental layers -----
libs <- c('sdmpredictors', 'tidyverse','raster','sf');lapply(libs,library, character.only=TRUE)

# ----- download  selected  sets of environmental layers from BioORACLE ------------------

#  ----  explore list of layers  available from  BioORACLE ----------------
listOfLayers<-list_layers()

# --- subset list of of layers to the bottom and and surface --------------
lLayer_b = listOfLayers$layer_code[str_detect(string = listOfLayers$layer_code,pattern = "_bdmean")]
lLayer_b2 = lLayer_b[!str_detect(lLayer_b,'range|curve|tmax|tmin')]
lLayer_s = listOfLayers$layer_code[str_detect(string = listOfLayers$layer_code,pattern = "_ss")]
lLayer_s2 = lLayer_s[str_detect(lLayer_s,'BO2_')]%>%.[!str_detect(., 'range|curve|tmax|tmin|ice')]

#  --- get data from BioORACLE  repository and store in the destination directory ----------------------

#  ---- Note: if the files already exsist in  the directory the function just loads the raser files --------
albot_layer <- load_layers(layercodes = lLayer_b2 , equalarea=FALSE, rasterstack=TRUE,
                            datadir = 'data/raw/env/bottom/') 
alsurf_layer <- load_layers(layercodes = lLayer_s2 , equalarea=FALSE, rasterstack=TRUE,
                           datadir = 'data/raw/env/surface/') 

# ------ proccess downloaded environmental layers and subset to the EEZ region -------------------------



# ------ write processed  environmental layers to a file ------------------







