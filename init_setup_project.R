################## Prepare and setup project ############################################
libs_ini <- c('purrr','here');lapply(libs_ini,library, character.only=T)

dir_list = c('R',here('data/raw/bio'),here('data/raw/env'),here('data/tidy/bio'),here('data/tidy/env/sa'),
             here('data/sel_env/sa'),
             'report',here('data/tidy/env/bclme'),
             'figures',here('analysis/ensemble_out_sa/all_species_ens'),
             here('analysis/biodiversity'),here('analysis/distribution'),
             here('spatial/shape'),'temp_raster',here::here('bathy_derived_indices/sa'),
             here::here('bathy_derived_indices/bclme'))
cust_dir_create <- function(x){ifelse(dir.exists(x),dir.exists(x),dir.create(x,recursive = TRUE))}

map(dir_list,cust_dir_create)
