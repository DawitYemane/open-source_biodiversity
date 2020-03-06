
stacked_biodiversity <- function(sdm_list,approach = "standard",method = 'bSSDM',rep.B=1000,Env = env_rast_sel,
                                 verbose=TRUE){
  
  stack_method_list <-c('bSSDM','pSSDM','Bernoulli','MaximumLikelihood','PRR.pSSDM')
  diversity.map=NULL
  enms = NULL
  
  if(approach=='standard'){
    if (method == "bSSDM") {
      # Threshold and sum (Calabrese et al, 2014)
      if (verbose)
        cat("\n Local species richness computed by thresholding and then summing. \n")
      diversity.map <- sum(stack(lapply(sdm_list, function(x)
        reclassify(x@projection,
                   c(-Inf,x@evaluation$threshold,0,x@evaluation$threshold,Inf,1))
      )),na.rm = TRUE)
    }
    
    if (method == "pSSDM") {
      # Individual probabilities sum (Calabrese et al, 2014)
      if (verbose)
        cat("\n Local species richness computed by summing individual probabilities. \n")
      diversity.map <- sum(stack(lapply(sdm_list, function(x) x@projection)),na.rm = TRUE)
    }
    
    if (method == "Bernoulli") {
      # Random Bernoulli distribution (Calabrese et al, 2014)
      if (verbose)
        cat("\n Local species richness computed by drawing repeatedly from a Bernoulli distribution. \n")
      proba <- stack(lapply(sdm_list, function(x) x@projection))
      diversity.map <- calc(proba,fun = function(...) {
        x <- c(...)
        x[is.na(x)] <- 0
        return(rbinom(lengths(x), rep.B, x))
      }, forcefun = FALSE)
      # diversity.map <- sum(diversity.map)/length(sdm_list)/rep.B # original function
      diversity.map <- sum(diversity.map)/rep.B#length(sdm_list)/rep.B
    }
    
    if (method == "MaximumLikelihood") {
      # Maximum likelihood (Calabrese et al, 2014)
      if (verbose)
        cat("\n Local species richness computed by maximum likelihood adjustment. \n")
      diversity.map <- sum(stack(lapply(sdm_list, function(x)
        reclassify(x@projection,
                   c(-Inf,x@evaluation$threshold,0,x@evaluation$threshold,Inf,1))
      )),na.rm = T)
      
      Richness <- .richness(sdm_list)
      SSDM_Richness <- values(mask(diversity.map, Richness))
      SSDM_Richness <- SSDM_Richness[-which(is.na(SSDM_Richness))]
      Richness <- values(Richness)
      Richness <- Richness[-which(is.na(Richness))]
      fit <- lm(Richness ~ SSDM_Richness)
      #fit <- lm(Richness[Richness>0] ~ SSDM_Richness[Richness>0])
      a <- fit$coefficients[1]
      b <- fit$coefficients[2]
      diversity.map <- a + b * diversity.map
    }
    
    
    if (method == "PRR.pSSDM") {
      # Probability ranking with MEM (SESAM, D'Amen et al, 2015)
      if (verbose)
        cat("\n Local species richness computed by probability ranking from pSSDM. \n")
      diversity.map <- sum(stack(lapply(sdm_list, function(x) x@projection)),na.rm = TRUE)
      enms <- .PRR(sdm_list, diversity.map)
    }
    
    if (method == "PRR.MEM") {
      # Probability ranking with MEM (SESAM, D'Amen et al, 2015)
      if (verbose)
        cat("\n Local species richness computed by probability ranking from MEM. \n")
      diversity.map <- .MEM(sdm_list,Env)@projection
      enms <- .PRR(sdm_list, diversity.map)
    }
    
    
    return(list(
      diversity.map = diversity.map,
      enms = enms
    ))
  }
  if(approach=='spatial'){
    
  }
  
}

##### Internals ####

.richness <- function(sdm_list){
  Richness <- reclassify(sdm_list[[1]]@projection, c(-Inf, Inf, 0))
  for (i in seq_len(length(sdm_list)))
    Richness <- Richness + rasterize(
      SpatialPoints(sdm_list[[i]]@data[1:2]),
      Richness, field = sdm_list[[i]]@data$Presence,
      background = 0)
  if (all(values(Richness) %in% c(0, 1, NA)))
    stop("Observed Richness is always equal to 1, modelled richness can't be adjusted !")
  return(Richness)
}

.PRR <- function(sdm_list, Richness){
  # Readjust each enm binary map
  richnesses <- values(Richness)
  names(richnesses) <- seq_len(length(richnesses))
  richnesses <- as.list(richnesses)
  probabilities <- lapply(lapply(sdm_list, FUN = slot, name = "projection"),
                          values)
  probabilities <- lapply(probabilities, function(x) {
    names(x) <- rep(seq_len(length(probabilities[[1]])))
    return(x)
  })
  probabilities <- lapply(probabilities, `[`, names(probabilities[[1]]))
  probabilities <- apply(do.call(rbind, probabilities), 2, as.list)
  binaries <- lapply(lapply(sdm_list, FUN = slot, name = "binary"),
                     values)
  binaries <- lapply(binaries, function(x) {
    names(x) <- rep(seq_len(length(binaries[[1]])))
    return(x)
  })
  binaries <- lapply(binaries, `[`, names(binaries[[1]]))
  binaries <- apply(do.call(rbind, binaries), 2, as.list)
  binaries <- mapply(function(rich, probability, binary) {
    if (!is.na(rich)) {
      ord <- order(unlist(probability), decreasing = TRUE)
      binary <- unlist(binary)
      if (length(ord) <= rich) {
        binary[ord] <- 1
      } else {
        binary[ord[1:rich]] <- 1
        binary[ord[rich + seq_len(length(ord))]] <- 0
      }
      binary <- as.list(binary)
    }
    return(binary)
  }, rich = richnesses, probability = probabilities, binary = binaries,
  SIMPLIFY = FALSE)
  binaries <- lapply(binaries, `[`, names(binaries[[1]]))
  binaries <- apply(do.call(rbind, binaries), 2, as.list)
  binaries <- lapply(binaries, unlist)
  binaries <- lapply(binaries, unname)
  
  mapply(function(enm, binary) {
    values(enm@binary) <- binary
    return(enm)
  }, enm = sdm_list, binary = binaries, SIMPLIFY = FALSE)
}

### check this section of the code
.MEM <- function(obj, Env){
  occ <- data.frame(rasterToPoints(.richness(obj), function(x) x > 0))
  maxOcc <- max(occ$layer) # Reucing occ for algorithms
  occ$layer <- occ$layer/max(maxOcc)
  algo <- unlist(
    strsplit(obj@enms[[1]]@parameters$algorithms,
             ".", fixed = TRUE))[-1]
  if("MAXENT" %in% algo)
    algo <- algo[-which(algo == "MAXENT")]
  MEM <- ensemble_modelling(algorithms = algo,
                            Occurrences = occ, Env = Env, Xcol = "x",
                            Ycol = "y", Pcol = "layer", rep = obj@enms[[1]]@parameters$rep,
                            name = "MEM", cv = obj@enms[[1]]@parameters$cv,
                            cv.param = as.numeric(unlist(
                              strsplit(obj@enms[[1]]@parameters$cv.param,
                                       "|", fixed = TRUE))[-1]),
                            metric = obj@enms[[1]]@parameters$metric,
                            axes.metric = obj@enms[[1]]@parameters$axes.metric,
                            ensemble.metric = unlist(
                              strsplit(obj@enms[[1]]@parameters$ensemble.metric,
                                       ".", fixed = TRUE))[-1],
                            ensemble.thresh = as.numeric(unlist(
                              strsplit(obj@enms[[1]]@parameters$ensemble.thresh,
                                       "|", fixed = TRUE))[-1]),
                            uncertainty = FALSE,
                            weight = as.logical(obj@enms[[1]]@parameters$weight),
                            verbose = FALSE)
  MEM@projection <- MEM@projection*maxOcc
  return(MEM)
}