###############################################################################
# caseStudy_Nigeria.R                                                         #
#    Figures and tables for case study of vaccination coverage of first dose  #
#    of measles vaccination in age group [12,24) months.                      #
###############################################################################

## Libraries
  library(rgdal)
  library(SUMMER)
  library(INLA)
  library(haven)
  library(rgeos)
  library(ggplot2)
  library(maptools)
  library(geoR)
  library(gstat)
  library(raster)
  library(verification)
  library(scoringRules)
  library(MASS)

## Helper functions
  source('functions.R')
  source('complexModels.R')

## Download Worldpop population rasters if not available on disk
  # If download times out, increase timeout limit
  options(timeout=60*10)
  if(!file.exists("../Data/Nigeria_pop/nga_f_1_2018.tif")){
    download.file(url = "https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2018/NGA/nga_f_1_2018.tif",
                  destfile = "../Data/Nigeria_pop/nga_f_1_2018.tif")
  }
  options(timeout=60*10)
  if(!file.exists("../Data/Nigeria_pop/nga_m_1_2018.tif")){
    download.file(url = "https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2018/NGA/nga_m_1_2018.tif",
                  destfile = "../Data/Nigeria_pop/nga_m_1_2018.tif")
  }
  options(timeout=60*10)
  if(!file.exists("../Data/Nigeria_pop/nga_f_1_2006.tif")){
    download.file(url = "https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2006/NGA/nga_f_1_2006.tif",
                  destfile = "../Data/Nigeria_pop/nga_f_1_2006.tif")
  }
  options(timeout=60*10)
  if(!file.exists("../Data/Nigeria_pop/nga_m_1_2006.tif")){
    download.file(url = "https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2006/NGA/nga_m_1_2006.tif",
                  destfile = "../Data/Nigeria_pop/nga_m_1_2006.tif")
  }
  options(timeout=60*10)
  if(!file.exists("../Data/Nigeria_pop/nga_ppp_2017_UNadj.tif")){
    download.file(url = "https://data.worldpop.org/GIS/Population/Global_2000_2020/2017/NGA/nga_ppp_2017_UNadj.tif",
                  destfile = "../Data/Nigeria_pop/nga_ppp_2017_UNadj.tif")
  }
  if(!file.exists("../Data/Nigeria_pop/nga_ppp_2018_UNadj.tif")){
    download.file(url = "https://data.worldpop.org/GIS/Population/Global_2000_2020/2018/NGA/nga_ppp_2018_UNadj.tif",
                  destfile = "../Data/Nigeria_pop/nga_ppp_2018_UNadj.tif")
  }
  
## The only correct option here is 2017
  sampleFrame = 2017

## Load and prepare geography
  # Read admin2 map
  if(!file.exists("../Data/gadm36_NGA_shp/gadm36_NGA_shp"))
    nigeriaMap = readOGR(dsn = "../Data/gadm36_NGA_shp",
                         layer = "gadm36_NGA_2")
  
  # Remove Lake Chad (Water body)
  nigeriaMap = nigeriaMap[-160,]
  
  # Read admin1 map
  nigeriaMap_admin1 = readOGR(dsn = "../Data/gadm36_NGA_shp",
                              layer = "gadm36_NGA_1")
  
  # Get graph
  nameVec = nigeriaMap$NAME_1
  for(i in 1:length(nameVec)){
    nameVec[i] = paste(nameVec[i], ":", nigeriaMap$NAME_2[i], sep = "")
  }
  nigeriaGraph = getAmat(nigeriaMap, names = nameVec)
  nigeriaGraph_admin1 = getAmat(nigeriaMap_admin1, nigeriaMap_admin1$NAME_1)
  
  # Make INLA graph object
  nigeriaGraphINLA = inla.read.graph(nigeriaGraph)

## Load and prepare data
  if(TRUE){
    load("../Data/Nigeria_AGG_DHS/preparedDHSdata.RData")
  } else{ #NOT-RUN: This requires registering and dowloading DHS data yourself
    # Read coordinates
    corList = readOGR(dsn = "../Data/Nigeria_GPS",
                      layer = "NGGE7BFL")
    
    # Remove points without coordinates
    idxKeep = corList$LATNUM != 0
    corList = corList[idxKeep,]
    
    # Load children's re-code
    cData = read_dta("../Data/Nigeria_children/NGKR7AFL.DTA")
    
    # Organize data
    cData = subset(cData, b5 == 1)
    myData = data.frame(clusterIdx = cData$v001, householdIdx = cData$v002,
                        stratum = cData$v023, 
                        measles = (cData$h9 == 1 | cData$h9 == 2 | cData$h9 == 3),
                        age = cData$b19,
                        weight = cData$v005/1000000)
    myData = subset(myData, age <= 23 & age >=12)
    myData = subset(myData, is.na(measles) == FALSE)
    myData$measles = myData$measles + 0
    
    # Add geographic information
    smallGeo = data.frame(clusterIdx = corList$DHSCLUST, urban = corList$URBAN_RURA,
                          lon = corList$LONGNUM, lat = corList$LATNUM,
                          admin1 = corList$ADM1NAME)
    myData = merge(myData, smallGeo, by = "clusterIdx")
    
    # Use coordinates to place in admin2
    clusters = unique(myData[,c("clusterIdx", "lon", "lat")])
    gridP = data.frame(Longitude = as.vector(clusters$lon),
                       Latitude = as.vector(clusters$lat))
    coordinates(gridP) = ~ Longitude + Latitude
    proj4string(gridP) = proj4string(nigeriaMap)
    admin2 = over(gridP, nigeriaMap)
    
    # Replace NAs with nearest area
    idx = which(is.na(admin2$NAME_2))
    n = length(idx)
    nearestArea = list()
    distToNearest = c()
    for(i in 1:n){
      gDists = gDistance(gridP[idx[i],], nigeriaMap, byid = TRUE)
      tmpArea = nigeriaMap[which.min(gDists),]
      distToNearest[i] = min(gDists)
      admin2$NAME_2[idx[i]] = tmpArea$NAME_2
      admin2$NAME_1[idx[i]] = tmpArea$NAME_1
    }
    
    # Fix names
    for(i in 1:dim(admin2)[1]){
      admin2$NAME_2[i] = paste(admin2$NAME_1[i], ":", admin2$NAME_2[i], sep = "")
    }
    clusters = list(clusterIdx = clusters$clusterIdx, admin2 = admin2$NAME_2)
    
    # Add geographic information to dataset
    myData = merge(myData, clusters, by = "clusterIdx")
    
    # Assign admin2 id consistent with map
    myData$admin2Idx = 0
    for(i in 1:length(myData$admin2Idx)){
      myData$admin2Idx[i] = which(nameVec == myData$admin2[i])
    }
    
    # Make factor
    myData$admin2Fac = factor(myData$admin2, levels = nameVec)
    myData$admin1Fac = as.factor(myData$admin1)
    
    # Add number of trials
    myData$Ntrials = 1
    
    # NB: These are removed in stored data object in repository due to the fact
    # that they are at a too fine level to be shared publicly. Household index
    # is required to do direct estimation, but not in the rest of the code
    # myData = subset(myData, select = -c(householdIdx, age))
    
    # Save
    save(file = "../Data/Nigeria_AGG_DHS/preparedDHSdata.RData", myData)
  } # END NOT-RUN
  # TODO: Take this away
  load("../Data/Nigeria_AGG_DHS/DONOTCOMMIT_preparedDHSdata.RData")
  
## Make figure for Section 2
  # Make figure
  png('Figures/Figure1a.png', width = 1200, height = 1200)
  plot(nigeriaMap, asp = 1, lwd = 4)
  points(myData$lon, myData$lat, lwd = 5, col = 'red', pch = 4, cex = 1.5)
  plot(nigeriaMap, asp = 1, lwd = 4, add = TRUE)
  par(lwd = 6)
  legend('bottomright', pch = 4, col = 'red', legend = '2018', cex = 4, pt.cex = 3,
         title = "NDHS")
  par(lwd = 1)
  dev.off()

## Extract populations and urban/rural
  print("Extracting populations...")
  nigeriaPop = getPop(myData, sampleFrame = sampleFrame)
save.image("Partial_Population.RData")


## Extract covariates
  print("Extracting covariates")
  # Poverty
  povertyCov = getCovariate(fName = "../Data/preparedCovRasters/poverty.tif",
                           obsLoc = cbind(myData$lon, myData$lat))
  myData$poverty = povertyCov$obsX
  
  # log1p(Access)
  accessCov = getCovariate(fName = "../Data/preparedCovRasters/access.tif",
                           obsLoc = cbind(myData$lon, myData$lat))
  accessCov$raster = log1p(accessCov$raster)
  myData$lAccess = log1p(accessCov$obsX)
  
  # list of covariate rasters
  listCov = list()
  listCov[[1]] = povertyCov
  listCov[[2]] = accessCov
save.image("Partial_Covariates.RData")  
################################################################################
## Direct estimates ############################################################
################################################################################
  print("Computing direct estimates...")
  
  if(is.null(myData$householdIdx)){
    load("../Data/Nigeria_AGG_DHS/directEstimates.RData")
  } else{
    # Compute estimates
    direct.est = getNigeriaDirect(myData)
    direct.est.ur = getNigeriaDirect(myData, byUrban = TRUE)
    
    save(file = "../Data/Nigeria_AGG_DHS/directEstimates.RData", direct.est, direct.est.ur)
  }

################################################################################
## Smoothed direct #############################################################
################################################################################
  print("Computing Fay-Herriot estimates...")
  ## Set-up
    # Get areal covariates
    nameAdm1 = c()
    for(i in 1:774){
      nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
    }
    areaX = computeArealCov(popList = nigeriaPop,
                            nameAdm1 = nameAdm1,
                            listCov = listCov)
    
  ## IID
    iidPrior = list(prec = list(param = c(1, 0.05),
                                prior = "pc.prec"))
    res.SmoothDirect.iid = runArealDirect(myData = myData,
                                      nigeriaGraph_admin1 = nigeriaGraph_admin1,
                                      direct.est = direct.est,
                                      rEffect = "iid",
                                      xCov = NULL,
                                      rPrior = iidPrior)
    
  ## IID + Cov
    iidPrior = list(prec = list(param = c(1, 0.05),
                                prior = "pc.prec"))
    res.SmoothDirect.iid.cov = runArealDirect(myData = myData,
                                          nigeriaGraph_admin1 = nigeriaGraph_admin1,
                                          direct.est = direct.est,
                                          rEffect = "iid",
                                          xCov = areaX$admin1.cov,
                                          rPrior = iidPrior)
    
    ## BYM
    bym2prior = list(prec = list(param = c(1, 0.05)),
                     phi  = list(param = c(0.5, 0.5)))
    res.SmoothDirect.bym = runArealDirect(myData = myData,
                                      nigeriaGraph_admin1 = nigeriaGraph_admin1,
                                      direct.est = direct.est,
                                      rEffect = "bym",
                                      xCov = NULL,
                                      rPrior = bym2prior)
    
    ## BYM + Cov
    bym2prior = list(prec = list(param = c(1, 0.05)),
                     phi  = list(param = c(0.5, 0.5)))
    res.SmoothDirect.bym.cov = runArealDirect(myData = myData,
                                          nigeriaGraph_admin1 = nigeriaGraph_admin1,
                                          direct.est = direct.est,
                                          rEffect = "bym",
                                          xCov = areaX$admin1.cov,
                                          rPrior = bym2prior)

save.image("Partial_DirectSmooth.RData")
################################################################################
## Synthetic estimates #########################################################
################################################################################
  ## Synthetic 1: Logit simple
    print("Computing sythetic estimate (logit1)...")
    res.synthLogit1 = runSynthLogit(myData, nigeriaPop, nameVec, nSamp = 1000, listCov, sepUR = FALSE)
  
  ## Synthetic 2: Logit complex
    print("Computing sythetic estimate (logit2)...")
    res.synthLogit2 = runSynthLogit(myData, nigeriaPop, nameVec, nSamp = 1000, listCov, sepUR = TRUE)
    
  ## Synthetic 3: Linear simple
    print("Computing sythetic estimate (linear1)...")
    res.synthLinear1 = runSynthLinear(myData, nigeriaPop, nameVec, nSamp = 1000, listCov, sepUR = FALSE)
    
  ## Synthetic 4: Linear complex
    print("Computing sythetic estimate (linear2)...")
    res.synthLinear2 = runSynthLinear(myData, nigeriaPop, nameVec, nSamp = 1000, listCov, sepUR = TRUE)
save.image("Partial_Synth.RData")    

#######################
## Unit-level models ##
#######################
  nameAdm1 = c()
  for(i in 1:774){
    nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
  }
  ##############################################
  ## Unit-level models: No area-level effects ##
  ##############################################
    # Covariates + cluster effect
    print("Unit-level: Covariates + cluster")
    iidPrior  = list(prec = list(prior = "pc.prec",
                                 param = c(1, 0.05)))
    res.noSpace.cov = runUnitLevel(myData = myData,
                                   clustPrior = iidPrior,
                                   randomEffect = "none",
                                   admin2 = FALSE, 
                                   covarModel = TRUE,
                                   nameVec = nameVec,
                                   popList = nigeriaPop,
                                   listCov = listCov,
                                   nSamp = 1000)
    save.image("Partial_UnitLevel.RData")
  
  ##################################################
  ## Unit-level models: admin1 area-level effects ##
  ##################################################
    ###################################################
    ## Unit-level models: iid admin1 + No covariates ##
    ###################################################
      # Intercept + IID(admin1) + cluster effect
      print("Unit-level: Intercept + IID(admin1) + cluster")
      iidPrior  = list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.05)))
      res.admin1.iid.noCov = runUnitLevel(myData = myData,
                                          clustPrior = iidPrior,
                                          areaPrior = iidPrior,
                                          randomEffect = "iid",
                                          admin2 = FALSE, 
                                          covarModel = FALSE,
                                          nameVec = nameVec,
                                          popList = nigeriaPop,
                                          listCov = listCov,
                                          nSamp = 1000)
      save.image("Partial_UnitLevel.RData")
      
    ###################################################
    ## Unit-level models: BYM admin1 + No covariates ##
    ###################################################
      # Intercept + BYM2(admin1) + cluster effect
      print("Unit-level: Intercept + BYM2(admin1) + cluster")
      iidPrior  = list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.05)))
      bym2prior = list(prec = list(param = c(1, 0.05)),
                       phi  = list(param = c(0.5, 0.5)))
      res.admin1.bym.noCov = runUnitLevel(myData = myData,
                                          clustPrior = iidPrior,
                                          areaPrior = bym2prior,
                                          nigeriaGraph = nigeriaGraph_admin1,
                                          randomEffect = "bym2",
                                          admin2 = FALSE, 
                                          covarModel = FALSE,
                                          nameVec = nameVec,
                                          popList = nigeriaPop,
                                          listCov = listCov,
                                          nSamp = 1000)
      save.image("Partial_UnitLevel.RData")
    
    ################################################
    ## Unit-level models: iid admin1 + Covariates ##
    ################################################
      # Covariates + IID(admin1) + cluster effect
      print("Unit-level: Covariates + IID(admin1) + cluster")
      iidPrior  = list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.05)))
      res.admin1.iid.cov = runUnitLevel(myData = myData,
                                        clustPrior = iidPrior,
                                        areaPrior = iidPrior,
                                        nigeriaGraph = nigeriaGraph_admin1,
                                        randomEffect = "iid",
                                        admin2 = FALSE, 
                                        covarModel = TRUE,
                                        nameVec = nameVec,
                                        popList = nigeriaPop,
                                        listCov = listCov,
                                        nSamp = 1000)
      save.image("Partial_UnitLevel.RData")
    
    ################################################
    ## Unit-level models: BYM admin1 + Covariates ##
    ################################################
      # Covariates + BYM(admin1) + cluster effect
      print("Unit-level: Covariates + BYM(admin1) + cluster")
      iidPrior  = list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.05)))
      bym2prior = list(prec = list(param = c(1, 0.05)),
                       phi  = list(param = c(0.5, 0.5)))
      res.admin1.bym.cov = runUnitLevel(myData = myData,
                                        clustPrior = iidPrior,
                                        areaPrior = bym2prior,
                                        nigeriaGraph = nigeriaGraph_admin1,
                                        randomEffect = "bym2",
                                        admin2 = FALSE, 
                                        covarModel = TRUE,
                                        nameVec = nameVec,
                                        popList = nigeriaPop,
                                        listCov = listCov,
                                        nSamp = 1000)
      save.image("Partial_UnitLevel.RData")

  ##############################################################################
  ## Unit-level models: admin2 area-level effects ##############################
  ##############################################################################
    ###################################################
    ## Unit-level models: BYM admin1 + No covariates ##
    ###################################################
      # Intercept + BYM2(admin2) + cluster effect
      print("Unit-level: Intercept + BYM2(admin2) + cluster")
      iidPrior  = list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.05)))
      bym2prior = list(prec = list(param = c(1, 0.05)),
                       phi  = list(param = c(0.5, 0.5)))
      res.admin2.bym.noCov = runUnitLevel(myData = myData,
                                          clustPrior = iidPrior,
                                          areaPrior = bym2prior,
                                          nigeriaGraph = nigeriaGraph,
                                          randomEffect = "bym2",
                                          admin2 = TRUE, 
                                          covarModel = FALSE,
                                          nameVec = nameVec,
                                          popList = nigeriaPop,
                                          listCov = listCov,
                                          nSamp = 1000)
      save.image("Partial_UnitLevel.RData")

    ################################################
    ## Unit-level models: BYM admin1 + Covariates ##
    ################################################
      # Covariates + BYM(admin2) + cluster effect
      print("Unit-level: Covariates + BYM2(admin2) + cluster")
      iidPrior  = list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.05)))
      bym2prior = list(prec = list(param = c(1, 0.05)),
                       phi  = list(param = c(0.5, 0.5)))
      res.admin2.bym.cov = runUnitLevel(myData = myData,
                                        clustPrior = iidPrior,
                                        areaPrior = bym2prior,
                                        nigeriaGraph = nigeriaGraph,
                                        randomEffect = "bym2",
                                        admin2 = TRUE, 
                                        covarModel = TRUE,
                                        nameVec = nameVec,
                                        popList = nigeriaPop,
                                        listCov = listCov,
                                        nSamp = 1000)
      save.image("Partial_UnitLevel.RData")
      
  #### Old TEST
  # Set priors
  bym2prior = list(prec = list(param = c(1, 0.05)),
                   phi  = list(param = c(0.5, 0.5)))
  iidPrior  = list(prec = list(prior = "pc.prec",
                               param = c(1, 0.05)))
  res.admin1.old = oldAdmin1Test(myData = myData,
                                 nigeriaGraph_admin1 = nigeriaGraph_admin1,
                                 nigeriaPop = nigeriaPop,
                                 iidPrior = iidPrior,
                                 bym2prior = bym2prior)
  save.image("Partial_UnitLevel.RData")
      
################################################################################
## No space and covariate model ################################################
################################################################################
    print("Computing simple model...")
    # Set priors
    iidPrior  = list(prec = list(prior = "pc.prec",
                                 param = c(1, 0.05)))
    
    # Compute full estimate
    inla.fixed = getFixedLGM(myData = myData,
                             clustPrior = iidPrior)
    
    # Calculate estimates
    fixed.model = aggFixed(res.inla = inla.fixed,
                           popList = nigeriaPop,
                           myData = myData,
                           nSamp = 1000)
    fixed.model$clustSum = list(cIdx = myData$clusterIdx,
                                cInf = inla.fixed$summary.linear.predictor)
    
    # Compute hold-out estimates
    fixed.holdOut = fixed.model
    fixed.holdOutMarginals = inla.fixed$marginals.linear.predictor[1:nrow(myData)]
    for(i in 1:37){
      print("Take out region (admin1):")
      print(i)
      # Remove data from region i
      tmpData = myData
      idx = as.numeric(myData$admin1Fac) == i
      tmpData$measles[idx] = NA
      
      # Fit model
      inla.fixed.tmp = getFixedLGM(myData = tmpData,
                                   clustPrior = iidPrior)
      fixed.model.tmp = aggFixed(res.inla = inla.fixed.tmp,
                                 popList = nigeriaPop,
                                 myData = tmpData,
                                 nSamp = 1000)
      
      # Extract estimate
      fixed.holdOut$overD.ur[(i-1)*2+c(1,2),] = fixed.model.tmp$overD.ur[(i-1)*2+c(1,2),]
      fixed.holdOut$overD[i,]                 = fixed.model.tmp$overD[i,]
      fixed.holdOut$samples$p.overD[i,]       = fixed.model.tmp$samples$p.overD[i,]
      fixed.holdOut$samples$pRur.overD[i,]    = fixed.model.tmp$samples$pRur.overD[i,]
      fixed.holdOut$samples$pUrb.overD[i,]    = fixed.model.tmp$samples$pUrb.overD[i,]
      
      fixed.holdOut$clustSum$cInf[idx,] = inla.fixed.tmp$summary.linear.predictor[idx,]
      
      fixed.holdOutMarginals[idx] = inla.fixed.tmp$marginals.linear.predictor[idx]
    }
    save.image("Partial_UnitLevel.RData")
    
################################################################################
## GRF Models ##################################################################
################################################################################
  ## SPDE
    print('Computing SPDE')
    
    # Mesh
    loc.domain = spsample(nigeriaMap, 1000, type = "random")@coords
    mesh = inla.mesh.2d(loc.domain = loc.domain, max.edge = 0.25, offset = -0.15)
    png('Figures/Nigeria_mesh.png', width = 1200, height = 800)
    plot(mesh, asp = 1, main = "")
    plot(nigeriaMap, add = TRUE, lwd = 3)
    points(myData$lon, myData$lat, col = "red")
    dev.off()
    
    # Priors
    prior.range = c(3, 0.5)
    prior.sigma = c(0.5, 0.5)
    iidPrior  = list(prec = list(prior = "pc.prec",
                                 param = c(1, 0.05)))
    
    # Fit SPDE
    inla.spde = getContLGM(myData = myData,
                           mesh = mesh,
                           prior.range = prior.range,
                           prior.sigma = prior.sigma,
                           clustPrior = iidPrior)
    
    # Aggregate estimates
    nameAdm1 = c()
    for(i in 1:774){
      nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
    }
    allLevels.spde = aggSPDE(res.inla = inla.spde,
                             popList = nigeriaPop,
                             myData = myData,
                             nameAdm1 = nameAdm1,
                             nSamp = 1000)
    numData = dim(myData)[1]
    allLevels.spde$clustSum = list(cIdx = myData$clusterIdx,
                                   cInf = inla.spde$summary.linear.predictor[1:numData,])
    
    # Compute hold-out estimates
    allLevels.spde.holdOut = allLevels.spde
    spde.holdOutMarginals = inla.spde$marginals.linear.predictor[1:nrow(myData)]
    for(i in 1:37){
      print("Take out region (SPDE):")
      print(i)
      print(Sys.time())
      
      # Remove data from region i
      tmpData = myData
      idx = as.numeric(myData$admin1Fac) == i
      tmpData$measles[idx] = NA
      
      # Fit model
      inla.spde.tmp = getContLGM(myData = tmpData,
                                 mesh = mesh,
                                 prior.range = prior.range,
                                 prior.sigma = prior.sigma,
                                 clustPrior = iidPrior)
      
      # Only estimate necessary regions
      unNameAdm1 = unique(nameAdm1)
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      allLevels.spde.tmp = aggSPDE(res.inla = inla.spde.tmp,
                                   popList = nigeriaPop,
                                   myData = tmpData,
                                   nameAdm1 = nameAdm1,
                                   nSamp = 1000,
                                   onlyAdm2 = idxAdm2,
                                   onlyAdm1 = c(i))
      
      # Extract estimate
      allLevels.spde.holdOut$overD.ur[2*(i-1)+c(1,2),]           = allLevels.spde.tmp$overD.ur[2*(i-1)+c(1,2),]
      allLevels.spde.holdOut$overD[i,]              = allLevels.spde.tmp$overD[i,]
      for(k in idxAdm2){
        allLevels.spde.holdOut$admin2.overD[k,]       = allLevels.spde.tmp$admin2.overD[k,]
        allLevels.spde.holdOut$admin2.overD.ur[(k-1)*2+c(1,2),]    = allLevels.spde.tmp$admin2.overD.ur[(k-1)*2+c(1,2),]
      }
      allLevels.spde.holdOut$samples$p.overD[i,]    = allLevels.spde.tmp$samples$p.overD[i,]
      allLevels.spde.holdOut$samples$pRur.overD[i,] = allLevels.spde.tmp$samples$pRur.overD[i,]
      allLevels.spde.holdOut$samples$pUrb.overD[i,] = allLevels.spde.tmp$samples$pUrb.overD[i,]
      
      idxNum = which(idx == TRUE)
      allLevels.spde.holdOut$clustSum$cInf[idxNum,] = inla.spde.tmp$summary.linear.predictor[idxNum,]
      spde.holdOutMarginals[idxNum] = inla.spde.tmp$marginals.linear.predictor[idxNum]
    }
    
    
    
  ## SPDE+Cov
    print('Computing SPDE+Cov')
    
    # Fit SPDE
    inla.spdeCov = getContLGM(myData = myData,
                           mesh = mesh,
                           useCov = TRUE,
                           prior.range = prior.range,
                           prior.sigma = prior.sigma,
                           clustPrior = iidPrior)
    
    # Aggregate estimates
    nameAdm1 = c()
    for(i in 1:774){
      nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
    }
    allLevels.spdeCov = aggSPDE(res.inla = inla.spdeCov,
                             popList = nigeriaPop,
                             myData = myData,
                             nameAdm1 = nameAdm1,
                             nSamp = 1000,
                             listCov = listCov)
    numData = dim(myData)[1]
    allLevels.spdeCov$clustSum = list(cIdx = myData$clusterIdx,
                                   cInf = inla.spdeCov$summary.linear.predictor[1:numData,])
    
    # Compute hold-out estimates
    allLevels.spdeCov.holdOut = allLevels.spdeCov
    spdeCov.holdOutMarginals = inla.spdeCov$marginals.linear.predictor[1:nrow(myData)]
    for(i in 1:37){
      print("Take out region (SPDE):")
      print(i)
      print(Sys.time())
      
      # Remove data from region i
      tmpData = myData
      idx = as.numeric(myData$admin1Fac) == i
      tmpData$measles[idx] = NA
      
      # Fit model
      inla.spdeCov.tmp = getContLGM(myData = tmpData,
                                 mesh = mesh,
                                 useCov = TRUE,
                                 prior.range = prior.range,
                                 prior.sigma = prior.sigma,
                                 clustPrior = iidPrior)
      
      # Only estimate necessary regions
      unNameAdm1 = unique(nameAdm1)
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      allLevels.spdeCov.tmp = aggSPDE(res.inla = inla.spdeCov.tmp,
                                   popList = nigeriaPop,
                                   myData = tmpData,
                                   nameAdm1 = nameAdm1,
                                   nSamp = 1000,
                                   onlyAdm2 = idxAdm2,
                                   onlyAdm1 = c(i),
                                   listCov = listCov)
      
      # Extract estimate
      allLevels.spdeCov.holdOut$overD.ur[2*(i-1)+c(1,2),]           = allLevels.spdeCov.tmp$overD.ur[2*(i-1)+c(1,2),]
      allLevels.spdeCov.holdOut$overD[i,]              = allLevels.spdeCov.tmp$overD[i,]
      for(k in idxAdm2){
        allLevels.spdeCov.holdOut$admin2.overD[k,]       = allLevels.spdeCov.tmp$admin2.overD[k,]
        allLevels.spdeCov.holdOut$admin2.overD.ur[(k-1)*2+c(1,2),]    = allLevels.spdeCov.tmp$admin2.overD.ur[(k-1)*2+c(1,2),]
      }
      allLevels.spdeCov.holdOut$samples$p.overD[i,]    = allLevels.spdeCov.tmp$samples$p.overD[i,]
      allLevels.spdeCov.holdOut$samples$pRur.overD[i,] = allLevels.spdeCov.tmp$samples$pRur.overD[i,]
      allLevels.spdeCov.holdOut$samples$pUrb.overD[i,] = allLevels.spdeCov.tmp$samples$pUrb.overD[i,]
      
      idxNum = which(idx == TRUE)
      allLevels.spdeCov.holdOut$clustSum$cInf[idxNum,] = inla.spdeCov.tmp$summary.linear.predictor[idxNum,]
      spdeCov.holdOutMarginals[idxNum] = inla.spdeCov.tmp$marginals.linear.predictor[idxNum]
    }
save.image("Partial_SPDE.RData")

## 10-fold cross validation
newData = data.frame(cIdx = unique(myData$clusterIdx))
for(i in 1:nrow(newData)){
  newData$y[i] = sum(myData$measles[myData$clusterIdx == newData$cIdx[i]])
  newData$n[i] = sum(myData$clusterIdx == newData$cIdx[i])
  newData$idx[i] = which.max(myData$clusterIdx == newData$cIdx[i])
}

cvIdx = data.frame(cIdx = newData$cIdx)
nCl = nrow(cvIdx)
for(i in 1:10){
  cvIdx[,i+1] = FALSE
  cvIdx[floor(nCl/10*(i-1)+1):floor(nCl/10*i),i+1] = TRUE
}
set.seed(14324)
idxShift = sample.int(n = nrow(cvIdx), size = nrow(cvIdx))


templateObj = data.frame(cIdx = newData$cIdx,
                        mean = NA,
                        sd = NA)
cvSummary = list(NoSpace = templateObj,
                 NoSpace.cov = templateObj,
                 A1.iid = templateObj,
                 A1.iid.cov = templateObj,
                 A1.bym = templateObj,
                 A1.bym.cov = templateObj,
                 A2.bym = templateObj,
                 A2.bym.cov = templateObj,
                 GRF = templateObj,
                 GRF.cov = templateObj)
cvMarginal = list(NoSpace = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  NoSpace.cov = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  A1.iid = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  A1.iid.cov = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  A1.bym = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  A1.bym.cov = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  A2.bym = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  A2.bym.cov = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  GRF = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                  GRF.cov = inla.spde$marginals.linear.predictor[1:nrow(myData)])

 
for(cvFold in 1:10){
  print("Take out fold:")
  print(cvFold)
  
  # Remove data from fold cvFold
  tmpData = myData
  tmpFold = cvIdx$cIdx[idxShift[which(cvIdx[,cvFold+1] == TRUE)]]
  idx = myData$clusterIdx%in%tmpFold
  tmpData$measles[idx] = NA
  
  # Priors
  iidPrior  = list(prec = list(prior = "pc.prec",
                               param = c(1, 0.05)))
  bym2prior = list(prec = list(param = c(1, 0.05)),
                   phi  = list(param = c(0.5, 0.5)))
  
  # No space -- no cov
  cv.inla.NoSpace = getFixedLGM(myData = tmpData,
                                clustPrior = iidPrior)
  
  # Covariates + cluster effect
  cv.NoSpace.cov = runUnitLevel(myData = tmpData,
                                clustPrior = iidPrior,
                                randomEffect = "none",
                                admin2 = FALSE, 
                                covarModel = TRUE,
                                nameVec = nameVec,
                                popList = nigeriaPop,
                                listCov = listCov,
                                nSamp = 1000,
                                holdOut = FALSE)$fit
  
  # Intercept + IID(admin1) + cluster effect
  cv.admin1.iid = runUnitLevel(myData = tmpData,
                               clustPrior = iidPrior,
                               areaPrior = iidPrior,
                               randomEffect = "iid",
                               admin2 = FALSE, 
                               covarModel = FALSE,
                               nameVec = nameVec,
                               popList = nigeriaPop,
                               listCov = listCov,
                               nSamp = 1000,
                               holdOut = FALSE)$fit
  
  # Intercept + BYM2(admin1) + cluster effect
  cv.admin1.bym = runUnitLevel(myData = tmpData,
                               clustPrior = iidPrior,
                               areaPrior = bym2prior,
                               nigeriaGraph = nigeriaGraph_admin1,
                               randomEffect = "bym2",
                               admin2 = FALSE, 
                               covarModel = FALSE,
                               nameVec = nameVec,
                               popList = nigeriaPop,
                               listCov = listCov,
                               nSamp = 1000,
                               holdOut = FALSE)$fit
  
  
  # Covariates + IID(admin1) + cluster effect
  cv.admin1.iid.cov = runUnitLevel(myData = tmpData,
                                   clustPrior = iidPrior,
                                   areaPrior = iidPrior,
                                   nigeriaGraph = nigeriaGraph_admin1,
                                   randomEffect = "iid",
                                   admin2 = FALSE, 
                                   covarModel = TRUE,
                                   nameVec = nameVec,
                                   popList = nigeriaPop,
                                   listCov = listCov,
                                   nSamp = 1000,
                                   holdOut = FALSE)$fit
  
  # Covariates + BYM(admin1) + cluster effect
  cv.admin1.bym.cov = runUnitLevel(myData = tmpData,
                                   clustPrior = iidPrior,
                                   areaPrior = bym2prior,
                                   nigeriaGraph = nigeriaGraph_admin1,
                                   randomEffect = "bym2",
                                   admin2 = FALSE, 
                                   covarModel = TRUE,
                                   nameVec = nameVec,
                                   popList = nigeriaPop,
                                   listCov = listCov,
                                   nSamp = 1000,
                                   holdOut = FALSE)$fit
  
  # Intercept + BYM2(admin2) + cluster effect
  cv.admin2.bym = runUnitLevel(myData = tmpData,
                               clustPrior = iidPrior,
                               areaPrior = bym2prior,
                               nigeriaGraph = nigeriaGraph,
                               randomEffect = "bym2",
                               admin2 = TRUE, 
                               covarModel = FALSE,
                               nameVec = nameVec,
                               popList = nigeriaPop,
                               listCov = listCov,
                               nSamp = 1000,
                               holdOut = FALSE)$fit
  
  # Covariates + BYM(admin2) + cluster effect
  cv.admin2.bym.cov = runUnitLevel(myData = tmpData,
                                   clustPrior = iidPrior,
                                   areaPrior = bym2prior,
                                   nigeriaGraph = nigeriaGraph,
                                   randomEffect = "bym2",
                                   admin2 = TRUE, 
                                   covarModel = TRUE,
                                   nameVec = nameVec,
                                   popList = nigeriaPop,
                                   listCov = listCov,
                                   nSamp = 1000,
                                   holdOut = FALSE)$fit
  
  
  # SPDE
  prior.range = c(3, 0.5)
  prior.sigma = c(0.5, 0.5)
  cv.spde = getContLGM(myData = tmpData,
                       mesh = mesh,
                       prior.range = prior.range,
                       prior.sigma = prior.sigma,
                       clustPrior = iidPrior)
  
  # SPDE + cov
  # Fit SPDE
  cv.spde.cov = getContLGM(myData = tmpData,
                           mesh = mesh,
                           useCov = TRUE,
                           prior.range = prior.range,
                           prior.sigma = prior.sigma,
                           clustPrior = iidPrior)
  
  # Extract results
  idxNum = which(idx == TRUE)
  uniqueCl = unique(myData$clusterIdx[idxNum])
  for(i in 1:length(cvResults$admin1$cIdx)){
    if(cvResults$admin1$cIdx[i]%in%uniqueCl){
      idxRep = newData$idx[i]
      cvSummary$NoSpace$mean[i] = cv.NoSpace$summary.linear.predictor$mean[idxRep]
      cvSummary$NoSpace.cov$mean[i] = cv.NoSpace.cov$summary.linear.predictor$mean[idxRep]
      cvSummary$A1.iid$mean[i] = cv.admin1.iid$summary.linear.predictor$mean[idxRep]
      cvSummary$A1.iid.cov$mean[i] = cv.admin1.iid.cov$summary.linear.predictor$mean[idxRep]
      cvSummary$A1.bym$mean[i] = cv.admin1.bym$summary.linear.predictor$mean[idxRep]
      cvSummary$A1.bym.cov$mean[i] = cv.admin1.bym.cov$summary.linear.predictor$mean[idxRep]
      cvSummary$A2.bym$mean[i] = cv.admin2.bym$summary.linear.predictor$mean[idxRep]
      cvSummary$A2.bym.cov$mean[i] = cv.admin2.bym.cov$summary.linear.predictor$mean[idxRep]
      cvSummary$GRF$mean[i] = cv.spde$summary.linear.predictor$mean[idxRep]
      cvSummary$GRF.cov$mean[i] = cv.spde.cov$summary.linear.predictor$mean[idxRep]
      
      cvSummary$NoSpace$sd[i] = cv.NoSpace$summary.linear.predictor$sd[idxRep]
      cvSummary$NoSpace.cov$sd[i] = cv.NoSpace.cov$summary.linear.predictor$sd[idxRep]
      cvSummary$A1.iid$sd[i] = cv.admin1.iid$summary.linear.predictor$sd[idxRep]
      cvSummary$A1.iid.cov$sd[i] = cv.admin1.iid.cov$summary.linear.predictor$sd[idxRep]
      cvSummary$A1.bym$sd[i] = cv.admin1.bym$summary.linear.predictor$sd[idxRep]
      cvSummary$A1.bym.cov$sd[i] = cv.admin1.bym.cov$summary.linear.predictor$sd[idxRep]
      cvSummary$A2.bym$sd[i] = cv.admin2.bym$summary.linear.predictor$sd[idxRep]
      cvSummary$A2.bym.cov$sd[i] = cv.admin2.bym.cov$summary.linear.predictor$sd[idxRep]
      cvSummary$GRF$sd[i] = cv.spde$summary.linear.predictor$sd[idxRep]
      cvSummary$GRF.cov$sd[i] = cv.spde.cov$summary.linear.predictor$sd[idxRep]
      
      cvMarginal$NoSpace[[i]] = cv.NoSpace$marginals.linear.predictor[[idxRep]]
      cvMarginal$NoSpace.cov[[i]] = cv.NoSpace.cov$marginals.linear.predictor[[idxRep]]
      cvMarginal$A1.iid[[i]] = cv.admin1.iid$marginals.linear.predictor[[idxRep]]
      cvMarginal$A1.iid.cov[[i]] = cv.admin1.iid.cov$marginals.linear.predictor[[idxRep]]
      cvMarginal$A1.bym[[i]] = cv.admin1.bym$marginals.linear.predictor[[idxRep]]
      cvMarginal$A1.bym.cov[[i]] = cv.admin1.bym.cov$marginals.linear.predictor[[idxRep]]
      cvMarginal$A2.bym[[i]] = cv.admin2.bym$marginals.linear.predictor[[idxRep]]
      cvMarginal$A2.bym.cov[[i]] = cv.admin2.bym.cov$marginals.linear.predictor[[idxRep]]
      cvMarginal$GRF[[i]] = cv.spde$marginals.linear.predictor[[idxRep]]
      cvMarginal$GRF.cov[[i]] = cv.spde.cov$marginals.linear.predictor[[idxRep]]
    }
  }
}


## Comparing different "nugget" assumptions
  # Absolute differences
  setEPS()
  postscript("Figures/compareAggAssumptions_abs.eps")
  par(mar = c(5, 4.5, 4, 2) + 0.1)
  plot(admin1.bym$overD$p_Med, admin1.bym$meas$p_Med-admin1.bym$overD$p_Med, 
       lwd = 3,
       cex.axis = 2, cex.lab = 2, cex = 1.5,
       xlab = "Overdispersion", ylab = "Absolute difference",
       bty = "l")
  abline(h = 0, lwd = 3)
  points(admin1.bym$overD$p_Med, admin1.bym$real$p_Med-admin1.bym$overD$p_Med,
         col = "red", lwd = 3, pch = 4, cex = 1.5)
  legend('bottomright', legend = c("Measurement error", "True signal"),
         col = c("black", "red"), pch = c(1, 4), cex = 1.5, lwd = 3, lty = 0)
  dev.off()
  
  postscript('Figures/compareNugAssumption_relative.eps')
  par(mar = c(5, 4.5, 4, 2) + 0.1)
  plot(admin1.bym$overD$p_Med, (admin1.bym$meas$p_Med-admin1.bym$overD$p_Med)/admin1.bym$overD$p_Med*100, 
       lwd = 3,
       cex.axis = 2, cex.lab = 2, cex = 1.5,
       xlab = "Overdispersion", ylab = "Relative difference (%)",
       bty = "l")
  abline(h = 0, lwd = 3)
  points(admin1.bym$overD$p_Med, (admin1.bym$real$p_Med-admin1.bym$overD$p_Med)/admin1.bym$overD$p_Med*100,
         col = "red", lwd = 3, pch = 4, cex = 1.5)
  legend('bottomright', legend = c("Measurement error", "True signal"),
         col = c("black", "red"), pch = c(1, 4), cex = 1.5, lwd = 3, lty = 0)
  dev.off()
  
## Scores (admin1)
  # Extract results
  noSpace.cov = res.noSpace.cov$est.holdOut
  
  admin1.iid.noCov = res.admin1.iid.noCov$est.holdOut
  admin1.iid.cov = res.admin1.iid.cov$est.holdOut
  admin1.bym.noCov = res.admin1.bym.noCov$est.holdOut
  admin1.bym.cov = res.admin1.bym.cov$est.holdOut
  admin1.old = res.admin1.old$est.holdOut
  
  admin2.bym.noCov = res.admin2.bym.noCov$est.holdOut
  admin2.bym.cov = res.admin2.bym.cov$est.holdOut
  
  synthLogit1 = res.synthLogit1$res.holdOut
  synthLogit2 = res.synthLogit2$res.holdOut
  synthLinear1 = res.synthLinear1$res.holdOut
  synthLinear2 = res.synthLinear2$res.holdOut
  
  SmoothDirect.iid = res.SmoothDirect.iid$res.holdOut
  SmoothDirect.iid.cov = res.SmoothDirect.iid.cov$res.holdOut
  SmoothDirect.bym = res.SmoothDirect.bym$res.holdOut
  SmoothDirect.bym.cov = res.SmoothDirect.bym.cov$res.holdOut
  
  # Hold-out state MSE
  mse.fixed = mean((logit(direct.est$p_Med)-logit(fixed.holdOut$overD$p_Med))^2)
  mse.noSpace.cov = mean((logit(direct.est$p_Med)-logit(noSpace.cov$admin1$p_Med))^2)
  
  mse.admin1.iid.noCov = mean((logit(direct.est$p_Med)-logit(admin1.iid.noCov$admin1$p_Med))^2)
  mse.admin1.iid.cov = mean((logit(direct.est$p_Med)-logit(admin1.iid.cov$admin1$p_Med))^2)
  mse.admin1.bym.noCov = mean((logit(direct.est$p_Med)-logit(admin1.bym.noCov$admin1$p_Med))^2)
  mse.admin1.bym.cov = mean((logit(direct.est$p_Med)-logit(admin1.bym.cov$admin1$p_Med))^2)

  mse.admin2.bym.noCov = mean((logit(direct.est$p_Med)-logit(admin2.bym.noCov$admin1$p_Med))^2)
  mse.admin2.bym.cov = mean((logit(direct.est$p_Med)-logit(admin2.bym.cov$admin1$p_Med))^2)
  
  mse.spde.noCov = mean((logit(direct.est$p_Med)-logit(allLevels.spde.holdOut$overD$p_Med))^2)
  mse.spde.cov = mean((logit(direct.est$p_Med)-logit(allLevels.spdeCov.holdOut$overD$p_Med))^2)
  
  mse.synthLogit1 = mean((logit(direct.est$p_Med)-logit(synthLogit1$admin1$p_Med))^2)
  mse.synthLogit2 = mean((logit(direct.est$p_Med)-logit(synthLogit2$admin1$p_Med))^2)
  mse.synthLinear1 = mean((logit(direct.est$p_Med)-logit(synthLinear1$admin1$p_Med))^2)
  mse.synthLinear2 = mean((logit(direct.est$p_Med)-logit(synthLinear2$admin1$p_Med))^2)
  
  mse.SmoothDirect.iid = mean((logit(direct.est$p_Med)-SmoothDirect.iid$logitP)^2)
  mse.SmoothDirect.iid.cov = mean((logit(direct.est$p_Med)-SmoothDirect.iid.cov$logitP)^2)
  mse.SmoothDirect.bym = mean((logit(direct.est$p_Med)-SmoothDirect.bym$logitP)^2)
  mse.SmoothDirect.bym.cov = mean((logit(direct.est$p_Med)-SmoothDirect.bym.cov$logitP)^2)
  
  mse.logit = cbind(mse.SmoothDirect.iid, mse.SmoothDirect.iid.cov, mse.SmoothDirect.bym, mse.SmoothDirect.bym.cov,
                    mse.synthLogit1, mse.synthLogit2, mse.synthLinear1, mse.synthLinear2,
                    mse.fixed, mse.noSpace.cov,
                    mse.admin1.iid.noCov, mse.admin1.iid.cov, mse.admin1.bym.noCov, mse.admin1.bym.cov,
                    mse.admin2.bym.noCov, mse.admin2.bym.cov,
                    mse.spde.noCov, mse.spde.cov)
  colnames(mse.logit) = c("FH-IID", "FH-IID-Cov", "FH-BYM", "FH-BYM-Cov",
                          "S-Logit", "S-Logit2", "S-Lin1", "S-Lin2",
                          "NoSpace-UR", "NoSpace-Cov",
                          "A1-IID-noCov", "A1-IID-Cov", "A1-BYM-noCov", "A1-BYM-Cov",
                          "A2-BYM-noCov", "A2-BYM-Cov",
                          "GRF-noCov", "GRF-Cov")
  
  # Hold-out CRPS
  compCRPS = function(est.holdOut, direct.est){
    mu = rowMeans(logit(est.holdOut$samples$p))
    stdDev = sqrt(direct.est$se^2 + apply(logit(est.holdOut$samples$p), 1, var))
    crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
    return(list(crps = crps.logit.fixed$CRPS,
                ign = crps.logit.fixed$IGN))
  }
  
  mu = logit(fixed.holdOut$overD$p_Med)
  stdDev = sqrt(direct.est$se^2 + apply(logit(fixed.holdOut$samples$p.overD), 1, var))
  crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  crps.fixed = list(crps = crps.logit.fixed$CRPS,
                           ign = crps.logit.fixed$IGN)
    
  crps.noSpace.cov = compCRPS(noSpace.cov, direct.est)
  
  crps.admin1.iid.noCov = compCRPS(admin1.iid.noCov, direct.est)
  crps.admin1.iid.cov = compCRPS(admin1.iid.cov, direct.est)
  crps.admin1.bym.noCov = compCRPS(admin1.bym.noCov, direct.est)
  crps.admin1.bym.cov = compCRPS(admin1.bym.cov, direct.est)
  
  crps.admin2.bym.noCov = compCRPS(admin2.bym.noCov, direct.est)
  crps.admin2.bym.cov = compCRPS(admin2.bym.cov, direct.est)
  
  mu = logit(allLevels.spde.holdOut$overD$p_Med)
  stdDev = sqrt(direct.est$se^2 + apply(logit(allLevels.spde.holdOut$samples$p.overD), 1, var))
  crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  crps.spde.noCov = list(crps = crps.logit.fixed$CRPS,
                         ign = crps.logit.fixed$IGN)
  
  mu = logit(allLevels.spdeCov.holdOut$overD$p_Med)
  stdDev = sqrt(direct.est$se^2 + apply(logit(allLevels.spdeCov.holdOut$samples$p.overD), 1, var))
  crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  crps.spde.cov = list(crps = crps.logit.fixed$CRPS,
                         ign = crps.logit.fixed$IGN)
  
  crps.synthLogit1 = compCRPS(synthLogit1, direct.est)
  crps.synthLogit2 = compCRPS(synthLogit2, direct.est)
  crps.synthLinear1 = compCRPS(synthLinear1, direct.est)
  crps.synthLinear2 = compCRPS(synthLinear2, direct.est)
  
  mu = SmoothDirect.iid$logitP
  stdDev = sqrt(direct.est$se^2 + SmoothDirect.iid$sd^2)
  crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  crps.SmoothDirect.iid = list(crps = crps.logit.fixed$CRPS,
                               ign = crps.logit.fixed$IGN)
  
  mu = SmoothDirect.iid.cov$logitP
  stdDev = sqrt(direct.est$se^2 + SmoothDirect.iid.cov$sd^2)
  crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  crps.SmoothDirect.iid.cov = list(crps = crps.logit.fixed$CRPS,
                                   ign = crps.logit.fixed$IGN)
  
  mu = SmoothDirect.bym$logitP
  stdDev = sqrt(direct.est$se^2 + SmoothDirect.bym$sd^2)
  crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  crps.SmoothDirect.bym = list(crps = crps.logit.fixed$CRPS,
                               ign = crps.logit.fixed$IGN)
  
  mu = SmoothDirect.bym.cov$logitP
  stdDev = sqrt(direct.est$se^2 + SmoothDirect.bym.cov$sd^2)
  crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  crps.SmoothDirect.bym.cov = list(crps = crps.logit.fixed$CRPS,
                                   ign = crps.logit.fixed$IGN)
  
  crps.logit = cbind(crps.SmoothDirect.iid$crps, crps.SmoothDirect.iid.cov$crps, crps.SmoothDirect.bym$crps, crps.SmoothDirect.bym.cov$crps,
                    crps.synthLogit1$crps, crps.synthLogit2$crps, crps.synthLinear1$crps, crps.synthLinear2$crps,
                    crps.fixed$crps, crps.noSpace.cov$crps,
                    crps.admin1.iid.noCov$crps, crps.admin1.iid.cov$crps, crps.admin1.bym.noCov$crps, crps.admin1.bym.cov$crps,
                    crps.admin2.bym.noCov$crps, crps.admin2.bym.cov$crps,
                    crps.spde.noCov$crps, crps.spde.cov$crps)
  colnames(crps.logit) = c("FH-IID", "FH-IID-Cov", "FH-BYM", "FH-BYM-Cov",
                           "S-Logit", "S-Logit2", "S-Lin1", "S-Lin2",
                           "NoSpace-UR", "NoSpace-Cov",
                           "A1-IID-noCov", "A1-IID-Cov", "A1-BYM-noCov", "A1-BYM-Cov",
                           "A2-BYM-noCov", "A2-BYM-Cov",
                           "GRF-noCov", "GRF-Cov")
  
  ign.logit = cbind(crps.SmoothDirect.iid$ign, crps.SmoothDirect.iid.cov$ign, crps.SmoothDirect.bym$ign, crps.SmoothDirect.bym.cov$ign,
                     crps.synthLogit1$ign, crps.synthLogit2$ign, crps.synthLinear1$ign, crps.synthLinear2$ign,
                     crps.fixed$ign, crps.noSpace.cov$ign,
                     crps.admin1.iid.noCov$ign, crps.admin1.iid.cov$ign, crps.admin1.bym.noCov$ign, crps.admin1.bym.cov$ign,
                     crps.admin2.bym.noCov$ign, crps.admin2.bym.cov$ign,
                    crps.spde.noCov$ign, crps.spde.cov$ign)
  colnames(ign.logit) = c("FH-IID", "FH-IID-Cov", "FH-BYM", "FH-BYM-Cov",
                          "S-Logit", "S-Logit2", "S-Lin1", "S-Lin2",
                          "NoSpace-UR", "NoSpace-Cov",
                          "A1-IID-noCov", "A1-IID-Cov", "A1-BYM-noCov", "A1-BYM-Cov",
                          "A2-BYM-noCov", "A2-BYM-Cov",
                          "GRF-noCov", "GRF-Cov")
  
  compCov = function(fixed.holdOut, direct.est){
    return(sum((fixed.holdOut$admin1$p_Low < direct.est$p_Med) & (fixed.holdOut$admin1$p_Upp > direct.est$p_Med))/37*100)
  }
  
  cov.fixed = sum((fixed.holdOut$overD$p_Low < direct.est$p_Med) & (fixed.holdOut$overD$p_Upp > direct.est$p_Med))/37*100
  cov.noSpace.cov = compCov(noSpace.cov, direct.est)
  
  cov.admin1.iid.noCov = compCov(admin1.iid.noCov, direct.est)
  cov.admin1.iid.cov = compCov(admin1.iid.cov, direct.est)
  cov.admin1.bym.noCov = compCov(admin1.bym.noCov, direct.est)
  cov.admin1.bym.cov = compCov(admin1.bym.cov, direct.est)
  
  cov.admin2.bym.noCov = compCov(admin2.bym.noCov, direct.est)
  cov.admin2.bym.cov = compCov(admin2.bym.cov, direct.est)
  
  cov.spde.noCov = sum((allLevels.spde.holdOut$overD$p_Low < direct.est$p_Med) & (allLevels.spde.holdOut$overD$p_Upp > direct.est$p_Med))/37*100
  cov.spde.cov = sum((allLevels.spdeCov.holdOut$overD$p_Low < direct.est$p_Med) & (allLevels.spdeCov.holdOut$overD$p_Upp > direct.est$p_Med))/37*100
  
  cov.synthLogit1 = compCov(synthLogit1, direct.est)
  cov.synthLogit2 = compCov(synthLogit2, direct.est)
  cov.synthLinear1 = compCov(synthLinear1, direct.est)
  cov.synthLinear2 = compCov(synthLinear2, direct.est)
  
  cov.SmoothDirect.iid = compCov(list(admin1 = SmoothDirect.iid), direct.est)
  cov.SmoothDirect.iid.cov = compCov(list(admin1 = SmoothDirect.iid.cov), direct.est)
  cov.SmoothDirect.bym = compCov(list(admin1 = SmoothDirect.bym), direct.est)
  cov.SmoothDirect.bym.cov = compCov(list(admin1 = SmoothDirect.bym.cov), direct.est)
  
  coverAdm1 = cbind(cov.SmoothDirect.iid, cov.SmoothDirect.iid.cov, cov.SmoothDirect.bym, cov.SmoothDirect.bym.cov,
                    cov.synthLogit1, cov.synthLogit2, cov.synthLinear1, cov.synthLinear2,
                    cov.fixed, cov.noSpace.cov,
                    cov.admin1.iid.noCov, cov.admin1.iid.cov, cov.admin1.bym.noCov, cov.admin1.bym.cov,
                    cov.admin2.bym.noCov, cov.admin2.bym.cov,
                    cov.spde.noCov, cov.spde.cov)
  colnames(coverAdm1) =  c("FH-IID", "FH-IID-Cov", "FH-BYM", "FH-BYM-Cov",
                           "S-Logit", "S-Logit2", "S-Lin1", "S-Lin2",
                           "NoSpace-UR", "NoSpace-Cov",
                           "A1-IID-noCov", "A1-IID-Cov", "A1-BYM-noCov", "A1-BYM-Cov",
                           "A2-BYM-noCov", "A2-BYM-Cov",
                           "GRF-noCov", "GRF-Cov")

  
  # Full table
  scores.logit = rbind(mse.logit, crps.logit, ign.logit, coverAdm1)
  row.names(scores.logit) = c("MSE", "CRPS", "Log-score","Coverage")
  scores.logit = t(scores.logit)
  print(round(scores.logit, 2))
  
## Scores (cluster level)
  getScores = function(y, n, inObj){
    scVec = matrix(NA, nrow = length(n), ncol = 5)
    for(i in 1:nrow(scVec)){
      etaSam = rnorm(50000, mean = inObj$mean[i], sd = inObj$sd[i])
      pSam = 1/(1+exp(-etaSam))
      yNew = rbinom(50000,size = n[i], prob = pSam)
      pSam = yNew/n[i]
      scVec[i,1] = crps_sample(y[i]/n[i], pSam)
      scVec[i,2] = (y[i]/n[i]-mean(pSam))^2/var(pSam) + log(var(pSam))
      scVec[i,3] = (y[i]/n[i]-mean(pSam))^2
      scVec[i,4] = (y[i]/n[i]-median(pSam))^2
      tmpProb = sum(y[i] == yNew)/length(yNew)
      scVec[i,5] = -log(tmpProb)
    }
    colnames(scVec) = c("CRPS", "DSS", "MSE (pred = mean)", "MSE (pred = median)", "Log-score")
    
    return(scVec)
  }
  
  getScoresMarg = function(y, n, inMarg){
    scVec = matrix(NA, nrow = length(n), ncol = 3)
    for(i in 1:nrow(scVec)){
      etaSam = inla.rmarginal(50000, inMarg[[i]])
      pSam = 1/(1+exp(-etaSam))
      yNew = rbinom(50000,size = n[i], prob = pSam)
      pSam = yNew/n[i]
      scVec[i,1] = crps_sample(y[i]/n[i], pSam)
      scVec[i,2] = (y[i]/n[i]-mean(pSam))^2
      tmpProb = sum(y[i] == yNew)/length(yNew)
      scVec[i,3] = -log(tmpProb)
    }
    colnames(scVec) = c("CRPS", "MSE (pred = mean)", "Log-score")
    
    return(scVec)
  }
  newData = data.frame(cIdx = newData$cIdx)
  for(i in 1:nrow(newData)){
    newData$y[i] = sum(myData$measles[myData$clusterIdx == newData$cIdx[i]])
    newData$n[i] = sum(myData$clusterIdx == newData$cIdx[i])
    newData$idx[i] = which.max(myData$clusterIdx == newData$cIdx[i])
  }
  admin1.cScore = getScoresMarg(newData$y, newData$n, cvResults$admin1.marg)
  admin2.cScore = getScoresMarg(newData$y, newData$n, cvResults$admin2.marg)
  spde.cScore   = getScoresMarg(newData$y, newData$n, cvResults$spde.marg)
  spdeCov.cScore   = getScoresMarg(newData$y, newData$n, cvResults$spdeCov.marg)
  noSpace.cScore   = getScoresMarg(newData$y, newData$n, cvResults$noSpace.marg)
  
  allClusterScores = rbind(colMeans(noSpace.cScore),
                           colMeans(admin1.cScore),
                           colMeans(admin2.cScore),
                           colMeans(spde.cScore),
                           colMeans(spdeCov.cScore))
  row.names(allClusterScores) = c("noSpace", "BYM (admin1)", "BYM (admin2)", "GRF", "GRF+Cov")
  print(round(allClusterScores, digits = 3))
  
## Some score plots
  # Bad places more bad with IGN scores
  plot(crps.logit.admin1$ign, crps.logit.admin2$ign)
  abline(a = 0, b = 1, col = "red")
  
  # Not so with CRPSs
  plot(crps.logit.admin1$crps, crps.logit.admin2$crps)
  abline(a = 0, b = 1, col = "red")
  
  # Worst location
# plot(xx, dnorm(xx, mean = direct.est$logitP[37], sd = direct.est$se[37]), type = "l")
#  lines(xx, dnorm(xx, mean = mu1, sqrt(apply(logit(admin1.bym.holdOut$samples$p.overD), 1, var)[37])), col = "blue")
#  lines(xx, dnorm(xx, mean = mu2, sqrt(apply(logit(admin2.bym.holdOut$samples$p.overD), 1, var)[37])), col = "red")
#  legend('topright', legend = c("direct", "admin1", "admin2"), col = c("black", "blue", "red"), lwd = 2)
#  abline(v = direct.est$logitP[37], lwd = 2)

## Make figures of estimates and uncertainty
  # Compare predictive distributions using full data
  postscript('Figures/nigeria_predDist_full.eps')
  idx = order(direct.est$p_Med)
  plot(1:37, direct.est$p_Med[idx], ylim = c(0,1), cex = 0.5, lwd = 2, xlim = c(0.5, 10.5),
       xlab = "State", ylab = "Vaccination coverage")
  for(i in 1:37)
    lines(c(i,i), c(direct.est$p_Low[idx[i]], direct.est$p_Upp[idx[i]]), col = "green", lwd = 1.5)
  points(1:37, direct.est$p_Med[idx], ylim = c(0,1), cex = 0.5, lwd = 2)
  for(i in 1:37)
    lines(0.2 + c(i,i), c(admin2.bym$overD$p_Low[idx[i]], admin2.bym$overD$p_Upp[idx[i]]), col = "blue", lwd = 1.5)
  points(0.2+(1:37), admin2.bym$overD$p_Med[idx], cex = 0.5, lwd = 2, col = "red")
  for(i in 1:37)
    lines(0.4 + c(i,i), c(admin1.bym$overD$p_Low[idx[i]], admin1.bym$overD$p_Upp[idx[i]]), col = "pink", lwd = 1.5)
  points(0.4+(1:37), admin1.bym$overD$p_Med[idx], cex = 0.5, lwd = 2, col = "grey")
  legend('topleft', legend = c("Direct", "Admin1", "Admin2"), col = c("black", "grey", "red"), pch = 1)
  for(i in 1:37)
    lines(0.6 + c(i,i), c(smooth.direct$p_Low[idx[i]], smooth.direct$p_Upp[idx[i]]), col = "yellow", lwd = 1.5)
  points(0.6+(1:37), smooth.direct$p_Med[idx], cex = 0.5, lwd = 2, col = "purple")
  legend('topleft', legend = c("Direct", "Admin1", "Admin2", "smooth"), col = c("black", "grey", "red", "purple"), pch = 1)
  dev.off()
  
  # Compare hold-out predictive distributions
  postscript('Figures/nigeria_predDist_holdout.eps')
  idx = order(direct.est$p_Med)
  plot(1:37, direct.est$p_Med[idx], ylim = c(0,1), cex = 0.5, lwd = 2, xlim = c(0.5, 10.5),
       xlab = "State", ylab = "Vaccination coverage")
  for(i in 1:37)
    lines(c(i,i), c(direct.est$p_Low[idx[i]], direct.est$p_Upp[idx[i]]), col = "green", lwd = 1.5)
  points(1:37, direct.est$p_Med[idx], ylim = c(0,1), cex = 0.5, lwd = 2)
  for(i in 1:37)
    lines(0.25 + c(i,i), c(admin2.bym.holdOut$overD$p_Low[idx[i]], admin2.bym.holdOut$overD$p_Upp[idx[i]]), col = "blue", lwd = 1.5)
  points(0.25+(1:37), admin2.bym.holdOut$overD$p_Med[idx], cex = 0.5, lwd = 2, col = "red")
  for(i in 1:37)
    lines(0.4 + c(i,i), c(admin1.bym.holdOut$overD$p_Low[idx[i]], admin1.bym.holdOut$overD$p_Upp[idx[i]]), col = "pink", lwd = 1.5)
  points(0.4+(1:37), admin1.bym.holdOut$overD$p_Med[idx], cex = 0.5, lwd = 2, col = "grey")
  for(i in 1:37)
    lines(0.6 + c(i,i), c(smooth.direct.holdOut$p_Low[idx[i]], smooth.direct.holdOut$p_Upp[idx[i]]), col = "yellow", lwd = 1.5)
  points(0.6+(1:37), smooth.direct.holdOut$p_Med[idx], cex = 0.5, lwd = 2, col = "purple")
  legend('topleft', legend = c("Direct", "Admin1", "Admin2", "smooth"), col = c("black", "grey", "red", "purple"), pch = 1)
  dev.off()
  
  ## Make admin1 plots
    # Color limits
    colLim = c(min(direct.est$p_Med,
                   smooth.direct$p_Med,
                   synthEst$admin1$p_Med,
                   admin1.bym$overD$p_Med,
                   admin2.bym$overD$p_Med,
                   allLevels.spde$overD$p_Med,
                   allLevels.spdeCov$overD$p_Med),
               max(direct.est$p_Med,
                   smooth.direct$p_Med,
                   synthEst$admin1$p_Med,
                   admin1.bym$overD$p_Med,
                   admin2.bym$overD$p_Med,
                   allLevels.spde$overD$p_Med,
                   allLevels.spdeCov$overD$p_Med))
    colLim2 = c(min(direct.est$p_Upp-direct.est$p_Low,
                    smooth.direct$p_Upp-smooth.direct$p_Low,
                    synthEst$admin1$p_Upp-synthEst$admin1$p_Low,
                    admin1.bym$overD$p_Upp-admin1.bym$overD$p_Low,
                    admin2.bym$overD$p_Upp-admin2.bym$overD$p_Low,
                    allLevels.spde$overD$p_Upp-allLevels.spde$overD$p_Low,
                    allLevels.spdeCov$overD$p_Upp-allLevels.spdeCov$overD$p_Low),
                max(direct.est$p_Upp-direct.est$p_Low,
                    smooth.direct$p_Upp-smooth.direct$p_Low,
                    synthEst$admin1$p_Upp-synthEst$admin1$p_Low,
                    admin1.bym$overD$p_Upp-admin1.bym$overD$p_Low,
                    admin2.bym$overD$p_Upp-admin2.bym$overD$p_Low,
                    allLevels.spde$overD$p_Upp-allLevels.spde$overD$p_Low,
                    allLevels.spdeCov$overD$p_Upp-allLevels.spdeCov$overD$p_Low))
    
    # Direct estimate
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_direct.png', 
                estVal = direct.est$p_Med, 
                graph = nigeriaMap_admin1,
                colLim = colLim)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_direct_unc.png', 
                estVal = direct.est$p_Upp-direct.est$p_Low, 
                graph = nigeriaMap_admin1,
                colLim = colLim2,
                leg = "CI width")
    
    # Smoothed direct
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_smoothDirect.png', 
                estVal = smooth.direct$p_Med, 
                graph = nigeriaMap_admin1,
                colLim = colLim)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_smoothDirect_unc.png', 
                estVal = smooth.direct$p_Upp-smooth.direct$p_Low, 
                graph = nigeriaMap_admin1,
                colLim = colLim2,
                leg = "CI width")
    
    # Smoothed direct
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_synthetic.png', 
                estVal = synthEst$admin1$p_Med, 
                graph = nigeriaMap_admin1,
                colLim = colLim)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_synthetic_unc.png', 
                estVal = synthEst$admin1$p_Upp-synthEst$admin1$p_Low, 
                graph = nigeriaMap_admin1,
                colLim = colLim2,
                leg = "CI width")
    
    # BYM (admin1)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_bymAdm1.png', 
                estVal = admin1.bym$overD$p_Med, 
                graph = nigeriaMap_admin1,
                colLim = colLim)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_bymAdm1_unc.png', 
                estVal = admin1.bym$overD$p_Upp-admin1.bym$overD$p_Low, 
                graph = nigeriaMap_admin1,
                colLim = colLim2,
                leg = "CI width")
    
    # BYM (admin2)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_bymAdm2.png', 
                estVal = admin2.bym$overD$p_Med, 
                graph = nigeriaMap_admin1,
                colLim = colLim)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_bymAdm2_unc.png', 
                estVal = admin2.bym$overD$p_Upp-admin2.bym$overD$p_Low, 
                graph = nigeriaMap_admin1,
                colLim = colLim2,
                leg = "CI width")
    
    # SPDE
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_spde.png', 
                estVal = allLevels.spde$overD$p_Med, 
                graph = nigeriaMap_admin1,
                colLim = colLim)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_spde_unc.png', 
                estVal = allLevels.spde$overD$p_Upp-allLevels.spde$overD$p_Low, 
                graph = nigeriaMap_admin1,
                colLim = colLim2,
                leg = "CI width")
    
    # SPDE+Cov
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_spdeCov.png', 
                estVal = allLevels.spdeCov$overD$p_Med, 
                graph = nigeriaMap_admin1,
                colLim = colLim)
    plotAreaCol(fName = 'Figures/nig_admin1_MCV1_spdeCov_unc.png', 
                estVal = allLevels.spdeCov$overD$p_Upp-allLevels.spdeCov$overD$p_Low, 
                graph = nigeriaMap_admin1,
                colLim = colLim2,
                leg = "CI width")
    
    
  ## Admin2   
    # Color limits
    colLim3 = c(min(synthEst$admin2$p_Med,
                    admin2.bym$admin2.overD$p_Med,
                    allLevels.spde$admin2.overD$p_Med,
                    allLevels.spdeCov$admin2.overD$p_Med),
                max(synthEst$admin2$p_Med,
                    admin2.bym$admin2.overD$p_Med,
                    allLevels.spde$admin2.overD$p_Med,
                    allLevels.spdeCov$admin2.overD$p_Med))
    colLim4 = c(min(synthEst$admin2$p_Upp-synthEst$admin2$p_Low,
                    admin2.bym$admin2.overD$p_Upp-admin2.bym$admin2.overD$p_Low,
                    allLevels.spde$admin2.overD$p_Upp-allLevels.spde$admin2.overD$p_Low,
                    allLevels.spdeCov$admin2.overD$p_Upp-allLevels.spdeCov$admin2.overD$p_Low),
                max(synthEst$admin2$p_Upp-synthEst$admin2$p_Low,
                    admin2.bym$admin2.overD$p_Upp-admin2.bym$admin2.overD$p_Low,
                    allLevels.spde$admin2.overD$p_Upp-allLevels.spde$admin2.overD$p_Low,
                    allLevels.spdeCov$admin2.overD$p_Upp-allLevels.spdeCov$admin2.overD$p_Low))
    
    # Synthetic
    plotAreaCol(fName = 'Figures/nig_admin2_MCV1_synthetic.png', 
                estVal = synthEst$admin2$p_Med, 
                graph = nigeriaMap,
                colLim = colLim3)
    plotAreaCol(fName = 'Figures/nig_admin2_MCV1_synthetic_unc.png', 
                estVal = synthEst$admin2$p_Upp-synthEst$admin2$p_Low, 
                graph = nigeriaMap,
                colLim = colLim4,
                leg = "CI width")
    
    # BYM (admin2)
    plotAreaCol(fName = 'Figures/nig_admin2_MCV1_bymAdm2.png', 
                estVal = admin2.bym$admin2.overD$p_Med, 
                graph = nigeriaMap,
                colLim = colLim3)
    plotAreaCol(fName = 'Figures/nig_admin2_MCV1_bymAdm2_unc.png', 
                estVal = admin2.bym$admin2.overD$p_Upp-admin2.bym$admin2.overD$p_Low, 
                graph = nigeriaMap,
                colLim = colLim4,
                leg = "CI width")
    
    # SPDE
    plotAreaCol(fName = 'Figures/nig_admin2_MCV1_spde.png', 
                estVal = allLevels.spde$admin2.overD$p_Med, 
                graph = nigeriaMap,
                colLim = colLim3)
    plotAreaCol(fName = 'Figures/nig_admin2_MCV1_spde_unc.png', 
                estVal = allLevels.spde$admin2.overD$p_Upp-allLevels.spde$admin2.overD$p_Low, 
                graph = nigeriaMap,
                colLim = colLim4,
                leg = "CI width")
    # SPDE+Cov
    plotAreaCol(fName = 'Figures/nig_admin2_MCV1_spdeCov.png', 
                estVal = allLevels.spdeCov$admin2.overD$p_Med, 
                graph = nigeriaMap,
                colLim = colLim3)
    plotAreaCol(fName = 'Figures/nig_admin2_MCV1_spdeCov_unc.png', 
                estVal = allLevels.spdeCov$admin2.overD$p_Upp-allLevels.spdeCov$admin2.overD$p_Low, 
                graph = nigeriaMap,
                colLim = colLim4,
                leg = "CI width")
    
#  save.image(paste('EverythingDone_', sampleFrame, '.RData', sep = ""))
  
## Plot of covariates
  # Poverty
  png('Figures/Nigeria_poverty.png', width = 1200, height = 1200)
  par(mar = c(7.5, 5, 4, 2)+0.1,  lwd=3)
  plot(listCov[[1]]$raster, xlab = "Longitude (°)", ylab = "Latitude (°)",
       cex.lab = 3, cex.axis = 3,
       axis.args = list(cex.axis = 2, lwd = 4),
       legend.width=3,
       legend.shrink = 0.8,
       legend.mar = 7)
  plot(nigeriaMap, add = TRUE, lwd = 3)
  dev.off()
  
  # Access
  png('Figures/Nigeria_Log1p_Access.png', width = 1200, height = 1200)
  par(mar = c(7.5, 5, 4, 2)+0.1,  lwd=3)
  xyCor = xyFromCell(listCov[[2]]$raster, cell = 1:ncell(listCov[[2]]$raster))
  gridP = data.frame(Longitude = xyCor[,1],
                     Latitude = xyCor[,2])
  coordinates(gridP) = ~ Longitude + Latitude
  proj4string(gridP) = proj4string(nigeriaMap)
  tmpRes = over(gridP, nigeriaMap)
  tmpResIdx = which(is.na(tmpRes[,1]))
  allVals = getValues(listCov[[2]]$raster)
  allVals[tmpResIdx] = NA;
  newRaster = setValues(listCov[[2]]$raster, values = allVals)
  plot(newRaster, xlab = "Longitude (°)", ylab = "Latitude (°)",
       cex.lab = 3, cex.axis = 3,
       axis.args = list(cex.axis = 2, lwd = 4),
       legend.width=3,
       legend.shrink = 0.8,
       legend.mar = 7)
  plot(nigeriaMap, add = TRUE, lwd = 3)
  dev.off()
  
  # Urbanicity
  png('Figures/Nigeria_urbanicity.png', width = 1200, height = 1200)
  urbRaster = raster("../Data/Nigeria_pop/nga_ppp_2018_UNadj.tif")
  allVals = getValues(urbRaster)
  allVals[!is.na(allVals)] = 0
  for(i in 1:774){
    idxU = nigeriaPop$popAdm2.2018[[i]][[1]][nigeriaPop$idxUrb[[i]]]
    idxR = nigeriaPop$popAdm2.2018[[i]][[1]][!nigeriaPop$idxUrb[[i]]]
    allVals[idxU] = 1
    allVals[idxR] = 0
  }
  newRaster = setValues(x = urbRaster,
                        values = allVals)
  par(mar = c(7.5, 5, 4, 2)+0.1,  lwd=3)
  plot(newRaster, xlab = "Longitude (°)", ylab = "Latitude (°)",
       cex.lab = 3, cex.axis = 3,
       axis.args = list(cex.axis = 2, lwd = 4),
       legend.width=3,
       legend.shrink = 0.8,
       legend.mar = 7)
  plot(nigeriaMap, add = TRUE, lwd = 3)
  dev.off()
