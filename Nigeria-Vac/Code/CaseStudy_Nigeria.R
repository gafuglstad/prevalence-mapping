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

## Helper functions
  source('functions.R')

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

## Direct estimates
  print("Computing direct estimates...")
  
  if(is.null(myData$householdIdx)){
    load("../Data/Nigeria_AGG_DHS/directEstimates.RData")
  } else{
    # Compute estimates
    direct.est = getNigeriaDirect(myData)
    direct.est.ur = getNigeriaDirect(myData, byUrban = TRUE)
    
    save(file = "../Data/Nigeria_AGG_DHS/directEstimates.RData", direct.est, direct.est.ur)
  }

## Smoothed direct
  print("Computing smoothed direct estimates...")
  # Set prior
  bym2prior = list(prec = list(param = c(1, 0.05)),
                   phi  = list(param = c(0.5, 0.5)))
  
  # Compute full estimate
  smooth.direct = getNigeriaSmoothDirect(direct.est = direct.est,
                                         nigeriaGraph = nigeriaGraph_admin1,
                                         bym2prior = bym2prior)
  
  # Compute hold-out estimates
  smooth.direct.holdOut = smooth.direct
  for(i in 1:dim(smooth.direct.holdOut)[1]){
    # Remove data from region i
    tmpData = direct.est
    tmpData$logitP[i] = NA
    tmpData$se[i] = 1
    
    # Fit model
    tmpRes = getNigeriaSmoothDirect(direct.est = tmpData,
                                    nigeriaGraph = nigeriaGraph_admin1,
                                    bym2prior = bym2prior)
    
    # Extract estimate
    smooth.direct.holdOut[i,] = tmpRes[i,]
  }

## Intercept + urban/rural
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


## BYM on admin1
  print("Computing BYM on admin1...")
  # Set priors
  bym2prior = list(prec = list(param = c(1, 0.05)),
                   phi  = list(param = c(0.5, 0.5)))
  iidPrior  = list(prec = list(prior = "pc.prec",
                               param = c(1, 0.05)))
  
  # Compute full estimate
  inla.admin1 = getAreaLGM(myData = myData,
                           nigeriaGraph = nigeriaGraph_admin1,
                           bym2prior = bym2prior,
                           clustPrior = iidPrior)
  
  # Calculate estimates
  admin1.bym = aggBYM_admin1(res.inla = inla.admin1,
                             popList = nigeriaPop,
                             myData = myData,
                             nSamp = 1000)
  admin1.bym$clustSum = list(cIdx = myData$clusterIdx,
                             cInf = inla.admin1$summary.linear.predictor)
  
  # Compute hold-out estimates
  admin1.bym.holdOut = admin1.bym
  admin1.bym.holdOutMarginals = inla.admin1$marginals.linear.predictor[1:nrow(myData)]
  for(i in 1:length(admin1.bym$meas$admin1)){
    print("Take out region (admin1):")
    print(i)
    # Remove data from region i
    tmpData = myData
    idx = as.numeric(myData$admin1Fac) == i
    tmpData$measles[idx] = NA
    
    # Fit model
    inla.admin1.tmp = getAreaLGM(myData = tmpData,
                                 nigeriaGraph = nigeriaGraph_admin1,
                                 bym2prior = bym2prior,
                                 clustPrior = iidPrior)
    admin1.bym.tmp = aggBYM_admin1(res.inla = inla.admin1.tmp,
                                   popList = nigeriaPop,
                                   myData = tmpData,
                                   nSamp = 1000)
    
    # Extract estimate
    admin1.bym.holdOut$meas.ur[(i-1)*2+c(1,2),]  = admin1.bym.tmp$meas.ur[(i-1)*2+c(1,2),]
    admin1.bym.holdOut$meas[i,]     = admin1.bym.tmp$meas[i,]
    admin1.bym.holdOut$real.ur[(i-1)*2+c(1,2),]  = admin1.bym.tmp$real.ur[(i-1)*2+c(1,2),]
    admin1.bym.holdOut$real[i,]     = admin1.bym.tmp$real[i,]
    admin1.bym.holdOut$overD.ur[(i-1)*2+c(1,2),] = admin1.bym.tmp$overD.ur[(i-1)*2+c(1,2),]
    admin1.bym.holdOut$overD[i,]    = admin1.bym.tmp$overD[i,]
    admin1.bym.holdOut$samples$p.overD[i,] = admin1.bym.tmp$samples$p.overD[i,]
    admin1.bym.holdOut$samples$pRur.overD[i,] = admin1.bym.tmp$samples$pRur.overD[i,]
    admin1.bym.holdOut$samples$pUrb.overD[i,] = admin1.bym.tmp$samples$pUrb.overD[i,]
    
    admin1.bym.holdOut$clustSum$cInf[idx,] = inla.admin1.tmp$summary.linear.predictor[idx,]
    
    admin1.bym.holdOutMarginals[idx] = inla.admin1.tmp$marginals.linear.predictor[idx]
  }
  save.image('everythingAdmin1.RData')

## BYM on admin2
  print("Computing BYM on admin2...")
  # Set priors
  bym2prior = list(prec = list(param = c(1, 0.05)),
                   phi  = list(param = c(0.5, 0.5)))
  iidPrior  = list(prec = list(prior = "pc.prec",
                               param = c(1, 0.05)))
  
  # Compute full estimate
  inla.admin2 = getAreaLGM(myData = myData,
                           nigeriaGraph = nigeriaGraph,
                           bym2prior = bym2prior,
                           clustPrior = iidPrior,
                           admin2 = TRUE)
  
  # Calculate estimates
  nameAdm1 = c()
  for(i in 1:774){
    nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
  }
  admin2.bym = aggBYM_admin2(res.inla = inla.admin2,
                             popList = nigeriaPop,
                             nameAdm1 = nameAdm1,
                             myData = myData,
                             nSamp = 1000)
  admin2.bym$clustSum = list(cIdx = myData$clusterIdx,
                             cInf = inla.admin2$summary.linear.predictor)
  
  # Compute hold-out estimates
  admin2.bym.holdOut = admin2.bym
  admin2.bym.holdOutMarginals = inla.admin2$marginals.linear.predictor[1:nrow(myData)]
  for(i in 1:37){
    print("Take out region (admin2):")
    print(i)
    print(Sys.time())
    # Remove data from region i
    tmpData = myData
    idx = as.numeric(myData$admin1Fac) == i
    tmpData$measles[idx] = NA
    
    # Fit model
    inla.admin2.tmp = getAreaLGM(myData = tmpData,
                                 nigeriaGraph = nigeriaGraph,
                                 bym2prior = bym2prior,
                                 clustPrior = iidPrior,
                                 admin2 = TRUE)
    admin2.bym.tmp = aggBYM_admin2(res.inla = inla.admin2.tmp,
                                   popList = nigeriaPop,
                                   nameAdm1 = nameAdm1,
                                   myData = tmpData,
                                   nSamp = 1000)
    
    # Wich admin2 have been observed
    unNameAdm1 = unique(nameAdm1)
    idxAdm2 = which(nameAdm1 == unNameAdm1[i])
    
    # Extract estimate
    admin2.bym.holdOut$overD.ur[(i-1)*2+c(1,2),]        = admin2.bym.tmp$overD.ur[(i-1)*2+c(1,2),]
    admin2.bym.holdOut$overD[i,]           = admin2.bym.tmp$overD[i,]
    for(k in idxAdm2){
      admin2.bym.holdOut$admin2.overD[k,]    = admin2.bym.tmp$admin2.overD[k,]
      admin2.bym.holdOut$admin2.overD.ur[(k-1)*2+c(1,2),] = admin2.bym.tmp$admin2.overD.ur[(k-1)*2+c(1,2),]
    }
    admin2.bym.holdOut$samples$p.overD[i,] = admin2.bym.tmp$samples$p.overD[i,]
    admin2.bym.holdOut$samples$pRur.overD[i,] = admin2.bym.tmp$samples$pRur.overD[i,]
    admin2.bym.holdOut$samples$pUrb.overD[i,] = admin2.bym.tmp$samples$pUrb.overD[i,]
    
    admin2.bym.holdOut$clustSum$cInf[idx,] = inla.admin2.tmp$summary.linear.predictor[idx,]
    admin2.bym.holdOutMarginals[idx] = inla.admin2.tmp$marginals.linear.predictor[idx]
  }

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
idxShift = sample.int(n = nrow(cvIdx), size = nrow(cvIdx))


cvResults = list(admin1 = data.frame(cIdx = newData$cIdx,
                                     mean = NA,
                                     sd = NA),
                 admin1.marg = inla.admin1$marginals.linear.predictor[1:nrow(myData)],
                 admin2 = data.frame(cIdx = newData$cIdx,
                               mean = NA,
                               sd = NA),
                 admin2.marg = inla.admin2$marginals.linear.predictor[1:nrow(myData)],
                 spde = data.frame(cIdx = newData$cIdx,
                             mean = NA,
                             sd = NA),
                 spde.marg = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                 spdeCov = data.frame(cIdx = newData$cIdx,
                                      mean = NA,
                                      sd = NA),
                 spdeCov.marg = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                 noSpace = data.frame(cIdx = newData$cIdx,
                                      mean = NA,
                                      sd = NA),
                 noSpace.marg = inla.spde$marginals.linear.predictor[1:nrow(myData)],
                 folds = cvIdx)
for(cvFold in 1:10){
  print("Take out fold:")
  print(cvFold)
  
  # Remove data from fold cvFold
  tmpData = myData
  tmpFold = cvIdx$cIdx[idxShift[which(cvIdx[,cvFold+1] == TRUE)]]
  idx = myData$clusterIdx%in%tmpFold
  tmpData$measles[idx] = NA
  
  # Fit model
  inla.admin1.tmp = getAreaLGM(myData = tmpData,
                               nigeriaGraph = nigeriaGraph_admin1,
                               bym2prior = bym2prior,
                               clustPrior = iidPrior)

  inla.admin2.tmp = getAreaLGM(myData = tmpData,
                               nigeriaGraph = nigeriaGraph,
                               bym2prior = bym2prior,
                               clustPrior = iidPrior,
                               admin2 = TRUE)

  inla.spde.tmp = getContLGM(myData = tmpData,
                             mesh = mesh,
                             prior.range = prior.range,
                             prior.sigma = prior.sigma,
                             clustPrior = iidPrior)
  inla.spdeCov.tmp = getContLGM(myData = tmpData,
                                mesh = mesh,
                                useCov = TRUE,
                                prior.range = prior.range,
                                prior.sigma = prior.sigma,
                                clustPrior = iidPrior)
  
  inlaRes = getFixedLGM(tmpData, iidPrior)

  idxNum = which(idx == TRUE)
  uniqueCl = unique(myData$clusterIdx[idxNum])
  for(i in 1:length(cvResults$admin1$cIdx)){
    if(cvResults$admin1$cIdx[i]%in%uniqueCl){
      idxRep = newData$idx[i]
      cvResults$admin1$mean[i] = inla.admin1.tmp$summary.linear.predictor$mean[idxRep]
      cvResults$admin2$mean[i] = inla.admin2.tmp$summary.linear.predictor$mean[idxRep]
      cvResults$spde$mean[i]   = inla.spde.tmp$summary.linear.predictor$mean[idxRep]
      cvResults$spdeCov$mean[i]   = inla.spdeCov.tmp$summary.linear.predictor$mean[idxRep]
      cvResults$noSpace$mean[i]   = inlaRes$summary.linear.predictor$mean[idxRep]
      
      cvResults$admin1$sd[i] = inla.admin1.tmp$summary.linear.predictor$sd[idxRep]
      cvResults$admin2$sd[i] = inla.admin2.tmp$summary.linear.predictor$sd[idxRep]
      cvResults$spde$sd[i]   = inla.spde.tmp$summary.linear.predictor$sd[idxRep]
      cvResults$spdeCov$sd[i]   = inla.spdeCov.tmp$summary.linear.predictor$sd[idxRep]
      cvResults$noSpace$sd[i]   = inlaRes$summary.linear.predictor$sd[idxRep]
      
      
      cvResults$admin1.marg[[i]] = inla.admin1.tmp$marginals.linear.predictor[[idxRep]]
      cvResults$admin2.marg[[i]] = inla.admin2.tmp$marginals.linear.predictor[[idxRep]]
      cvResults$spde.marg[[i]]   = inla.spde.tmp$marginals.linear.predictor[[idxRep]]
      cvResults$spdeCov.marg[[i]]   = inla.spdeCov.tmp$marginals.linear.predictor[[idxRep]]
      cvResults$noSpace.marg[[i]]   = inlaRes$marginals.linear.predictor[[idxRep]]
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
  # Hold-out state MSE
  mse.fixed = mean((direct.est$p_Med-fixed.holdOut$overD$p_Med)^2)
  mse.smooth = mean((direct.est$p_Med-smooth.direct.holdOut$p_Med)^2)
  mse.admin1 = mean((direct.est$p_Med-admin1.bym.holdOut$overD$p_Med)^2)
  mse.admin2 = mean((direct.est$p_Med-admin2.bym.holdOut$overD$p_Med)^2)
  mse.spde = mean((direct.est$p_Med-allLevels.spde.holdOut$overD$p_Med)^2)
  mse.spdeCov = mean((direct.est$p_Med-allLevels.spdeCov.holdOut$overD$p_Med)^2)
  
  mse.logit.fixed = mean((logit(direct.est$p_Med)-logit(fixed.holdOut$overD$p_Med))^2)
  mse.logit.smooth = mean((logit(direct.est$p_Med)-logit(smooth.direct.holdOut$p_Med))^2)
  mse.logit.admin1 = mean((logit(direct.est$p_Med)-logit(admin1.bym.holdOut$overD$p_Med))^2)
  mse.logit.admin2 = mean((logit(direct.est$p_Med)-logit(admin2.bym.holdOut$overD$p_Med))^2)
  mse.logit.spde = mean((logit(direct.est$p_Med)-logit(allLevels.spde.holdOut$overD$p_Med))^2)
  mse.logit.spdeCov = mean((logit(direct.est$p_Med)-logit(allLevels.spdeCov.holdOut$overD$p_Med))^2)
  
  mse.logit = cbind(mse.logit.fixed, mse.logit.smooth, mse.logit.admin1, mse.logit.admin2, mse.logit.spde, mse.logit.spdeCov)
  colnames(mse.logit) = c("NoSpace", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  # Mean absolute error
  mae.logit.fixed = mean(abs(logit(direct.est$p_Med)-logit(fixed.holdOut$overD$p_Med)))
  mae.logit.smooth = mean(abs(logit(direct.est$p_Med)-logit(smooth.direct.holdOut$p_Med)))
  mae.logit.admin1 = mean(abs(logit(direct.est$p_Med)-logit(admin1.bym.holdOut$overD$p_Med)))
  mae.logit.admin2 = mean(abs(logit(direct.est$p_Med)-logit(admin2.bym.holdOut$overD$p_Med)))
  mae.logit.spde = mean(abs(logit(direct.est$p_Med)-logit(allLevels.spde.holdOut$overD$p_Med)))
  mae.logit.spdeCov = mean(abs(logit(direct.est$p_Med)-logit(allLevels.spdeCov.holdOut$overD$p_Med)))
  
  mae.logit = cbind(mae.logit.fixed, mae.logit.smooth, mae.logit.admin1, mae.logit.admin2, mae.logit.spde, mae.logit.spdeCov)
  colnames(mae.logit) = c("Fixed", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  # Hold-out CRPS
  mu = rowMeans(logit(fixed.holdOut$samples$p.overD))
  stdDev = sqrt(direct.est$se^2 + apply(logit(fixed.holdOut$samples$p.overD), 1, var))
  crps.logit.fixed = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  ds.logit.fixed = mean(((direct.est$logitP-mu)/stdDev)^2 + log(stdDev^2))
  
  mu = smooth.direct.holdOut$logitP
  stdDev = sqrt(direct.est$se^2 + smooth.direct.holdOut$sd^2)
  crps.logit.smooth = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  ds.logit.smooth = mean(((direct.est$logitP-mu)/stdDev)^2 + log(stdDev^2))
  
  mu = rowMeans(logit(admin1.bym.holdOut$samples$p.overD))
  stdDev = sqrt(direct.est$se^2 + apply(logit(admin1.bym.holdOut$samples$p.overD), 1, var))
  crps.logit.admin1 = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  ds.logit.admin1 = mean(((direct.est$logitP-mu)/stdDev)^2 + log(stdDev^2))

  mu = rowMeans(logit(admin2.bym.holdOut$samples$p.overD))
  stdDev = sqrt(direct.est$se^2 + apply(logit(admin2.bym.holdOut$samples$p.overD), 1, var))
  crps.logit.admin2 = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  ds.logit.admin2 = mean(((direct.est$logitP-mu)/stdDev)^2 + log(stdDev^2))

  mu = rowMeans(logit(allLevels.spde.holdOut$samples$p.overD))
  stdDev = sqrt(direct.est$se^2 + apply(logit(allLevels.spde.holdOut$samples$p.overD), 1, var))
  crps.logit.spde = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  ds.logit.spde = mean(((direct.est$logitP-mu)/stdDev)^2 + log(stdDev^2))
  
  mu = rowMeans(logit(allLevels.spdeCov.holdOut$samples$p.overD))
  stdDev = sqrt(direct.est$se^2 + apply(logit(allLevels.spdeCov.holdOut$samples$p.overD), 1, var))
  crps.logit.spdeCov = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  ds.logit.spdeCov = mean(((direct.est$logitP-mu)/stdDev)^2 + log(stdDev^2))
  
  crps.logit = cbind(crps.logit.fixed$crps,
                     crps.logit.smooth$crps,
                     crps.logit.admin1$crps,
                     crps.logit.admin2$crps,
                     crps.logit.spde$crps,
                     crps.logit.spdeCov$crps)
  colnames(crps.logit) = c("NoSpace", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  ign.logit = cbind(crps.logit.fixed$ign,
                    crps.logit.smooth$ign,
                    crps.logit.admin1$ign,
                    crps.logit.admin2$ign,
                    crps.logit.spde$ign,
                    crps.logit.spdeCov$ign)
  colnames(ign.logit) = c("NoSpace", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  ds.logit = cbind(ds.logit.fixed,
                   ds.logit.smooth,
                   ds.logit.admin1,
                   ds.logit.admin2,
                   ds.logit.spde,
                   ds.logit.spdeCov)
  colnames(ds.logit) = c("NoSpace", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  cov.noS = sum((fixed.holdOut$overD$p_Low < direct.est$p_Med) & (fixed.holdOut$overD$p_Upp > direct.est$p_Med))/37*100
  cov.smooth = sum((smooth.direct.holdOut$p_Low < direct.est$p_Med) & (smooth.direct.holdOut$p_Upp > direct.est$p_Med))/37*100
  cov.a1 = sum((admin1.bym.holdOut$overD$p_Low < direct.est$p_Med) & (admin1.bym.holdOut$overD$p_Upp > direct.est$p_Med))/37*100
  cov.a2 = sum((admin2.bym.holdOut$overD$p_Low < direct.est$p_Med) & (admin2.bym.holdOut$overD$p_Upp > direct.est$p_Med))/37*100
  cov.cont = sum((allLevels.spde.holdOut$overD$p_Low < direct.est$p_Med) & (allLevels.spde.holdOut$overD$p_Upp > direct.est$p_Med))/37*100
  cov.covC = sum((allLevels.spdeCov.holdOut$overD$p_Low < direct.est$p_Med) & (allLevels.spdeCov.holdOut$overD$p_Upp > direct.est$p_Med))/37*100
  coverAdm1 = cbind(cov.noS,
                    cov.smooth,
                    cov.a1,
                    cov.a2,
                    cov.cont,
                    cov.covC)
  colnames(coverAdm1) = colnames(ds.logit) = c("NoSpace", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  # Full table
  scores.logit = rbind(mse.logit, mae.logit, colMeans(crps.logit), colMeans(ign.logit), ds.logit, coverAdm1)
  row.names(scores.logit) = c("MSE", "MAE", "CRPS", "Log-score", "DSS", "Coverage")
  scores.logit = t(scores.logit)
  print(round(scores.logit, 2))
  mu = rowMeans(logit(allLevels.spdeCov.holdOut$samples$p.overD))
  stdDev = sqrt(direct.est$se^2 + apply(logit(allLevels.spdeCov.holdOut$samples$p.overD), 1, var))
  crps.logit.spdeCov = verification:::crps(direct.est$logitP, cbind(mu, stdDev))
  ds.logit.spdeCov = mean(((direct.est$logitP-mu)/stdDev)^2 + log(stdDev^2))
  
  crps.logit = cbind(crps.logit.fixed$crps,
                     crps.logit.smooth$crps,
                     crps.logit.admin1$crps,
                     crps.logit.admin2$crps,
                     crps.logit.spde$crps,
                     crps.logit.spdeCov$crps)
  colnames(crps.logit) = c("NoSpace", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  ign.logit = cbind(crps.logit.fixed$ign,
                    crps.logit.smooth$ign,
                    crps.logit.admin1$ign,
                    crps.logit.admin2$ign,
                    crps.logit.spde$ign,
                    crps.logit.spdeCov$ign)
  colnames(ign.logit) = c("NoSpace", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  ds.logit = cbind(ds.logit.fixed,
                   ds.logit.smooth,
                   ds.logit.admin1,
                   ds.logit.admin2,
                   ds.logit.spde,
                   ds.logit.spdeCov)
  colnames(ds.logit) = c("NoSpace", "Smooth", "BYM (admin1)", "BYM (admin2)", "SPDE", "SPDE+Cov")
  
  # Full table
  scores.logit = rbind(mse.logit, mae.logit, colMeans(crps.logit), colMeans(ign.logit), ds.logit)
  row.names(scores.logit) = c("MSE", "MAE", "CRPS", "Log-score", "DSS")
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
                   admin1.bym$overD$p_Med,
                   admin2.bym$overD$p_Med,
                   allLevels.spde$overD$p_Med,
                   allLevels.spdeCov$overD$p_Med),
               max(direct.est$p_Med,
                   smooth.direct$p_Med,
                   admin1.bym$overD$p_Med,
                   admin2.bym$overD$p_Med,
                   allLevels.spde$overD$p_Med,
                   allLevels.spdeCov$overD$p_Med))
    colLim2 = c(min(direct.est$p_Upp-direct.est$p_Low,
                    smooth.direct$p_Upp-smooth.direct$p_Low,
                    admin1.bym$overD$p_Upp-admin1.bym$overD$p_Low,
                    admin2.bym$overD$p_Upp-admin2.bym$overD$p_Low,
                    allLevels.spde$overD$p_Upp-allLevels.spde$overD$p_Low,
                    allLevels.spdeCov$overD$p_Upp-allLevels.spdeCov$overD$p_Low),
                max(direct.est$p_Upp-direct.est$p_Low,
                    smooth.direct$p_Upp-smooth.direct$p_Low,
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
    colLim3 = c(min(admin2.bym$admin2.overD$p_Med,
                   allLevels.spde$admin2.overD$p_Med,
                   allLevels.spdeCov$admin2.overD$p_Med),
               max(admin2.bym$admin2.overD$p_Med,
                   allLevels.spde$admin2.overD$p_Med,
                   allLevels.spdeCov$admin2.overD$p_Med))
    colLim4 = c(min(admin2.bym$admin2.overD$p_Upp-admin2.bym$admin2.overD$p_Low,
                    allLevels.spde$admin2.overD$p_Upp-allLevels.spde$admin2.overD$p_Low,
                    allLevels.spdeCov$admin2.overD$p_Upp-allLevels.spdeCov$admin2.overD$p_Low),
                max(admin2.bym$admin2.overD$p_Upp-admin2.bym$admin2.overD$p_Low,
                    allLevels.spde$admin2.overD$p_Upp-allLevels.spde$admin2.overD$p_Low,
                    allLevels.spdeCov$admin2.overD$p_Upp-allLevels.spdeCov$admin2.overD$p_Low))
    
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
    
  save.image(paste('EverythingDone_', sampleFrame, '.RData', sep = ""))
  
  

  
### Figures for talk
  # Calculate predictive distributions of direct estimates
  smoothDir = matrix(NA, nrow = 37, ncol = 3)
  bymAdm1   = matrix(NA, nrow = 37, ncol = 3)
  bymAdm2   = matrix(NA, nrow = 37, ncol = 3)
  for(i in 1:37){
    # Smooth
    tmpSamp = rnorm(1000, mean = smooth.direct.holdOut$logitP[i], sd = smooth.direct.holdOut$sd[i])+rnorm(1000, mean = 0, sd = direct.est$se[i])
    pSmooth = expit(tmpSamp)
    smoothDir[i,] = quantile(pSmooth, c(0.025, 0.5, 0.975))
    
    # Admin1
    tmpSamp = logit(admin1.bym.holdOut$samples$p.overD[i,]) + rnorm(1000, mean = 0, sd = direct.est$se[i])
    pAdm1   = expit(tmpSamp)
    bymAdm1[i,] = quantile(pAdm1, c(0.025, 0.5, 0.975))
    
    # Admin1
    tmpSamp = logit(admin2.bym.holdOut$samples$p.overD[i,]) + rnorm(1000, mean = 0, sd = direct.est$se[i])
    pAdm2   = expit(tmpSamp)
    bymAdm2[i,] = quantile(pAdm2, c(0.025, 0.5, 0.975))
  }
  
  # Compare hold-out predictive distributions
  postscript('Figures/nigeria_predict_directEstimate.eps')
  idx = order(direct.est$p_Med)
  plot(NULL, NULL, ylim = c(0,1), cex = 0.5, lwd = 2.5, xlim = c(0.5, 37.5),
       xlab = "State", ylab = "Vaccination coverage")
  for(i in 1:37)
    lines(0.25 + c(i,i), c(bymAdm2[idx[i],1], bymAdm2[idx[i], 3]), col = "grey", lwd = 2)
  points(0.25+(1:37), bymAdm2[idx, 2], cex = 0.5, lwd = 2, col = "blue")
  for(i in 1:37)
    lines(-0.25 + c(i,i), c(bymAdm1[idx[i],1], bymAdm1[idx[i], 3]), col = "green", lwd = 2)
  points(-0.25+(1:37), bymAdm1[idx, 2], cex = 0.5, lwd = 2, col = "red")
#  for(i in 1:37)
#    lines(-0.25 + c(i,i), c(smooth.direct.holdOut$p_Low[idx[i]], smooth.direct.holdOut$p_Upp[idx[i]]), col = "yellow", lwd = 1.5)
#  points(-0.25+(1:37), smooth.direct.holdOut$p_Med[idx], cex = 0.5, lwd = 2, col = "blue")
  for(i in 1:37){
    lines(i + c(-0.25, 0.25), direct.est$p_Med[idx[i]]*c(1,1), lwd = 2, col = "black")
  }
  legend('topleft', legend = c("Direct", "Admin1", "Admin2"), col = c("black", "red", "blue"), pch = 1)

  setEPS()
  postscript('Figures/nigeria_predict_comp.eps')
  idx = order(crps.logit.admin2$ign)
  plot(NULL, NULL, ylim = c(-0.8,0.8), cex = 1, lwd = 2.5, xlim = c(0.5, 37.5), cex.lab = 2, cex.axis = 2,
       xlab = "State", ylab = "Error in prevalence")
  points(0.2+(1:37), bymAdm2[idx, 2]-direct.est$p_Med[idx], cex = 1.25, lwd = 3, col = "blue")
  for(i in 1:37)
    lines(0.2 + c(i,i), c(bymAdm2[idx[i],1], bymAdm2[idx[i], 3])-direct.est$p_Med[idx[i]], col = "green", lwd = 4)
  points(-0.2+(1:37), bymAdm1[idx, 2]-direct.est$p_Med[idx], cex = 1.25, lwd = 3, col = "red")
  for(i in 1:37)
    lines(-0.2 + c(i,i), c(bymAdm1[idx[i],1], bymAdm1[idx[i], 3])-direct.est$p_Med[idx[i]], col = "black", lwd = 4, lty = 1)
  abline(a = 0, b = 0, lwd = 3, col = "black")
  #  for(i in 1:37)
  #    lines(-0.25 + c(i,i), c(smooth.direct.holdOut$p_Low[idx[i]], smooth.direct.holdOut$p_Upp[idx[i]]), col = "yellow", lwd = 1.5)
  #  points(-0.25+(1:37), smooth.direct.holdOut$p_Med[idx], cex = 0.5, lwd = 2, col = "blue")
  legend('topleft', legend = c("Spat-A1", "Spat-A2"), col = c("red", "blue"), pch = 1, cex = 1.5, lty = 0, lwd = 2)
  dev.off()
  
  postscript('Figures/nigeria_predict_Scores.eps')
  plot(1:37, crps.logit.admin1$ign[idx], ylab = "Log-score", xlab = "State", lwd = 3,
       cex.axis = 2, cex.lab = 2, col = "red", ylim = c(-0.5, 5), cex = 1.5)
  points(1:37, crps.logit.admin2$ign[idx],  lwd = 2, col = "blue", cex  = 1.5)
  legend('topleft', legend = c("Spat-A1", "Spat-A2"), col = c("red", "blue"), pch = 1, cex = 1.5, lty = 0, lwd = 3)
  dev.off()

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
  
  
  
### ISBA 2021 figures
  # National averages
    ## Direct estimates
    my.svydesign <- svydesign(id= ~clusterIdx + householdIdx,
                              strata=~stratum, nest=T, 
                              weights=~weight, data=myData)
    direct.national = svyglm(measles~1,
                             design = my.svydesign,
                             family = quasibinomial)
    nat.est = data.frame(method = "direct",
                         logitP = direct.national$coefficients,
                         se = sqrt(diag(vcov(direct.national))))
    nat.est$lowQ[1] = expit(nat.est$logitP[1] - qnorm(0.975)*nat.est$se[1])
    nat.est$medQ[1] = expit(nat.est$logitP[1])
    nat.est$higQ[1] = expit(nat.est$logitP[1] + qnorm(0.975)*nat.est$se[1])
    
    ## National
    samp = fixed.model$samples$p.overD
    nat.samp = rep(0, dim(samp)[2])
    totPop = 0
    for(i in 1:37){
      # Get indicies of admin2 areas
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      totUrb = 0
      totRur = 0
      for(k in idxAdm2){
        idxUrb = nigeriaPop$idxUrb[[k]]
        currUrb = sum(nigeriaPop$popAdm2.2018[[k]][[2]][idxUrb], na.rm = TRUE)
        currRur = sum(nigeriaPop$popAdm2.2018[[k]][[2]][!idxUrb], na.rm = TRUE)
        totUrb = totUrb + currUrb
        totRur = totRur + currRur
        
        nat.samp = nat.samp + samp[i,]*(currUrb+currRur)
        totPop = totPop + (currUrb+currRur)
      }
    }
    nat.samp = nat.samp/totPop
    
    nat.est.tmp = data.frame(method = "National",
                             logitP = logit(median(nat.samp)),
                             se = sd(logit(nat.samp)),
                             lowQ = quantile(nat.samp, c(0.025)),
                             medQ = quantile(nat.samp, c(0.500)), 
                             higQ = quantile(nat.samp, c(0.975)))
    nat.est = rbind(nat.est, nat.est.tmp)
    
    ## Admin1
    samp = admin1.bym$samples$p.overD
    nat.samp = rep(0, dim(samp)[2])
    totPop = 0
    for(i in 1:37){
      # Get indicies of admin2 areas
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      totUrb = 0
      totRur = 0
      for(k in idxAdm2){
        idxUrb = nigeriaPop$idxUrb[[k]]
        currUrb = sum(nigeriaPop$popAdm2.2018[[k]][[2]][idxUrb], na.rm = TRUE)
        currRur = sum(nigeriaPop$popAdm2.2018[[k]][[2]][!idxUrb], na.rm = TRUE)
        totUrb = totUrb + currUrb
        totRur = totRur + currRur
        
        nat.samp = nat.samp + samp[i,]*(currUrb+currRur)
        totPop = totPop + (currUrb+currRur)
      }
    }
    nat.samp = nat.samp/totPop
    
    nat.est.tmp = data.frame(method = "Admin1",
                             logitP = logit(median(nat.samp)),
                             se = sd(logit(nat.samp)),
                             lowQ = quantile(nat.samp, c(0.025)),
                             medQ = quantile(nat.samp, c(0.500)), 
                             higQ = quantile(nat.samp, c(0.975)))
    nat.est = rbind(nat.est, nat.est.tmp)
  
    ## Admin2
    samp = admin2.bym$samples$p.overD
    nat.samp = rep(0, dim(samp)[2])
    totPop = 0
    for(i in 1:37){
      # Get indicies of admin2 areas
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      totUrb = 0
      totRur = 0
      for(k in idxAdm2){
        idxUrb = nigeriaPop$idxUrb[[k]]
        currUrb = sum(nigeriaPop$popAdm2.2018[[k]][[2]][idxUrb], na.rm = TRUE)
        currRur = sum(nigeriaPop$popAdm2.2018[[k]][[2]][!idxUrb], na.rm = TRUE)
        totUrb = totUrb + currUrb
        totRur = totRur + currRur
        
        nat.samp = nat.samp + samp[i,]*(currUrb+currRur)
        totPop = totPop + (currUrb+currRur)
      }
    }
    nat.samp = nat.samp/totPop
    
    nat.est.tmp = data.frame(method = "Admin2",
                             logitP = logit(median(nat.samp)),
                             se = sd(logit(nat.samp)),
                             lowQ = quantile(nat.samp, c(0.025)),
                             medQ = quantile(nat.samp, c(0.500)), 
                             higQ = quantile(nat.samp, c(0.975)))
    nat.est = rbind(nat.est, nat.est.tmp)
    
    ## SPDE
    samp = allLevels.spde$samples$p.overD
    nat.samp = rep(0, dim(samp)[2])
    totPop = 0
    popAdm1 = rep(0, 37)
    for(i in 1:37){
      # Get indicies of admin2 areas
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      totUrb = 0
      totRur = 0
      for(k in idxAdm2){
        idxUrb = nigeriaPop$idxUrb[[k]]
        currUrb = sum(nigeriaPop$popAdm2.2018[[k]][[2]][idxUrb], na.rm = TRUE)
        currRur = sum(nigeriaPop$popAdm2.2018[[k]][[2]][!idxUrb], na.rm = TRUE)
        totUrb = totUrb + currUrb
        totRur = totRur + currRur
        
        nat.samp = nat.samp + samp[i,]*(currUrb+currRur)
        totPop = totPop + (currUrb+currRur)
      }
      popAdm1[i] = totUrb+totRur
    }
    nat.samp = nat.samp/totPop
    
    nat.est.tmp = data.frame(method = "SPDE",
                             logitP = logit(median(nat.samp)),
                             se = sd(logit(nat.samp)),
                             lowQ = quantile(nat.samp, c(0.025)),
                             medQ = quantile(nat.samp, c(0.500)), 
                             higQ = quantile(nat.samp, c(0.975)))
    nat.est = rbind(nat.est, nat.est.tmp)
    
    
    
    ## SPDE+Cov
    samp = allLevels.spdeCov$samples$p.overD
    nat.samp = rep(0, dim(samp)[2])
    totPop = 0
    popAdm1 = rep(0, 37)
    for(i in 1:37){
      # Get indicies of admin2 areas
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      totUrb = 0
      totRur = 0
      for(k in idxAdm2){
        idxUrb = nigeriaPop$idxUrb[[k]]
        currUrb = sum(nigeriaPop$popAdm2.2018[[k]][[2]][idxUrb], na.rm = TRUE)
        currRur = sum(nigeriaPop$popAdm2.2018[[k]][[2]][!idxUrb], na.rm = TRUE)
        totUrb = totUrb + currUrb
        totRur = totRur + currRur
        
        nat.samp = nat.samp + samp[i,]*(currUrb+currRur)
        totPop = totPop + (currUrb+currRur)
      }
      popAdm1[i] = totUrb+totRur
    }
    nat.samp = nat.samp/totPop
    
    nat.est.tmp = data.frame(method = "SPDE+Cov",
                             logitP = logit(median(nat.samp)),
                             se = sd(logit(nat.samp)),
                             lowQ = quantile(nat.samp, c(0.025)),
                             medQ = quantile(nat.samp, c(0.500)), 
                             higQ = quantile(nat.samp, c(0.975)))
    nat.est = rbind(nat.est, nat.est.tmp)
    
    ## Figure
    setEPS()
    postscript("NationalLevel.eps")
    par(las = 2,mar = c(7,6.5,4,2)+0.1)
    plot(NULL, NULL, xlim = c(0, 5), ylim = c(51, 57), xlab = "", ylab = "Vaccination coverage\n",
         cex.axis = 2, cex.lab = 2, axes = FALSE)
    axis(1, at = 0:5, labels = c("Direct", "National", "Admin1", "Admin2", "SPDE", "SPDE+Cov"), cex.axis = 2)
    axis(2, cex.axis = 2, cex.lab = 2)
    lines(c(0,0), c(nat.est$lowQ[1], nat.est$higQ[1])*100, lwd = 5, col = "blue")
    points(0, nat.est$medQ[1]*100, lwd = 5)
    lines(c(1,1), c(nat.est$lowQ[2], nat.est$higQ[2])*100, lwd = 5, col = "green")
    points(1, nat.est$medQ[2]*100, lwd = 5)
    lines(c(2,2), c(nat.est$lowQ[3], nat.est$higQ[3])*100, lwd = 5, col = "green")
    points(2, nat.est$medQ[3]*100, lwd = 5)
    lines(c(3,3), c(nat.est$lowQ[4], nat.est$higQ[4])*100, lwd = 5, col = "green")
    points(3, nat.est$medQ[4]*100, lwd = 5)
    lines(c(4,4), c(nat.est$lowQ[5], nat.est$higQ[5])*100, lwd = 5, col = "green")
    points(4, nat.est$medQ[5]*100, lwd = 5)
    lines(c(5,5), c(nat.est$lowQ[6], nat.est$higQ[6])*100, lwd = 5, col = "green")
    points(5, nat.est$medQ[6]*100, lwd = 5)
    dev.off()
