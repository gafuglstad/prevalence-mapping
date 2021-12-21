###############################################################################
# complexModels.R                                                             #
#    Code for models and CV.                                                  #
# author: Geir-Arne Fuglstad <geir-arne.fuglstad@ntnu.no>                     #
###############################################################################

# Smooth direct model
runSmoothDirect = function(myData, nigeriaGraph_admin1, direct.est){
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
  
  return(list(res = smooth.direct,
              res.holdOut = smooth.direct.holdOut))
}

# Synthetic models on logit-scale
runSynthLogit = function(myData, nigeriaPop, nameVec, nSamp, listCov, sepUR){
  nameAdm1 = c()
  for(i in 1:774){
    nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
  }
  synthEst.logit1 = getNigeriaSyntheticLogit(myData,
                                             popList = nigeriaPop,
                                             nameAdm1 = nameAdm1,
                                             nSamp = nSamp,
                                             listCov = listCov,
                                             sepUR = sepUR)
  
  # Compute hold-out estimates
  synthEst.logit1.holdOut = synthEst.logit1
  for(i in 1:37){
    print("Take out region (Synthetic 1):")
    print(i)
    print(Sys.time())
    
    # Remove data from region i
    tmpData = myData
    idx = as.numeric(myData$admin1Fac) == i
    tmpData$measles[idx] = NA
    
    # Only estimate necessary regions
    unNameAdm1 = unique(nameAdm1)
    idxAdm2 = which(nameAdm1 == unNameAdm1[i])
    
    # Fit model
    synth.tmp = getNigeriaSyntheticLogit(tmpData,
                                         popList = nigeriaPop,
                                         nameAdm1 = nameAdm1,
                                         nSamp = 1000,
                                         listCov = listCov,
                                         onlyAdm2 = idxAdm2,
                                         onlyAdm1 = c(i),
                                         sepUR = sepUR)
    
    # Extract estimate
    synthEst.logit1.holdOut$admin1.ur[2*(i-1)+c(1,2),] = synth.tmp$admin1.ur[2*(i-1)+c(1,2),]
    synthEst.logit1.holdOut$admin1[i,]                 = synth.tmp$admin1[i,]
    for(k in idxAdm2){
      synthEst.logit1.holdOut$admin2[k,]                    = synth.tmp$admin2[k,]
      synthEst.logit1.holdOut$admin2.ur[(k-1)*2+c(1,2),]    = synth.tmp$admin2.ur[(k-1)*2+c(1,2),]
    }
    synthEst.logit1.holdOut$samples$p[i,]    = synth.tmp$samples$p[i,]
    synthEst.logit1.holdOut$samples$pRur[i,] = synth.tmp$samples$pRur[i,]
    synthEst.logit1.holdOut$samples$pUrb[i,] = synth.tmp$samples$pUrb[i,]
  }
  
  return(list(res = synthEst.logit1,
              res.holdOut = synthEst.logit1.holdOut))
}

# Synthetic models on linear scale
runSynthLinear = function(myData, nigeriaPop, nameVec, nSamp, listCov, sepUR){
  nameAdm1 = c()
  for(i in 1:774){
    nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
  }
  synthEst.linear1 = getNigeriaSyntheticLinear(myData,
                                               popList = nigeriaPop,
                                               nameAdm1 = nameAdm1,
                                               nSamp = nSamp,
                                               listCov = listCov,
                                               sepUR = sepUR)
  
  # Compute hold-out estimates
  synthEst.linear1.holdOut = synthEst.linear1
  for(i in 1:37){
    print("Take out region (Synthetic 3):")
    print(i)
    print(Sys.time())
    
    # Remove data from region i
    tmpData = myData
    idx = as.numeric(myData$admin1Fac) == i
    tmpData$measles[idx] = NA
    
    # Only estimate necessary regions
    unNameAdm1 = unique(nameAdm1)
    idxAdm2 = which(nameAdm1 == unNameAdm1[i])
    
    # Fit model
    synth.tmp = getNigeriaSyntheticLinear(tmpData,
                                          popList = nigeriaPop,
                                          nameAdm1 = nameAdm1,
                                          nSamp = 1000,
                                          listCov = listCov,
                                          onlyAdm2 = idxAdm2,
                                          onlyAdm1 = c(i),
                                          sepUR = sepUR)
    
    # Extract estimate
    synthEst.linear1.holdOut$admin1.ur[2*(i-1)+c(1,2),] = synth.tmp$admin1.ur[2*(i-1)+c(1,2),]
    synthEst.linear1.holdOut$admin1[i,]                 = synth.tmp$admin1[i,]
    for(k in idxAdm2){
      synthEst.linear1.holdOut$admin2[k,]                    = synth.tmp$admin2[k,]
      synthEst.linear1.holdOut$admin2.ur[(k-1)*2+c(1,2),]    = synth.tmp$admin2.ur[(k-1)*2+c(1,2),]
    }
    synthEst.linear1.holdOut$samples$p[i,]    = synth.tmp$samples$p[i,]
    synthEst.linear1.holdOut$samples$pRur[i,] = synth.tmp$samples$pRur[i,]
    synthEst.linear1.holdOut$samples$pUrb[i,] = synth.tmp$samples$pUrb[i,]
  }
  return(list(res = synthEst.linear1,
              res.holdOut = synthEst.linear1.holdOut))
}

runUnitLevel = function(myData,
                        nigeriaGraph = NULL,
                        popList,
                        nameVec,
                        areaPrior = NULL,
                        clustPrior, 
                        randomEffect,
                        nSamp = 1000,
                        admin2 = FALSE, 
                        covarModel = FALSE,
                        listCov){
  # Fit unit-level model
  inla.fit = getUnitLevelModel(myData = myData,
                               nigeriaGraph = nigeriaGraph,
                               areaPrior = areaPrior,
                               clustPrior = clustPrior, 
                               randomEffect = randomEffect,
                               admin2 = admin2, 
                               covarModel = covarModel)
  
  # Aggregate estimates
  nameAdm1 = c()
  for(i in 1:774){
    nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
  }
  unNameAdm1 = unique(nameAdm1)
  adm2Toadm1 = rep(NA, 774)
  for(i in 1:37){
    tmpIdx = which(nameAdm1 == unNameAdm1[i])
    adm2Toadm1[tmpIdx] = i
  }
  if(admin2){
    numAreas = 774
  }else{
    numAreas = 37
  }
  inla.agg = aggUnitLevel(res.inla = inla.fit,
                          popList = popList,
                          myData = myData,
                          nameAdm1 = nameAdm1,
                          nSamp = nSamp,
                          listCov = listCov,
                          numAreas = numAreas,
                          randEff = randomEffect,
                          covarModel = covarModel,
                          adm2Toadm1 = adm2Toadm1,
                          areaIsAdmin2 = admin2)
  
  numData = dim(myData)[1]
  inla.agg$clustSum = list(cIdx = myData$clusterIdx,
                           cInf = inla.fit$summary.linear.predictor[1:numData,])
  
  # Compute hold-out estimates
  inla.agg.houldOut = inla.agg
  inla.holdOutMarginals = inla.fit$marginals.linear.predictor[1:nrow(myData)]
  for(i in 1:37){
    print("Take out region (SPDE):")
    print(i)
    print(Sys.time())
    
    # Remove data from region i
    tmpData = myData
    idx = as.numeric(myData$admin1Fac) == i
    tmpData$measles[idx] = NA
    
    # Fit model
    inla.tmp = getUnitLevelModel(myData = tmpData,
                                 nigeriaGraph = nigeriaGraph,
                                 areaPrior = areaPrior,
                                 clustPrior = clustPrior, 
                                 randomEffect = randomEffect,
                                 admin2 = admin2, 
                                 covarModel = covarModel)
    
    # Only estimate necessary regions
    unNameAdm1 = unique(nameAdm1)
    idxAdm2 = which(nameAdm1 == unNameAdm1[i])
    inla.agg.tmp = aggUnitLevel = function(res.inla = inla.tmp,
                                           popList = popList,
                                           myData = myData,
                                           nameAdm1 = nameAdm1,
                                           nSamp = nSamp,
                                           listCov = listCov,
                                           numAreas = numAreas,
                                           randEff = randomEffect,
                                           covarModel = covarModel,
                                           adm2Toadm1 = adm2Toadm1,
                                           areaIsAdmin2 = admin2,
                                           onlyAdm2 = idxAdm2,
                                           onlyAdm1 = c(i))
    
    # Extract estimate
    inla.agg.houldOut$admin1.ur[2*(i-1)+c(1,2),]           = inla.agg.tmp$admin1.ur[2*(i-1)+c(1,2),]
    inla.agg.houldOut$admin1[i,]              = inla.agg.tmp$admin1[i,]
    for(k in idxAdm2){
      inla.agg.houldOut$admin2[k,]       = inla.agg.tmp$admin2[k,]
      inla.agg.houldOut$admin2.ur[(k-1)*2+c(1,2),]    = inla.agg.tmp$admin2.ur[(k-1)*2+c(1,2),]
    }
    inla.agg.houldOut$samples$p[i,]    = inla.agg.tmp$samples$p[i,]
    inla.agg.houldOut$samples$pRur[i,] = inla.agg.tmp$samples$pRur[i,]
    inla.agg.houldOut$samples$pUrb[i,] = inla.agg.tmp$samples$pUrb[i,]
    
    idxNum = which(idx == TRUE)
    inla.agg.houldOut$clustSum$cInf[idxNum,] = inla.tmp$summary.linear.predictor[idxNum,]
    inla.holdOutMarginals[idxNum] = inla.tmp$marginals.linear.predictor[idxNum]
  }
    
  return(list(fit = inla.fit,
              est = inla.agg,
              est.holdOut = inla.agg.houldOut,
              marginal.holdOut = inla.holdOutMarginals))
}
