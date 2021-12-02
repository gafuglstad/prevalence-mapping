###############################################################################
# functions.R                                                                 #
#    Help functions for Nigeria MCV1 case study.                              #
# author: Geir-Arne Fuglstad <geir-arne.fuglstad@ntnu.no>                     #
###############################################################################

# Required libraries
library(survey)
library(INLA)
library(ggplot2)
library(raster)
library(readxl)


# Extract covariates from an appropriately organized raster
getCovariate = getPop = function(fName, obsLoc){
  # Covariate raster
  covRaster = raster(fName)
  
  # Extract covariate value at observation locations
  obsX = raster::extract(x = covRaster, y = obsLoc)
  
  return(list(obsX = obsX,
              raster = covRaster))
}


# Compute direct estimates
getNigeriaDirect = function(myData, byUrban = FALSE){
  # Get direct stratum-specific direct estimates
  my.svydesign <- svydesign(id= ~clusterIdx + householdIdx,
                            strata=~stratum, nest=T, 
                            weights=~weight, data=myData)
  if(byUrban){
    direct.res = svyglm(measles~(-1)+as.factor(stratum),
                        design = my.svydesign,
                        family = quasibinomial)
    strat.admin1 = rep(unique(myData$admin1), each = 2)
    urb.admin1 = rep(c("U", "R"), 37)
    direct.est = data.frame(logitP = direct.res$coefficients,
                            se = sqrt(diag(vcov(direct.res))),
                            admin1 = strat.admin1,
                            urb = urb.admin1)
  }else{
    direct.res = svyglm(measles~(-1)+admin1Fac,
                        design = my.svydesign,
                        family = quasibinomial)
    strat.admin1 = levels(myData$admin1Fac)
    direct.est = data.frame(logitP = direct.res$coefficients,
                            se = sqrt(diag(vcov(direct.res))), 
                            admin1 = strat.admin1)
  }
  
  # Get confidence intervals
  k = qnorm(0.975)
  direct.est$logitP_Upp = direct.est$logitP + k*direct.est$se
  direct.est$logitP_Low = direct.est$logitP - k*direct.est$se
  
  # Map to probability
  direct.est$p_Low = 1/(1+exp(-direct.est$logitP_Low))
  direct.est$p_Med = 1/(1+exp(-direct.est$logitP))
  direct.est$p_Upp = 1/(1+exp(-direct.est$logitP_Upp))
  
  # Sort admin1
  if(byUrban){
    idx = order(direct.est$admin1, direct.est$urb)
  }else{
    idx = order(direct.est$admin1)
  }
  direct.est = direct.est[idx,]
  
  return(direct.est)
}

# Synthetic unit-level estimation
getNigeriaSynthetic = function(myData, popList, nameAdm1, nSamp = 1000, onlyAdm2 = 1:774, onlyAdm1 = 1:37, listCov){
  ## Step 1: Design-based estimates
    # Introduce dummy variable
    myData$urbanDummy = (myData$urban == "U")+0
    
    # Get direct stratum-specific direct estimates
    my.svydesign <- svydesign(id= ~clusterIdx + householdIdx,
                              strata=~stratum, nest=T, 
                              weights=~weight, data=myData)
    
    # Fit GLM
    syntGLM.res = svyglm(measles~urban + poverty+lAccess,
                         design = my.svydesign,
                         family = quasibinomial)
    
    # Get fit
    mu = syntGLM.res$coefficients
    Sig = vcov(syntGLM.res)
    
    # Get samples
    beta = mvrnorm(n = nSamp, mu = mu, Sigma = Sig)
  
  ## Step 2: Admin2 samples
    # Iterate through admin2 areas
    admin2 = matrix(NA, nrow = 774, ncol = nSamp)
    admin2.urb = matrix(NA, nrow = 774, ncol = nSamp)
    admin2.rur = matrix(NA, nrow = 774, ncol = nSamp)
    for(i in 1:774){
      if(!(i%in%onlyAdm2)){
        next
      }
      # Urban indicies & cell numbers from raster
      idxUrb = popList$idxUrb[[i]]
      cellNum = popList$popAdm2.2018[[i]][[1]][,1]
      
      # Get x and y coordinates
      pop2018 = raster("../Data/Nigeria_pop/nga_ppp_2018_UNadj.tif")
      xyCor = xyFromCell(pop2018, cell = cellNum)
      
      # Get covariate values
      Xdesign = matrix(NA, nrow = dim(xyCor)[1], ncol = 2)
      for(j in 1:length(listCov)){
        Xdesign[,j] = raster::extract(x = listCov[[j]]$raster,
                                      y = xyCor)
      }
      Xdesign = cbind(1, idxUrb+0, Xdesign)
      
      # Get populations
      urbPop = sum(popList$popAdm2.2018[[i]][[2]][idxUrb], na.rm = TRUE)
      rurPop = sum(popList$popAdm2.2018[[i]][[2]][!idxUrb], na.rm = TRUE)
      uProp = urbPop/(urbPop+rurPop)
      
      # Calculate samples
      localSamples = Xdesign%*%t(beta)
      localSamples = 1/(1+exp(-localSamples))
      localSamples = localSamples*kronecker(matrix(1, nrow = 1, ncol = nSamp), matrix(popList$popAdm2.2018[[i]][[2]], ncol = 1))
      
      # Save in data object
      admin2[i, ] = colSums(localSamples, na.rm = TRUE)/(urbPop+rurPop)
      admin2.urb[i,] = colSums(localSamples[idxUrb,, drop = FALSE], na.rm = TRUE)/urbPop
      admin2.rur[i,] = colSums(localSamples[!idxUrb,,drop = FALSE], na.rm = TRUE)/rurPop
    }
  
  ## Step 3: Admin1 samples
    # Iterate through admin2 areas
    admin1 = matrix(0, nrow = 37, ncol = nSamp)
    admin1.urb = matrix(0, nrow = 37, ncol = nSamp)
    admin1.rur = matrix(0, nrow = 37, ncol = nSamp)
    unNameAdm1 = unique(nameAdm1)
    for(i in 1:37){
      if(!(i%in%onlyAdm1)){
        next
      }
      # Get indicies of admin2 areas
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      totUrb = 0
      totRur = 0
      for(k in idxAdm2){
        idxUrb = popList$idxUrb[[k]]
        currUrb = sum(popList$popAdm2.2018[[k]][[2]][idxUrb], na.rm = TRUE)
        currRur = sum(popList$popAdm2.2018[[k]][[2]][!idxUrb], na.rm = TRUE)
        totUrb = totUrb + currUrb
        totRur = totRur + currRur
        
        admin1[i,] = admin1[i,] + admin2[k,]*(currUrb+currRur)
        if(currUrb > 0)
          admin1.urb[i,] = admin1.urb[i,] + admin2.urb[k,]*currUrb
        if(currRur > 0)
          admin1.rur[i,] = admin1.rur[i,] + admin2.rur[k,]*currRur
      }
      admin1[i,] = admin1[i,]/(totUrb+totRur)
      admin1.urb[i,] = admin1.urb[i,]/totUrb
      admin1.rur[i,] = admin1.rur[i,]/totRur
    }
    
  ## Calculate quantiles
  pAdmin1     = matrix(0, nrow = 37, ncol = 3)
  pAdmin1.rur = matrix(0, nrow = 37, ncol = 3)
  pAdmin1.urb = matrix(0, nrow = 37, ncol = 3)
  for(i in 1:37){
    pAdmin1[i,] = quantile(admin1[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    pAdmin1.rur[i,] = quantile(admin1.rur[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    pAdmin1.urb[i,] = quantile(admin1.urb[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  
  pAdmin2     = matrix(0, nrow = 774, ncol = 3)
  pAdmin2.rur = matrix(0, nrow = 774, ncol = 3)
  pAdmin2.urb = matrix(0, nrow = 774, ncol = 3)
  for(i in 1:774){
    pAdmin2[i,] = quantile(admin2[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    pAdmin2.rur[i,] = quantile(admin2.rur[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    pAdmin2.urb[i,] = quantile(admin2.urb[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  
  ## Make result objects (admin1)
  synt1.overD.ur = data.frame(admin1 = rep(levels(myData$admin1Fac), 2),
                              urb = rep( c("R", "U"), each = 37),
                              p_Low = c(pAdmin1.rur[, 1], pAdmin1.urb[, 1]),
                              p_Med = c(pAdmin1.rur[, 2], pAdmin1.urb[, 2]),
                              p_Upp = c(pAdmin1.rur[, 3], pAdmin1.urb[, 3]))
  idx = order(synt1.overD.ur$admin1, synt1.overD.ur$urb)
  synt1.overD.ur = synt1.overD.ur[idx,]
  
  synt1.overD = data.frame(admin1 = levels(myData$admin1Fac),
                           p_Low = pAdmin1[, 1],
                           p_Med = pAdmin1[, 2],
                           p_Upp = pAdmin1[, 3])
  idx = order(synt1.overD$admin1)
  synt1.overD = synt1.overD[idx,]
  
  ## Make result objects (admin2)
  synt2.overD.ur = data.frame(admin2 = rep(levels(myData$admin2Fac), 2),
                              urb = rep( c("R", "U"), each = 774),
                              p_Low = c(pAdmin2.rur[, 1], pAdmin2.urb[, 1]),
                              p_Med = c(pAdmin2.rur[, 2], pAdmin2.urb[, 2]),
                              p_Upp = c(pAdmin2.rur[, 3], pAdmin2.urb[, 3]))
  idx = order(synt2.overD.ur$admin2, synt2.overD.ur$urb)
  synt2.overD.ur = synt2.overD.ur[idx,]
  
  synt2.overD = data.frame(admin2 = levels(myData$admin2Fac),
                           p_Low = pAdmin2[, 1],
                           p_Med = pAdmin2[, 2],
                           p_Upp = pAdmin2[, 3])
  idx = order(synt2.overD$admin2)
  synt2.overD = synt2.overD[idx,]
  
  ## Return results
  return(list(admin1.ur = synt1.overD.ur,
              admin1 = synt1.overD,
              admin2.ur = synt2.overD.ur,
              admin2 = synt2.overD,
              samples = list(p = admin1,
                             pRur = admin1.rur,
                             pUrb = admin1.urb)))
}

# Compute smoothed direct
getNigeriaSmoothDirect = function(direct.est, nigeriaGraph, bym2prior){
  # make dataobject
  smoothed.data = data.frame(y = direct.est$logitP,
                             se = direct.est$se,
                             idxAdmin1 = 1:37)
  
  # Formula
  smoothed.formula = y ~ 1 + f(idxAdmin1,
                               model = "bym2",
                               graph = nigeriaGraph,
                               hyper = bym2prior,
                               scale.model = TRUE)
  
  # Run model
  smoothed.inla.res = inla(formula = smoothed.formula,
                           family = "gaussian",
                           control.family = list(hyper = list(prec = list (initial = 0,
                                                                           fixed = TRUE))),
                           data = smoothed.data,
                           scale = 1/se^2,
                           control.predictor = list(compute = TRUE))
  
  # Extract estimates
  smooth.direct = data.frame(admin1 = direct.est$admin1,
                             logitP = smoothed.inla.res$summary.linear.predictor$`0.5quant`,
                             sd     = smoothed.inla.res$summary.linear.predictor$sd,
                             logitP_Low = smoothed.inla.res$summary.linear.predictor$`0.025quant`,
                             logitP_Upp = smoothed.inla.res$summary.linear.predictor$`0.975quant`)
  
  # Transform to probability scale
  smooth.direct$p_Med = 1/(1+exp(-smooth.direct$logitP))
  smooth.direct$p_Low = 1/(1+exp(-smooth.direct$logitP_Low))
  smooth.direct$p_Upp = 1/(1+exp(-smooth.direct$logitP_Upp))
  
  return(smooth.direct)
}

# Compute areal BYM model
getFixedLGM = function(myData, clustPrior){
  # Formula
  formula = measles ~ urban + f(clusterIdx,
                                model = "iid",
                                hyper = clustPrior)
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE, return.marginals.predictor = TRUE),
                  control.predictor = list(compute = TRUE))
  
  return(res.inla)
}

# Compute areal BYM model
getAreaLGM = function(myData, nigeriaGraph, bym2prior, clustPrior, admin2 = FALSE){
  # Formula
  if(!admin2){
    formula = measles ~ urban + f(as.numeric(admin1Fac), 
                                  model = "bym2",
                                  graph = nigeriaGraph,
                                  scale.model = TRUE) + 
      f(clusterIdx,
        model = "iid",
        hyper = clustPrior)
  }else{
    formula = measles ~ urban + f(as.numeric(admin2Fac), 
                                  model = "bym2",
                                  graph = nigeriaGraph,
                                  scale.model = TRUE) +
      f(clusterIdx,
        model = "iid",
        hyper = clustPrior)
  }
  
  # Fit
  res.inla = inla(formula = formula,
                  data = myData,
                  family = "binomial",
                  Ntrials = Ntrials,
                  control.fixed = list(prec = 1e-4, prec.intercept = 1e-4),
                  control.compute=list(config = TRUE, return.marginals.predictor = TRUE),
                  control.predictor = list(compute = TRUE))
  
  return(res.inla)
}

getContLGM = function(myData, mesh, useCov = FALSE, prior.range, prior.sigma, clustPrior){
  # SPDE
  spde = inla.spde2.pcmatern(mesh = mesh,
                             prior.range = prior.range,
                             prior.sigma = prior.sigma)
  
  # Mapping
  A = inla.spde.make.A(mesh = mesh,
                       loc = cbind(myData$lon, myData$lat))
  
  # stack
  if(!useCov){
    stk = inla.stack(data = list(measles = myData$measles,
                                 Ntrials = myData$Ntrials),
                     A = list(A,1),
                     effects = list(list(s = 1:spde$n.spde),
                                    list(urban = (myData$urban == "U")+0 ,
                                         intercept = 1,
                                         clusterIdx = myData$clusterIdx)))
  } else{
    stk = inla.stack(data = list(measles = myData$measles,
                                 Ntrials = myData$Ntrials),
                     A = list(A,1),
                     effects = list(list(s = 1:spde$n.spde),
                                    list(urban = (myData$urban == "U")+0 ,
                                         intercept = 1,
                                         poverty = myData$poverty,
                                         lAccess = myData$lAccess,
                                         clusterIdx = myData$clusterIdx)))
  }
    
  # Formula
  if(!useCov){
    formula_spde = measles ~ intercept + urban + f(s, model = spde) + f(clusterIdx, model = "iid", hyper = clustPrior) - 1
  } else{
    formula_spde = measles ~ intercept + urban + poverty + lAccess + f(s, model = spde) + f(clusterIdx, model = "iid", hyper = clustPrior) - 1
  }
  
  # Fit
  res.inla.spde = inla(formula = formula_spde,
                       data = inla.stack.data(stk),
                       family = "binomial",
                       Ntrials = Ntrials,
                       control.inla = list(int.strategy = "ccd"),
                       control.fixed = list(prec = 1e-4,
                                            prec.intercept = 1e-4),
                       control.predictor = list(A = inla.stack.A(stk), compute = TRUE),
                       control.compute = list(config = TRUE, return.marginals.predictor = TRUE))
  
  return(res.inla.spde)
}

aggSPDE = function(res.inla, popList, myData, nameAdm1, nSamp = 1000, onlyAdm2 = 1:774, onlyAdm1 = 1:37, listCov = NULL){
  ## Step 1: Extract posterior samples of quantities of interest
    # Sample from inla object
    post.sample.spde = inla.posterior.sample(n = nSamp, result = res.inla)
    
    # Extract different effects
    intIdx = res.inla$misc$configs$contents$start[5]
    urbIdx = res.inla$misc$configs$contents$start[6]
    if(!is.null(listCov)){
      covIdx = res.inla$misc$configs$contents$start[6 + (1:length(listCov))]
    }
    spaceIdx = res.inla$misc$configs$contents$start[3]:(res.inla$misc$configs$contents$start[4]-1)
    intSample = matrix(NA, nrow = 1, ncol = nSamp)
    urbSample = matrix(NA, nrow = 1, ncol = nSamp)
    if(!is.null(listCov)){
      covSample = matrix(NA, nrow = length(covIdx), ncol = nSamp)
    }
    spaceSample = matrix(NA, nrow = length(spaceIdx), ncol = nSamp)
    clustSD = matrix(NA, nrow = 1, ncol = nSamp)
    
    for(i in 1:nSamp){
      clustSD[1,i] = 1/sqrt(post.sample.spde[[i]]$hyperpar[3])
      intSample[1,i] = post.sample.spde[[i]]$latent[intIdx]
      urbSample[1,i] = post.sample.spde[[i]]$latent[urbIdx]
      if(!is.null(listCov)){
        covSample[,i] = post.sample.spde[[i]]$latent[covIdx]
      }
      spaceSample[,i] = post.sample.spde[[i]]$latent[spaceIdx]
    }

  ## Step 2: Admin2 samples
    # Iterate through admin2 areas
    admin2 = matrix(NA, nrow = 774, ncol = nSamp)
    admin2.urb = matrix(NA, nrow = 774, ncol = nSamp)
    admin2.rur = matrix(NA, nrow = 774, ncol = nSamp)
    for(i in 1:774){
      if(!(i%in%onlyAdm2)){
        next
      }
      # Urban indicies & cell numbers from raster
      idxUrb = popList$idxUrb[[i]]
      cellNum = popList$popAdm2.2018[[i]][[1]][,1]
      
      # Get x and y coordinates
      pop2018 = raster("../Data/Nigeria_pop/nga_ppp_2018_UNadj.tif")
      xyCor = xyFromCell(pop2018, cell = cellNum)
      
      # Get covariate values
      if(!is.null(listCov)){
        Xdesign = matrix(NA, nrow = dim(xyCor)[1], ncol = 2)
        for(j in 1:length(listCov)){
          Xdesign[,j] = raster::extract(x = listCov[[j]]$raster,
                                        y = xyCor)
        }
      }
      
      # Map spatial effect to these coordinates
      Amap = inla.spde.make.A(mesh = mesh, loc = xyCor)
      
      # Full resolution
      idxSub = 1:dim(Amap)[1]
      idxSub = idxSub[!is.na(popList$popAdm2.2018[[i]][[2]][idxSub])]
      if(!is.null(listCov)){
        for(j in 1:length(listCov)){
          idxSub = idxSub[!is.na(Xdesign[idxSub,j])]
        }
      }
      inSample = rep(FALSE, dim(Amap)[1])
      inSample[idxSub] = TRUE
      Amap = Amap[idxSub,]
      if(!is.null(listCov)){
        Xdesign = Xdesign[idxSub,]
      }
      
      # Get populations
      urbPop = sum(popList$popAdm2.2018[[i]][[2]][inSample & idxUrb], na.rm = TRUE)
      rurPop = sum(popList$popAdm2.2018[[i]][[2]][inSample & !idxUrb], na.rm = TRUE)
      uProp = urbPop/(urbPop+rurPop)
      
      # Calculate samples
      localSamples = Amap%*%spaceSample + kronecker(matrix(urbSample, nrow = 1, ncol = nSamp), matrix(idxUrb[idxSub]+0, ncol = 1)) + kronecker(matrix(intSample, nrow = 1, ncol = nSamp), matrix(1, nrow = length(idxSub), ncol = 1))
      localSamples = localSamples + matrix(rnorm(length(localSamples)), nrow = dim(localSamples)[1], ncol = dim(localSamples)[2])*kronecker(matrix(clustSD, nrow = 1, ncol = nSamp), matrix(1, nrow = length(idxSub), ncol = 1))
      if(!is.null(listCov)){
        localSamples = localSamples + Xdesign%*%covSample
      }
      localSamples = 1/(1+exp(-localSamples))
      localSamples = localSamples*kronecker(matrix(1, nrow = 1, ncol = nSamp), matrix(popList$popAdm2.2018[[i]][[2]][idxSub], nrow = length(idxSub), ncol = 1))
      
      # Save in data object
      admin2[i, ] = colSums(localSamples, na.rm = TRUE)/(urbPop+rurPop)
      admin2.urb[i,] = colSums(localSamples[idxUrb[inSample],, drop = FALSE], na.rm = TRUE)/urbPop
      admin2.rur[i,] = colSums(localSamples[!idxUrb[inSample],,drop = FALSE], na.rm = TRUE)/rurPop
    }
    
  ## Step 3: Admin1 samples
    # Iterate through admin2 areas
    admin1 = matrix(0, nrow = 37, ncol = nSamp)
    admin1.urb = matrix(0, nrow = 37, ncol = nSamp)
    admin1.rur = matrix(0, nrow = 37, ncol = nSamp)
    unNameAdm1 = unique(nameAdm1)
    for(i in 1:37){
      if(!(i%in%onlyAdm1)){
        next
      }
      # Get indicies of admin2 areas
      idxAdm2 = which(nameAdm1 == unNameAdm1[i])
      totUrb = 0
      totRur = 0
      for(k in idxAdm2){
        idxUrb = popList$idxUrb[[k]]
        currUrb = sum(popList$popAdm2.2018[[k]][[2]][idxUrb], na.rm = TRUE)
        currRur = sum(popList$popAdm2.2018[[k]][[2]][!idxUrb], na.rm = TRUE)
        totUrb = totUrb + currUrb
        totRur = totRur + currRur
        
        admin1[i,] = admin1[i,] + admin2[k,]*(currUrb+currRur)
        if(currUrb > 0)
          admin1.urb[i,] = admin1.urb[i,] + admin2.urb[k,]*currUrb
        if(currRur > 0)
          admin1.rur[i,] = admin1.rur[i,] + admin2.rur[k,]*currRur
      }
      admin1[i,] = admin1[i,]/(totUrb+totRur)
      admin1.urb[i,] = admin1.urb[i,]/totUrb
      admin1.rur[i,] = admin1.rur[i,]/totRur
    }
    
  ## Calculate quantiles
    pAdmin1     = matrix(0, nrow = 37, ncol = 3)
    pAdmin1.rur = matrix(0, nrow = 37, ncol = 3)
    pAdmin1.urb = matrix(0, nrow = 37, ncol = 3)
    for(i in 1:37){
      pAdmin1[i,] = quantile(admin1[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      pAdmin1.rur[i,] = quantile(admin1.rur[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      pAdmin1.urb[i,] = quantile(admin1.urb[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    }
    
    pAdmin2     = matrix(0, nrow = 774, ncol = 3)
    pAdmin2.rur = matrix(0, nrow = 774, ncol = 3)
    pAdmin2.urb = matrix(0, nrow = 774, ncol = 3)
    for(i in 1:774){
      pAdmin2[i,] = quantile(admin2[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      pAdmin2.rur[i,] = quantile(admin2.rur[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
      pAdmin2.urb[i,] = quantile(admin2.urb[i,], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    }
    
 ## Make result objects (admin1)
    spat1.overD.ur = data.frame(admin1 = rep(levels(myData$admin1Fac), 2),
                                urb = rep( c("R", "U"), each = 37),
                                p_Low = c(pAdmin1.rur[, 1], pAdmin1.urb[, 1]),
                                p_Med = c(pAdmin1.rur[, 2], pAdmin1.urb[, 2]),
                                p_Upp = c(pAdmin1.rur[, 3], pAdmin1.urb[, 3]))
    idx = order(spat1.overD.ur$admin1, spat1.overD.ur$urb)
    spat1.overD.ur = spat1.overD.ur[idx,]
    
    spat1.overD = data.frame(admin1 = levels(myData$admin1Fac),
                             p_Low = pAdmin1[, 1],
                             p_Med = pAdmin1[, 2],
                             p_Upp = pAdmin1[, 3])
    idx = order(spat1.overD$admin1)
    spat1.overD = spat1.overD[idx,]
  
  ## Make result objects (admin2)
    spat2.overD.ur = data.frame(admin2 = rep(levels(myData$admin2Fac), 2),
                                urb = rep( c("R", "U"), each = 774),
                                p_Low = c(pAdmin2.rur[, 1], pAdmin2.urb[, 1]),
                                p_Med = c(pAdmin2.rur[, 2], pAdmin2.urb[, 2]),
                                p_Upp = c(pAdmin2.rur[, 3], pAdmin2.urb[, 3]))
    idx = order(spat2.overD.ur$admin2, spat2.overD.ur$urb)
    spat2.overD.ur = spat2.overD.ur[idx,]
    
    spat2.overD = data.frame(admin2 = levels(myData$admin2Fac),
                             p_Low = pAdmin2[, 1],
                             p_Med = pAdmin2[, 2],
                             p_Upp = pAdmin2[, 3])
    idx = order(spat2.overD$admin2)
    spat2.overD = spat2.overD[idx,]
    
  ## Return results
    return(list(overD.ur = spat1.overD.ur,
                overD = spat1.overD,
                admin2.overD.ur = spat2.overD.ur,
                admin2.overD = spat2.overD,
                samples = list(p.overD = admin1,
                               pRur.overD = admin1.rur,
                               pUrb.overD = admin1.urb)))
}

# Aggregate fixed effects model
aggFixed = function(res.inla, popList, myData, nSamp = 1000){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
 
  # Initialize storage
  etaUrb.overD = matrix(0, nrow = nSamp, ncol = 37)
  etaRur.overD = matrix(0, nrow = nSamp, ncol = 37)
  eta.overD    = matrix(0, nrow = nSamp, ncol = 37)
  
  # Get indicies in the sample
  intIdx  = res.inla$misc$configs$contents$start[3]
  urbIdx  = res.inla$misc$configs$contents$start[4]
  
  # Sample urban and rural
  nugSimStd = rnorm(1e5, mean = 0, sd = 1)
  for(i in 1:dim(etaUrb.overD)[1]){
    for(j in 1:37){
      # Nugget is measurement error
      etaUrb.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[urbIdx]
      etaRur.tmp = post.sample[[i]]$latent[intIdx]
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[1])
      etaUrb.overD[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      etaRur.overD[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))
    }
  }
  
  # Calculate combined samples
  for(j in 1:37){
    if(popList$prop.2018.child.admin1[j] < 1){
      eta.overD[,j] = logit(popList$prop.2018.child.admin1[j]*expit(etaUrb.overD[,j]) + (1-popList$prop.2018.child.admin1[j])*expit(etaRur.overD[,j]))
    } else{
      eta.overD[,j] = etaUrb.overD[,j]
    }
  }
  
  # Calculate quantiles
  etaQuantUrb.overD = matrix(0, nrow = 37, ncol = 3)
  etaQuantRur.overD = matrix(0, nrow = 37, ncol = 3)
  etaQuant.overD = matrix(0, nrow = 37, ncol = 3)
  for(i in 1:dim(etaUrb.overD)[2]){
    etaQuantUrb.overD[i,] = quantile(etaUrb.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuantRur.overD[i,] = quantile(etaRur.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuant.overD[i,]    = quantile(eta.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  pUrb.overD = 1/(1+exp(-etaQuantUrb.overD))
  pRur.overD = 1/(1+exp(-etaQuantUrb.overD))
  p.overD = 1/(1+exp(-etaQuant.overD))
  
  # OverD
  spat1.overD.ur = data.frame(admin1 = rep(levels(myData$admin1Fac), 2),
                              urb = rep( c("R", "U"), each = 37),
                              p_Low = c(pRur.overD[, 1], pUrb.overD[, 1]),
                              p_Med = c(pRur.overD[, 2], pUrb.overD[, 2]),
                              p_Upp = c(pRur.overD[, 3], pUrb.overD[, 3]))
  idx = order(spat1.overD.ur$admin1, spat1.overD.ur$urb)
  spat1.overD.ur = spat1.overD.ur[idx,]
  
  spat1.overD = data.frame(admin1 = levels(myData$admin1Fac),
                           p_Low = p.overD[, 1],
                           p_Med = p.overD[, 2],
                           p_Upp = p.overD[, 3])
  idx = order(spat1.overD$admin1)
  spat1.overD = spat1.overD[idx,]
  
  return(list(overD.ur = spat1.overD.ur,
              overD = spat1.overD,
              samples = list(p.overD = 1/(1+exp(-t(eta.overD))),
                             pRur.overD = 1/(1+exp(-t(etaRur.overD))),
                             pUrb.overD = 1/(1+exp(-t(etaUrb.overD))))))
}


# Aggregate on Admin1 level
aggBYM_admin1 = function(res.inla, popList, myData, nSamp = 1000){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Numbers of EAs in strata
  clustNumUrb = c(2452, 2006, 5492, 11715, 2008, 5126, 3949,
                  2820, 2761, 7798, 1955, 1657, 3053,
                  2293, 9529, 16957, 6874, 2621, 2548, 3090,
                  2106, 18409, 11911, 9774, 10006,
                  908, 2628, 1410, 9008, 7964, 12480,
                  9438, 25424, 7085, 8588, 19810, 22405)
  clustNumRur = c(1138, 20850, 10354, 4556, 7211, 18319, 11930,
                  9988, 17124, 16288, 7539, 8943, 11870,
                  18900, 12263, 19402, 26442, 14020, 10231, 13942,
                  9463, 3498, 1977, 4223, 9567,
                  16205, 6379, 14912, 9201, 4829, 12381,
                  2123, 0, 7408, 10667, 6097, 8701)
  clustNum = c(3590, 22856, 15846, 16271, 9219,23445,15879,
               12808,19885, 24086, 9494,10600,14923,
               21193, 21792, 36359, 33316, 16641, 12779, 17032,
               11569, 21907, 13888, 13997, 19573,
               17113, 9007, 16322, 18209, 12793, 24861,
               11561, 25424, 14493, 19255, 25907, 31106)
  idxOrder = c(15, 7, 23, 24, 26, 27, 32,
               2, 5, 8, 16, 35, 36,
               18, 19, 20, 21, 22, 34,37,
               1, 4, 11, 14, 17,
               3, 6, 9, 10, 12, 33,
               13, 25, 28, 29, 30, 31)
  clustNum[idxOrder] = clustNum
  clustNumUrb[idxOrder] = clustNumUrb
  clustNumRur[idxOrder] = clustNumRur
  
  # Initialize storage
  etaUrb.meas  = matrix(0, nrow = nSamp, ncol = 37)
  etaRur.meas  = matrix(0, nrow = nSamp, ncol = 37)
  eta.meas     = matrix(0, nrow = nSamp, ncol = 37)
  etaUrb.real  = matrix(0, nrow = nSamp, ncol = 37)
  etaRur.real  = matrix(0, nrow = nSamp, ncol = 37)
  eta.real     = matrix(0, nrow = nSamp, ncol = 37)
  etaUrb.overD = matrix(0, nrow = nSamp, ncol = 37)
  etaRur.overD = matrix(0, nrow = nSamp, ncol = 37)
  eta.overD    = matrix(0, nrow = nSamp, ncol = 37)
  
  # Get indicies in the sample
  spatIdx = res.inla$misc$configs$contents$start[2]
  intIdx  = res.inla$misc$configs$contents$start[4]
  urbIdx  = res.inla$misc$configs$contents$start[5]
  
  # Sample urban and rural
  nugSimStd = rnorm(1e5, mean = 0, sd = 1)
  for(i in 1:dim(etaUrb.meas)[1]){
    for(j in 1:37){
      # Nugget is measurement error
      etaUrb.meas[i,j] = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[urbIdx] + post.sample[[i]]$latent[spatIdx + j-1]
      etaRur.meas[i,j] = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spatIdx + j-1]
      
      # Nugget is real
      nugSimUrb = rnorm(clustNumUrb[j], mean = 0, sd = 1/sqrt(post.sample[[i]]$hyperpar[3]))
      nugSimRur = rnorm(clustNumRur[j], mean = 0, sd = 1/sqrt(post.sample[[i]]$hyperpar[3]))
      etaUrb.real[i,j] = logit(mean(expit(etaUrb.meas[i,j] + nugSimUrb)))
      etaRur.real[i,j] = logit(mean(expit(etaRur.meas[i,j] + nugSimRur)))
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[3])
      etaUrb.overD[i,j] = logit(mean(expit(etaUrb.meas[i,j] + nugSimStd*cSD)))
      etaRur.overD[i,j] = logit(mean(expit(etaRur.meas[i,j] + nugSimStd*cSD)))
    }
  }
  
  # Calculate combined samples
  for(j in 1:37){
    if(popList$prop.2018.child.admin1[j] < 1){
      eta.meas[,j]  = logit(popList$prop.2018.child.admin1[j]*expit(etaUrb.meas[,j])  + (1-popList$prop.2018.child.admin1[j])*expit(etaRur.meas[,j]))
      eta.real[,j]  = logit(popList$prop.2018.child.admin1[j]*expit(etaUrb.real[,j])  + (1-popList$prop.2018.child.admin1[j])*expit(etaRur.real[,j]))
      eta.overD[,j] = logit(popList$prop.2018.child.admin1[j]*expit(etaUrb.overD[,j]) + (1-popList$prop.2018.child.admin1[j])*expit(etaRur.overD[,j]))
    } else{
      eta.meas[,j]  = etaUrb.meas[,j]
      eta.real[,j]  = etaUrb.real[,j]
      eta.overD[,j] = etaUrb.overD[,j]
    }
  }
  
  # Calculate quantiles
  etaQuantUrb.meas = matrix(0, nrow = 37, ncol = 3)
  etaQuantRur.meas = matrix(0, nrow = 37, ncol = 3)
  etaQuant.meas = matrix(0, nrow = 37, ncol = 3)
  etaQuantRur.real = matrix(0, nrow = 37, ncol = 3)
  etaQuantUrb.real = matrix(0, nrow = 37, ncol = 3)
  etaQuant.real = matrix(0, nrow = 37, ncol = 3)
  etaQuantUrb.overD = matrix(0, nrow = 37, ncol = 3)
  etaQuantRur.overD = matrix(0, nrow = 37, ncol = 3)
  etaQuant.overD = matrix(0, nrow = 37, ncol = 3)
  for(i in 1:dim(etaUrb.meas)[2]){
    etaQuantUrb.meas[i,]  = quantile(etaUrb.meas[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuantRur.meas[i,]  = quantile(etaRur.meas[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuant.meas[i,]     = quantile(eta.meas[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuantRur.real[i,]  = quantile(etaRur.real[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuantUrb.real[i,]  = quantile(etaUrb.real[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuant.real[i,]     = quantile(eta.real[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuantUrb.overD[i,] = quantile(etaUrb.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuantRur.overD[i,] = quantile(etaRur.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuant.overD[i,]    = quantile(eta.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  pUrb.meas = 1/(1+exp(-etaQuantUrb.meas))
  pRur.meas = 1/(1+exp(-etaQuantRur.meas))
  p.meas = 1/(1+exp(-etaQuant.meas))
  pUrb.real = 1/(1+exp(-etaQuantUrb.real))
  pRur.real = 1/(1+exp(-etaQuantRur.real))
  p.real = 1/(1+exp(-etaQuant.real))
  pUrb.overD = 1/(1+exp(-etaQuantUrb.overD))
  pRur.overD = 1/(1+exp(-etaQuantUrb.overD))
  p.overD = 1/(1+exp(-etaQuant.overD))
  
  
  # Measurement error
  spat1.mError.ur = data.frame(admin1 = rep(levels(myData$admin1Fac), 2),
                               urb = rep( c("R", "U"), each = 37),
                               p_Low = c(pRur.meas[, 1], pUrb.meas[, 1]),
                               p_Med = c(pRur.meas[, 2], pUrb.meas[, 2]),
                               p_Upp = c(pRur.meas[, 3], pUrb.meas[, 3]),
                               numEAs = c(clustNumRur, clustNumUrb))
  idx = order(spat1.mError.ur$admin1, spat1.mError.ur$urb)
  spat1.mError.ur = spat1.mError.ur[idx,]
  
  spat1.mError = data.frame(admin1 = levels(myData$admin1Fac),
                            p_Low = p.meas[, 1],
                            p_Med = p.meas[, 2],
                            p_Upp = p.meas[, 3],
                            numEAs = clustNum)
  idx = order(spat1.mError$admin1)
  spat1.mError = spat1.mError[idx,]
  
  # Real
  spat1.real.ur = data.frame(admin1 = rep(levels(myData$admin1Fac), 2),
                             urb = rep( c("R", "U"), each = 37),
                             p_Low = c(pRur.real[, 1], pUrb.real[, 1]),
                             p_Med = c(pRur.real[, 2], pUrb.real[, 2]),
                             p_Upp = c(pRur.real[, 3], pUrb.real[, 3]),
                             numEAs = c(clustNumRur, clustNumUrb))
  idx = order(spat1.real.ur$admin1, spat1.real.ur$urb)
  spat1.real.ur = spat1.real.ur[idx,]
  
  spat1.real = data.frame(admin1 = levels(myData$admin1Fac),
                          p_Low = p.real[, 1],
                          p_Med = p.real[, 2],
                          p_Upp = p.real[, 3],
                          numEAs = clustNum)
  idx = order(spat1.real$admin1)
  spat1.real = spat1.real[idx,] 
  
  # OverD
  spat1.overD.ur = data.frame(admin1 = rep(levels(myData$admin1Fac), 2),
                              urb = rep( c("R", "U"), each = 37),
                              p_Low = c(pRur.overD[, 1], pUrb.overD[, 1]),
                              p_Med = c(pRur.overD[, 2], pUrb.overD[, 2]),
                              p_Upp = c(pRur.overD[, 3], pUrb.overD[, 3]),
                              numEAs = c(clustNumRur, clustNumUrb))
  idx = order(spat1.overD.ur$admin1, spat1.overD.ur$urb)
  spat1.overD.ur = spat1.overD.ur[idx,]
  
  spat1.overD = data.frame(admin1 = levels(myData$admin1Fac),
                           p_Low = p.overD[, 1],
                           p_Med = p.overD[, 2],
                           p_Upp = p.overD[, 3],
                           numEAs = clustNum)
  idx = order(spat1.overD$admin1)
  spat1.overD = spat1.overD[idx,]
  
  return(list(meas.ur = spat1.mError.ur,
              meas = spat1.mError,
              real.ur = spat1.real.ur,
              real = spat1.real,
              overD.ur = spat1.overD.ur,
              overD = spat1.overD,
              samples = list(p.overD = 1/(1+exp(-t(eta.overD))),
                             pRur.overD = 1/(1+exp(-t(etaRur.overD))),
                             pUrb.overD = 1/(1+exp(-t(etaUrb.overD))))))
}


# Aggregate on Admin2 level
aggBYM_admin2 = function(res.inla, popList, nameAdm1, myData, nSamp = 1000){
  # Draw posterior samples
  post.sample = inla.posterior.sample(n = nSamp, result = res.inla)
  
  # Initialize storage
  etaUrb.overD = matrix(0, nrow = nSamp, ncol = 774)
  etaRur.overD = matrix(0, nrow = nSamp, ncol = 774)
  eta.overD    = matrix(0, nrow = nSamp, ncol = 774)
  
  # Get indicies in the sample
  spatIdx = res.inla$misc$configs$contents$start[2]
  intIdx  = res.inla$misc$configs$contents$start[4]
  urbIdx  = res.inla$misc$configs$contents$start[5]
  
  # Sample urban and rural
  nugSimStd = rnorm(1e5, mean = 0, sd = 1)
  for(i in 1:dim(etaUrb.overD)[1]){
    for(j in 1:774){
      # Nugget is measurement error
      etaUrb.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[urbIdx] + post.sample[[i]]$latent[spatIdx + j-1]
      etaRur.tmp = post.sample[[i]]$latent[intIdx] + post.sample[[i]]$latent[spatIdx + j-1]
      
      # Nugget is overdispersion
      cSD = 1/sqrt(post.sample[[i]]$hyperpar[3])
      etaUrb.overD[i,j] = logit(mean(expit(etaUrb.tmp + nugSimStd*cSD)))
      etaRur.overD[i,j] = logit(mean(expit(etaRur.tmp + nugSimStd*cSD)))
    }
  }
  
  # Calculate combined samples
  for(j in 1:774){
    if(popList$prop.2018.child.admin2[j] < 1){
      eta.overD[,j] = logit(popList$prop.2018.child.admin2[j]*expit(etaUrb.overD[,j]) + (1-popList$prop.2018.child.admin2[j])*expit(etaRur.overD[,j]))
    } else{
      eta.overD[,j] = etaUrb.overD[,j]
    }
  }
  
  # Calculate quantiles
  etaQuantUrb.overD = matrix(0, nrow = 774, ncol = 3)
  etaQuantRur.overD = matrix(0, nrow = 774, ncol = 3)
  etaQuant.overD = matrix(0, nrow = 774, ncol = 3)
  for(i in 1:dim(etaUrb.overD)[2]){
    etaQuantUrb.overD[i,] = quantile(etaUrb.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuantRur.overD[i,] = quantile(etaRur.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    etaQuant.overD[i,]    = quantile(eta.overD[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  pUrb.overD = 1/(1+exp(-etaQuantUrb.overD))
  pRur.overD = 1/(1+exp(-etaQuantRur.overD))
  p.overD = 1/(1+exp(-etaQuant.overD))
  
  # OverD
  spat2.overD.ur = data.frame(admin2 = rep(levels(myData$admin2Fac), 2),
                              urb = rep( c("R", "U"), each = 774),
                              p_Low = c(pRur.overD[, 1], pUrb.overD[, 1]),
                              p_Med = c(pRur.overD[, 2], pUrb.overD[, 2]),
                              p_Upp = c(pRur.overD[, 3], pUrb.overD[, 3]))
  idx = order(spat2.overD.ur$admin2, spat2.overD.ur$urb)
  spat2.overD.ur = spat2.overD.ur[idx,]
  
  spat2.overD = data.frame(admin2 = levels(myData$admin2Fac),
                           p_Low = p.overD[, 1],
                           p_Med = p.overD[, 2],
                           p_Upp = p.overD[, 3])
  idx = order(spat2.overD$admin2)
  spat2.overD = spat2.overD[idx,]
  
  # Calculate admin1 samples
  admin1 = matrix(0, nrow = nSamp, ncol = 37)
  urb.admin1 = admin1
  rur.admin1 = admin1
  
  unNameAdm1 = sort(unique(nameAdm1))
  for(i in 1:37){
    idxAdm2 = which(nameAdm1 == unNameAdm1[i])
    totUrb = 0
    totRur = 0
    for(k in idxAdm2){
      idxUrb = popList$idxUrb[[k]]
      currUrb = sum(popList$popAdm2.2018[[k]][[2]][idxUrb], na.rm = TRUE)
      currRur = sum(popList$popAdm2.2018[[k]][[2]][!idxUrb], na.rm = TRUE)
      totUrb = totUrb + currUrb
      totRur = totRur + currRur
      
      admin1[,i]     = admin1[,i] + expit(eta.overD[,k])*(currUrb+currRur)
      if(currUrb > 0)
        urb.admin1[,i] = urb.admin1[,i] + expit(etaUrb.overD[,k])*currUrb
      if(currRur > 0)
        rur.admin1[,i] = rur.admin1[,i] + expit(etaRur.overD[,k])*currRur
    }
    admin1[,i] = admin1[,i]/(totUrb+totRur)
    urb.admin1[,i] = urb.admin1[,i]/totUrb
    rur.admin1[,i] = rur.admin1[,i]/totRur
  }
  
  # Calculate quantiles (admin1)
  pUrb.overD = matrix(0, nrow = 37, ncol = 3)
  pRur.overD = matrix(0, nrow = 37, ncol = 3)
  p.overD = matrix(0, nrow = 37, ncol = 3)
  for(i in 1:dim(admin1)[2]){
    pUrb.overD[i,] = quantile(urb.admin1[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    pRur.overD[i,] = quantile(rur.admin1[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
    p.overD[i,]    = quantile(admin1[,i], probs = c(0.025, 0.50, 0.975), na.rm = TRUE)
  }
  
  # OverD
  spat1.overD.ur = data.frame(admin1 = rep(levels(myData$admin1Fac), 2),
                              urb = rep( c("R", "U"), each = 37),
                              p_Low = c(pRur.overD[, 1], pUrb.overD[, 1]),
                              p_Med = c(pRur.overD[, 2], pUrb.overD[, 2]),
                              p_Upp = c(pRur.overD[, 3], pUrb.overD[, 3]))
  idx = order(spat1.overD.ur$admin1, spat1.overD.ur$urb)
  spat1.overD.ur = spat1.overD.ur[idx,]
  
  spat1.overD = data.frame(admin1 = levels(myData$admin1Fac),
                           p_Low = p.overD[, 1],
                           p_Med = p.overD[, 2],
                           p_Upp = p.overD[, 3])
  idx = order(spat1.overD$admin1)
  spat1.overD = spat1.overD[idx,]
  
  
  return(list(overD.ur = spat1.overD.ur,
              overD = spat1.overD,
              admin2.overD.ur = spat2.overD.ur,
              admin2.overD = spat2.overD,
              samples = list(p.overD = t(admin1),
                             pRur.overD = t(rur.admin1),
                             pUrb.overD = t(urb.admin1))))
}


# Extract population at risk
getPop = function(myData, sampleFrame = 2006){
  # Load geotiff files (full population)
  if(sampleFrame == 2006){
    pop2006 = raster("../Data/Nigeria_pop/nga_ppp_2006_UNadj.tif")
  }else{
    pop2006 = raster("../Data/Nigeria_pop/nga_ppp_2017_UNadj.tif")
  }
  pop2018 = raster("../Data/Nigeria_pop/nga_ppp_2018_UNadj.tif")
  
  # Load geotiff files for [1, 5) population by gender and year
  if(sampleFrame == 2006){
    pop2006.m1 = raster("../Data/Nigeria_pop/nga_m_1_2006.tif")
    pop2006.f1 = raster("../Data/Nigeria_pop/nga_f_1_2006.tif")
  } else{
    # I didn't include download 2017 children populations. We should take
    # this away. It is not used for anything.
    pop2006.m1 = raster("../Data/Nigeria_pop/nga_m_1_2006.tif")
    pop2006.f1 = raster("../Data/Nigeria_pop/nga_f_1_2006.tif")
  }
  pop2018.m1 = raster("../Data/Nigeria_pop/nga_m_1_2018.tif")
  pop2018.f1 = raster("../Data/Nigeria_pop/nga_f_1_2018.tif")
  
  # Extract cells inside each admin2 area
  popAdm2.2006 = list()
  for(i in 1:774){
    popAdm2.2006[[i]] = extract(x = pop2006, y = nigeriaMap[i,], cellnumbers = TRUE)
    tmpM1 = extract(x = pop2006.m1, y = nigeriaMap[i,], cellnumbers = TRUE)
    tmpM2 = extract(x = pop2006.f1, y = nigeriaMap[i,], cellnumbers = TRUE)
    popAdm2.2006[[i]][2] = list(childCount = tmpM1[[1]][,2] + tmpM2[[1]][,2])
  }
  popAdm2.2018 = list()
  for(i in 1:774){
    popAdm2.2018[[i]] = extract(x = pop2018, y = nigeriaMap[i,], cellnumbers = TRUE)
    tmpM1 = extract(x = pop2018.m1, y = nigeriaMap[i,], cellnumbers = TRUE)
    tmpM2 = extract(x = pop2018.f1, y = nigeriaMap[i,], cellnumbers = TRUE)
    popAdm2.2018[[i]][2] = list(childCount = tmpM1[[1]][,2] + tmpM2[[1]][,2])
  }
  
  # total census populations
  if(sampleFrame == 2006){
    totPopUrb = c(899703, 463094, 1110418, 1619155, 411089, 931288, 871623,
                  783977, 611908, 1387434, 539899, 355091, 492518,
                  452462, 2799079, 3925245, 1093024, 496745, 733481, 573709,
                  551090, 3387426, 1827862, 2282713, 1812748,
                  119472, 410562, 398369, 1920210, 1824233, 2412283,
                  1783711, 9112690, 1866997, 1608673, 2605526, 3969525)
    totPopRur = c(506682, 3790515, 2204910, 745731, 1459609, 3023387, 2335570,
                  2395523, 4039764, 2784113, 1825601, 1940734, 1828398,
                  3909329, 3315487, 5478986, 4709805, 2761395, 2968623, 2703142,
                  2293502, 790919, 348889, 986361, 2114782,
                  3782844, 1293325, 2495235, 2194864, 1409735, 2787892,
                  615371, 0, 1885127, 1852151, 810424, 1610771)
    idxOrder = c(15, 7, 23, 24, 26, 27, 32,
                 2, 5, 8, 16, 35, 36,
                 18, 19, 20, 21, 22, 34,37,
                 1, 4, 11, 14, 17,
                 3, 6, 9, 10, 12, 33,
                 13, 25, 28, 29, 30, 31)
    totPopUrb[idxOrder] = totPopUrb
    totPopRur[idxOrder] = totPopRur
    pop.census.2006 = data.frame(admin1 = levels(myData$admin1Fac), 
                                 urb = totPopUrb, 
                                 rur = totPopRur, 
                                 propUrb = totPopUrb/(totPopUrb+totPopRur))
  } else{
    tmpSampleFrame = read_excel("../Data/Nigeria_2017_proj/Nigeria_2017_admin1.xlsx")
    tmpSampleFrame = tmpSampleFrame[1:37,]
    pop.census.2006 = data.frame(admin1 = levels(myData$admin1Fac), 
                                 propUrb = as.vector(tmpSampleFrame$`Urban Pop%`))
  }
  
  # Find best thresholds in each admin1
  nameAdm1 = c()
  for(i in 1:774){
    nameAdm1 = c(nameAdm1, strsplit(nameVec[i], ":")[[1]][1])
  }
  unNameAdm1 = unique(nameAdm1)
  bestThr = rep(100, 37)
  bestProp = rep(0, 37)
  for(i in 1:37){
    low = 0
    high = 1e3
    idxAdm2 = which(nameAdm1 == unNameAdm1[i])
    for(j in 1:100){
      mid = (low+high)/2
      totUrb = 0
      totRur = 0
      for(k in idxAdm2){
        idxUrb = popAdm2.2006[[k]][[1]][,2] >= mid
        totUrb = totUrb + sum(popAdm2.2006[[k]][[1]][idxUrb,2], na.rm = TRUE)
        totRur = totRur + sum(popAdm2.2006[[k]][[1]][!idxUrb,2], na.rm = TRUE)
      }
      cProp = totUrb/(totUrb+totRur)
      if(cProp < pop.census.2006$propUrb[i]){
        high = mid
      } else{
        low = mid
      }
    }
    bestThr[i] = mid
    bestProp[i] = cProp
  }
  
  # Classify pixels
  idxUrb.admin2 = list()
  for(i in 1:37){
    idxAdm2 = which(nameAdm1 == unNameAdm1[i])
    for(k in idxAdm2){
      idxUrb = popAdm2.2006[[k]][[1]][,2] >= bestThr[i]
      idxUrb.admin2 = c(idxUrb.admin2, list(idxUrb = idxUrb))
    }
  }
  
  # Get proportions
  prop.2006.admin1 = rep(-1, 37)
  prop.2006.admin2 = rep(-1, 774)
  prop.2018.admin1 = rep(-1, 37)
  prop.2018.admin2 = rep(-1, 774)
  prop.2006.child.admin1 = rep(-1, 37)
  prop.2006.child.admin2 = rep(-1, 774)
  prop.2018.child.admin1 = rep(-1, 37)
  prop.2018.child.admin2 = rep(-1, 774)
  
  for(i in 1:37){
    idxAdm2 = which(nameAdm1 == unNameAdm1[i])
    totUrb.2006 = 0
    totRur.2006 = 0
    totUrb.2018 = 0
    totRur.2018 = 0
    totUrb.2006.child = 0
    totRur.2006.child = 0
    totUrb.2018.child = 0
    totRur.2018.child = 0
    for(k in idxAdm2){
      ## Full population
      currUrb = sum(popAdm2.2006[[k]][[1]][idxUrb.admin2[[k]],2], na.rm = TRUE)
      currRur = sum(popAdm2.2006[[k]][[1]][!idxUrb.admin2[[k]],2], na.rm = TRUE)
      prop.2006.admin2[k] = currUrb / (currUrb+currRur)
      totUrb.2006 = totUrb.2006 + currUrb
      totRur.2006 = totRur.2006 + currRur
      
      currUrb = sum(popAdm2.2018[[k]][[1]][idxUrb.admin2[[k]],2], na.rm = TRUE)
      currRur = sum(popAdm2.2018[[k]][[1]][!idxUrb.admin2[[k]],2], na.rm = TRUE)
      prop.2018.admin2[k] = currUrb / (currUrb+currRur)
      totUrb.2018 = totUrb.2018 + currUrb
      totRur.2018 = totRur.2018 + currRur
      
      ## Child population
      currUrb = sum(popAdm2.2006[[k]][[2]][idxUrb.admin2[[k]]], na.rm = TRUE)
      currRur = sum(popAdm2.2006[[k]][[2]][!idxUrb.admin2[[k]]], na.rm = TRUE)
      prop.2006.child.admin2[k] = currUrb / (currUrb+currRur)
      totUrb.2006.child = totUrb.2006.child + currUrb
      totRur.2006.child = totRur.2006.child + currRur
      
      currUrb = sum(popAdm2.2018[[k]][[2]][idxUrb.admin2[[k]]], na.rm = TRUE)
      currRur = sum(popAdm2.2018[[k]][[2]][!idxUrb.admin2[[k]]], na.rm = TRUE)
      prop.2018.child.admin2[k] = currUrb / (currUrb+currRur)
      totUrb.2018.child = totUrb.2018.child + currUrb
      totRur.2018.child = totRur.2018.child + currRur
    }
    prop.2006.admin1[i] = totUrb.2006/(totUrb.2006 + totRur.2006)
    prop.2018.admin1[i] = totUrb.2018/(totUrb.2018 + totRur.2018)
    prop.2006.child.admin1[i] = totUrb.2006.child/(totUrb.2006.child + totRur.2006.child)
    prop.2018.child.admin1[i] = totUrb.2018.child/(totUrb.2018.child + totRur.2018.child)
  }
  
  return(list(idxUrb = idxUrb.admin2,
              prop.2006.admin1 = prop.2006.admin1,
              prop.2018.admin1 = prop.2018.admin1,
              prop.2006.child.admin1 = prop.2006.child.admin1,
              prop.2018.child.admin1 = prop.2018.child.admin1,
              prop.2006.admin2 = prop.2006.admin2,
              prop.2018.admin2 = prop.2018.admin2,
              prop.2006.child.admin2 = prop.2006.child.admin2,
              prop.2018.child.admin2 = prop.2018.child.admin2,
              popAdm2.2018 = popAdm2.2018))
}

# Make areal plot with colours
plotAreaCol = function(fName, estVal, graph, colLim, leg = NULL){
  if(is.null(leg))
    leg = "MCV1"
  # Make png
  png(fName, width = 1700, height = 1200)
  
  # Set up data object for plotting
  nigeriaMapTmp = graph
  nigeriaMapTmp$MCV1 = estVal
  nigeria.df = merge(fortify(nigeriaMapTmp), as.data.frame(nigeriaMapTmp), by.x = "id", by.y = 0)
  nigeria.df$Longitude = nigeria.df$long
  nigeria.df$Latitude  = nigeria.df$lat
  
  # Plot
  map = ggplot() +
    geom_polygon(data = nigeria.df,
                 aes(x = Longitude, y = Latitude, group = group, fill = MCV1),
                 color = 'gray', size = .2)+
    scale_fill_viridis_c(direction = 1,
                         begin = 1,
                         end = 0,
                         limit = colLim,
                         name = leg) + 
    #scale_fill_continuous(trans = 'reverse',
    #                      type = "viridis",
    #                      limit = colLim) +
    coord_fixed() + 
    theme(text = element_text(size=50),
          legend.key.height = unit(4, 'cm'),
          legend.key.width  = unit(2, 'cm'))
  print(map)
  
  # Close graphical object
  dev.off()
}
