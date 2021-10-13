# Libraries
library(INLA)
library(rgdal)
library(SUMMER)
library(ggplot2)
library(dplyr)
library(maptools)
load("../Data-input/Malawi-Data.RData")
load("../Data-input/pop_mw_under1.RData")


#####################################################################
# Use the urban proportion from thresholding
#####################################################################
weights <- prop[, c("DistrictName", "year", "urban_prop")]
weights$rural_prop <- 1 - weights$urban_prop
colnames(weights) <- c("region", "years", "urban", "rural")

#####################################################################
## Process DHS data
#####################################################################
DHS2010.count <- subset(DHS2010.count, age == "0")
DHS2010.count$Y <- DHS2010.count$died
DHS2010.count$years <- DHS2010.count$time
DHS2010.count$region <- DHS2010.count$admin2
DHS2010.count$strata <- DHS2010.count$v025
DHS2010.count$cluster <- DHS2010.count$v001
DHS2010.count$survey <- 2

DHS2015.count <- subset(DHS2015.count, age == "0")
DHS2015.count$Y <- DHS2015.count$died
DHS2015.count$years <- DHS2015.count$time
DHS2015.count$region <- DHS2015.count$admin2
DHS2015.count$strata <- DHS2015.count$v025
DHS2015.count$cluster <- DHS2015.count$v001
DHS2015.count$survey <- 3


#####################################################################
## Fit three national estimates
#####################################################################
prop.nat <- prop %>% group_by(year) %>% 
					 summarize(urban = sum(urban), rural = sum(rural)) %>%
					 mutate(urban_prop = urban / (urban + rural))
weights.nat <- prop.nat[, c("year", "urban_prop")]
weights.nat$rural_prop <- 1 - weights.nat$urban_prop
colnames(weights.nat) <- c("years", "urban", "rural")

time.model.list <- c("ar1", "rw2", "rw1")
fit.bb.nat  <- proj.bb.nat  <- fit.bb.nat.strat <- proj.bb.nat.strat <- NULL
for(time.model in time.model.list){
	DHSboth <- rbind(DHS2010.count, DHS2015.count) 
	DHSboth$strata <- NA
	fit.bb.nat[[time.model]] <- smoothCluster(data = DHSboth, Amat = NULL, 
						family = "betabinomial", 
						strata.time.effect = TRUE,
						age.groups = c("0"), age.n = c(1), age.rw.group = c(1),
						year_label = 2000:2019, 
						time.model = time.model,
						survey.effect = TRUE, 
						common.trend = TRUE,
						verbose = FALSE)
	proj.bb.nat[[time.model]] <-  getSmoothed(fit.bb.nat[[time.model]], nsim = 1000, save.draws = TRUE) 
	 

	DHSboth <- rbind(DHS2010.count, DHS2015.count) 
	fit.bb.nat.strat[[time.model]] <- smoothCluster(data = DHSboth, Amat = NULL,
						family = "betabinomial", 
					    strata.time.effect = TRUE,
						age.groups = c("0"), age.n = c(1), age.rw.group = c(1),
						year_label = 2000:2019, 
						time.model = time.model, 
						common.trend = TRUE,
						survey.effect = TRUE)
	proj.bb.nat.strat[[time.model]] <-  getSmoothed(fit.bb.nat.strat[[time.model]], weight.strata = weights.nat, weight.frame = NULL, nsim = 1000, save.draws = TRUE) 

 }
save(fit.bb.nat, proj.bb.nat, fit.bb.nat.strat, proj.bb.nat.strat, file = "../Results/Malawi-Results-National-2000-2019-2surveys.RData")



#####################################################################
## Fit three subnational estimates
#####################################################################
time.model.list <- c("ar1", "rw2", "rw1")
fit.bb <- proj.bb <- fit.bb.strat <- proj.bb.strat <- NULL
for(time.model in time.model.list){
	DHSboth <- rbind(DHS2010.count, DHS2015.count)
	DHSboth$strata <- NA
	fit.bb[[time.model]] <- smoothCluster(data = DHSboth, Amat = malawiGraph, 
						family = "betabinomial", 
					    strata.time.effect = TRUE,
						age.groups = c("0"), age.n = c(1), age.rw.group = c(1),
						year_label = 2000:2019, 
						time.model = time.model, st.time.model = "ar1",
						pc.st.slope.u = 1, 
						pc.st.slope.alpha = 0.01,
						survey.effect = TRUE,  
						common.trend = TRUE,
						verbose = FALSE)
	proj.bb[[time.model]] <-  getSmoothed(fit.bb[[time.model]], nsim = 1000, save.draws = TRUE) 


	DHSboth <- rbind(DHS2010.count, DHS2015.count)
	fit.bb.strat[[time.model]] <- smoothCluster(data = DHSboth, Amat = malawiGraph, 
						family = "betabinomial", 
					    strata.time.effect = TRUE,
						age.groups = c("0"), age.n = c(1), age.rw.group = c(1),
						year_label = 2000:2019, 
						time.model = time.model, st.time.model = "ar1",
						pc.st.slope.u = 1, 
						pc.st.slope.alpha = 0.01,
						survey.effect = TRUE, 
						common.trend = TRUE)
	# check this work with multiple frames
	proj.bb.strat[[time.model]] <-  getSmoothed(fit.bb.strat[[time.model]], weight.strata = weights, weight.frame = NULL, nsim = 1000, save.draws = TRUE) 
}
save(fit.bb, proj.bb, fit.bb.strat, proj.bb.strat, file = "../Results/Malawi-Results-2000-2019-2surveys.RData")


#####################################################################
## Get Direct estimates
#####################################################################
direct <- getDirectList(births = list(`DHS 2010`=subset(DHS2010, age == "0"),
									 `DHS 2015`=subset(DHS2015, age == "0")),
						years = as.character(2000:2019), 
						regionVar = "admin2", timeVar = "time",
						clusterVar = "~v001", ageVar = "age", 
						weightsVar = "v005")
# Combine surveys
direct.comb <- aggregateSurvey(direct)
save(direct, direct.comb, file = "../Results/Malawi-Results-Direct.RData")


#####################################################################
## Get Direct estimates 5 years
#####################################################################
DHS2010_4yr <- DHS2010
DHS2010_4yr$time <- recode(DHS2010_4yr$time, "2000"="00-03", "2001"="00-03", 
											 "2002"="00-03", "2003"="00-03", 
											 "2004"="04-07", "2005"="04-07", 
											 "2006"="04-07", "2007"="04-07", 
											 "2008"="08-11", "2009"="08-11", 
											 "2010"="08-11", "2011"="08-11", 
											 "2012"="12-15", "2013"="12-15", 
											 "2014"="12-15", "2015"="12-15")
DHS2015_4yr <- DHS2015
DHS2015_4yr$time <- recode(DHS2015_4yr$time, "2000"="00-03", "2001"="00-03", 
											 "2002"="00-03", "2003"="00-03", 
											 "2004"="04-07", "2005"="04-07", 
											 "2006"="04-07", "2007"="04-07", 
											 "2008"="08-11", "2009"="08-11", 
											 "2010"="08-11", "2011"="08-11", 
											 "2012"="12-15", "2013"="12-15", 
											 "2014"="12-15", "2015"="12-15")


direct.4yr <- getDirectList(births = list(`DHS 2010`=subset(DHS2010_4yr, age == "0"),
									 `DHS 2015`=subset(DHS2015_4yr, age == "0")),
						years = c("00-03", "04-07", "08-11", "12-15"), 
						regionVar = "admin2", timeVar = "time",
						clusterVar = "~v001", ageVar = "age", 
						weightsVar = "v005")
# Combine surveys
direct.comb.4yr <- aggregateSurvey(direct.4yr)
save(direct.4yr, direct.comb.4yr, file = "../Results/Malawi-Results-Direct-4yr.RData")

#####################################################################
## Get Smooth Direct estimates
#####################################################################
fit.smooth.direct <- smoothDirect(data = direct.comb, Amat = malawiGraph, year_label = as.character(2000:2019), year_range = c(2000, 2019), time.model = 'ar1', st.time.model = "ar1", m = 1, type.st = 4, pc.alpha = 0.05, pc.u = 1)
smooth.direct <- getSmoothed(fit.smooth.direct)
# plot(smooth.direct, plot.CI = TRUE) + facet_wrap(~region)
save(fit.smooth.direct, smooth.direct, file = "../Results/Malawi-Results-Smooth-Direct.RData")
