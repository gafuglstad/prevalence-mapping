# Libraries
library(INLA)
library(rgdal)
library(SUMMER)
library(ggplot2)
library(dplyr)
library(maptools)
load("../Data/Malawi_FULL_DHS/Malawi-Full-Data.RData")


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
