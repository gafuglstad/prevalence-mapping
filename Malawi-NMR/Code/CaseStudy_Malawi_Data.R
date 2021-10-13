########################################################### 
##
## Processing input files for smoothing NMR
## 
## 	this script contains codes to 
##	 	1. download Malawi 2010, 2015 DHS 
##		2. process files for model fitting
##		3. process maps
##		4. process urban/rural proportions
##
############################################################


# load libraries
library(SUMMER)
library(ggplot2)
library(patchwork)
library(rdhs)
library(stringr)
library(dplyr)
library(rgdal)
library(maptools)
library(geosphere)
# download data files from the DHS website if not exist in local folder
# Registration with DHS is required to download the files directly from R
# See `?rdhs::set_rdhs_config` for details
if(!file.exists("../Data/DHS_downloads.RData")){
	sv <- dhs_surveys(countryIds = "MW", surveyType = "DHS", 
					  surveyYearStart = 2000)
	BR <- dhs_datasets(surveyIds = sv$SurveyId, fileFormat = "Flat", 
					   fileType = "BR")
	BRfiles <- get_datasets(BR$FileName, reformat=TRUE)
	GPS <- dhs_datasets(surveyIds = sv$SurveyId, fileFormat = "Flat", 
						fileType = "GE")
	GPSfiles <- get_datasets(GPS$FileName, reformat=TRUE)
	Surv2004 <- readRDS(BRfiles[[2]])
	Surv2010 <- readRDS(BRfiles[[3]])
	Surv2015 <- readRDS(BRfiles[[4]])
	DHS2004.geo <- readRDS(GPSfiles[[2]])
	DHS2010.geo <- readRDS(GPSfiles[[3]])
	DHS2015.geo <- readRDS(GPSfiles[[4]])
	save(Surv2004, Surv2010, Surv2015, DHS2004.geo, DHS2010.geo, DHS2015.geo, file = "../Data/DHS_downloads.RData")
}else{
	load("../Data/DHS_downloads.RData")
}


# Process data into birth records
# Double check strata for 2004 DHS, v023 is 11 large regions, or should we use v024 (3 regions)? v022 does not seem correct in this survey.
DHS2004 <- getBirths(data = Surv2004,
				  month.cut = c(1, 12, 24, 36, 48, 60),
                  year.cut = seq(2000, 2020, by = 1), strata = c("v023", "v025"))
DHS2010 <- getBirths(data = Surv2010,
				  month.cut = c(1, 12, 24, 36, 48, 60),
                  year.cut = seq(2000, 2020, by = 1), strata = "v022")
DHS2015 <- getBirths(data = Surv2015,
				  month.cut = c(1, 12, 24, 36, 48, 60),
                  year.cut = seq(2000, 2020, by = 1), strata = "v022")
# Subset the observations in partial year 2016 
DHS2015 <- subset(DHS2015, time != 2016)
 
# get the list of clusters for 2010 DHS
cluster.list <- data.frame(DHS2010.geo) %>% 
		 distinct(DHSCLUST, DHSREGNA, URBAN_RURA) %>%
		 mutate(admin2 = gsub(" - rural", "", DHSREGNA)) %>% 
		 mutate(admin2 = gsub(" - urban", "", admin2)) %>%
		 mutate(admin2 = recode(admin2, "nkhatabay" = "nkhata bay")) %>%
		 mutate(admin2 = recode(admin2, "nkhota kota" = "nkhotakota")) %>%
		 mutate(admin2 =  str_to_title(admin2)) %>%
		 mutate(urban = ifelse(URBAN_RURA=="U", 1, 0))  %>%
		 select(v001 = DHSCLUST, admin2, urban) 
DHS2010 <- DHS2010 %>% left_join(cluster.list)
vars <- c("v001", "v025", "admin2", "time", "age")
DHS2010.count <- getCounts(DHS2010[, c(vars, "died")], variables = 'died', 
				  by = vars, drop=TRUE)

# get the list of clusters for 2015 survey
cluster.list <- data.frame(DHS2015.geo) %>% 
		 mutate(v001 = DHSCLUST) %>%
		 left_join(DHS2015, by = "v001") %>% 
		 distinct(v001, v022, URBAN_RURA) %>%
		 mutate(admin2 = gsub(" - rural", "", v022)) %>% 
		 mutate(admin2 = gsub(" - urban", "", admin2)) %>%
		 mutate(admin2 = recode(admin2, "nkhatabay" = "nkhata bay")) %>%
		 mutate(admin2 = recode(admin2, "nkhota kota" = "nkhotakota")) %>%
		 mutate(admin2 =  str_to_title(admin2)) %>%
		 mutate(urban = ifelse(URBAN_RURA=="U", 1, 0))  %>%
		 select(v001, admin2, urban) 
DHS2015 <- DHS2015 %>% left_join(cluster.list)
DHS2015.count <- getCounts(DHS2015[, c(vars, "died")], variables = 'died', 
				  by = vars, drop=TRUE)
# get the map
shpFile = "../Data/shapefiles"
malawiMap = readOGR(shpFile, layer = "mwi_admbnda_adm2_nso_20181016", verbose = FALSE)

IDs = match(malawiMap$ADM2_EN, malawiMap$ADM2_EN[malawiMap$ADM2_EN %in% DHS2015$admin2]) - 1
IDs[malawiMap$ADM2_EN == "Blantyre City"] = IDs[malawiMap$ADM2_EN == "Blantyre"]
IDs[malawiMap$ADM2_EN == "Lilongwe City"] = IDs[malawiMap$ADM2_EN == "Lilongwe"]
IDs[malawiMap$ADM2_EN == "Mzuzu City"] = IDs[malawiMap$ADM2_EN == "Mzimba"]
IDs[malawiMap$ADM2_EN == "Zomba City"] = IDs[malawiMap$ADM2_EN == "Zomba"]
malawiJoin = unionSpatialPolygons(SpP = malawiMap,
                                    IDs = IDs)
malawiMap = malawiMap[malawiMap$ADM2_EN %in% DHS2015$admin2, ]
row.names(malawiMap) = as.character(0:27)
malawiMap = SpatialPolygonsDataFrame(malawiJoin, malawiMap@data)
malawiMap = malawiMap[order(malawiMap$ADM2_EN),]
malawiGraph <- getAmat(malawiMap, malawiMap$ADM2_EN)
# mapPlot(data=NULL, geo = malawiMap, by.geo = "ADM2_EN") + theme_dark()


## Now add 2004 data by mapping GPS to polygon

cluster.list <- data.frame(DHS2004.geo) %>% 
		 mutate(v001 = DHSCLUST) %>%
		 left_join(DHS2004, by = "v001") %>% 
		 distinct(v001, URBAN_RURA, LATNUM,  LONGNUM) %>%
		 mutate(urban = ifelse(URBAN_RURA=="U", 1, 0)) 
cluster.list <- cluster.list[cluster.list$LONGNUM > 0, ] 
cluster.list2 <- mapPoints(cluster.list, malawiMap, long = "LONGNUM", lat = "LATNUM", names = "ADM2_EN")
miss <- cluster.list2[is.na(cluster.list2$ADM2_EN), ]
polygons <- SpatialPolygons(malawiMap@polygons)
proj4string(polygons) <-  proj4string(malawiMap)
miss.nearest <- dist2Line(miss[,c("LONGNUM", "LATNUM")], polygons)
miss.nearest <- data.frame(miss.nearest)
miss.nearest$ADM2_EN <- malawiMap@data$ADM2_EN[as.numeric(miss.nearest$ID)]
for(i in 1:dim(miss)[1]){
	cluster.list2[cluster.list2$v001 == miss$v001[i], "ADM2_EN"] <- miss.nearest$ADM2_EN[i]
}
if(FALSE){
	plot(malawiMap)
	check <- subset(cluster.list2, v001 %in% miss$v001)
	text(check$LONGNUM, check$LATNUM, check$ADM2_EN, col = 2, cex = .5)
} 

DHS2004 <- DHS2004 %>% left_join(cluster.list2, by = "v001")%>% mutate(admin2 = ADM2_EN)
DHS2004 <- subset(DHS2004, !is.na(ADM2_EN))
DHS2004 <- subset(DHS2004, time != 2005)
vars <- c("v001", "v025", "admin2", "time", "age")
DHS2004.count <- getCounts(DHS2004[, c(vars, "died")], variables = 'died', 
				  by = vars, drop=TRUE)

# Save the data for further modeling
save(DHS2004, DHS2004.count, DHS2010, DHS2015, DHS2010.count, DHS2015.count, malawiMap, malawiGraph, file = "../Data/Malawi_FULL_DHS/Malawi-Full-Data.RData")

# Save the aggregated data for reproducibility
save(DHS2010.count, DHS2015.count, malawiMap, malawiGraph, file = "../Data/Malawi_AGG_DHS/Malawi-Agg-Data.RData")
