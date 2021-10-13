########################################################### 
##
## Processing the Urban/Rural fractions in Malawi
##
############################################################

library(INLA)
library(rgdal)
library(maptools)
library(raster)
library(spdep)
library(rgeos)
library(data.table)
library(broom)
library(ggplot2)
library(patchwork)
library(tidyr)
library(SUMMER)

##############################
# useful functions
##############################
match_loc_spdfvar <- function(loc_x, loc_y, spdf, varname){
    loc_df <- data.frame(x = loc_x, y = loc_y)
    coordinates(loc_df) <- ~x+y
    proj4string(loc_df) <- proj4string(spdf)
    return(over(loc_df, spdf)[, varname])
}

match_loc_raster <- function(loc_x, loc_y, raster){
    loc_df <- data.frame(x = loc_x, y = loc_y)
    return(raster::extract(x = raster, y = loc_df))
}

match_loc_raster_buffer <- function(loc_x, loc_y, raster, radius){
    loc_df <- data.frame(x = loc_x, y = loc_y)
    return(raster::extract(x = raster, y = loc_df, buffer = radius, fun = mean, na.rm = T))
}

# ##############################
# # get map points
# ##############################
data(MalawiMap)
geo <- MalawiMap

##############################
# get Census data
##############################
census <- read.csv("../data/Malawi_Pop_Frac/Table2-Census2008.csv")
census$Region <- gsub("        ", "", census$Region)
census$Region <- gsub(" City", "", census$Region)
census$Region[census$Region == "Mzuzu"] <- "Mzimba"
census <- aggregate(.~Region, data = census, FUN=sum)
# census$Region[census$Region %in% geo$NAME_1 == FALSE]
census$urban_prop <- census$Urban / census$Total
# 8.1% in Chiradzulu --> DHS 2015 says 8.3%
# 66.00% in Blantyre --> DHS 2015 says 65.5%
range(census$urban_prop)



##############################
#  download world pop raster
##############################
pop.raw.dir <-  "../Data/WorldPop-Population"
for(year in 2000:2019){
  for(age in c(0)){
    for(sex in c("f", "m")){
      file <- paste0(pop.raw.dir,'/mwi_', sex, '_', age, '_', year,'.tif')
      if(!file.exists(file)){
        url <- paste0("https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/", year, "/MWI/mwi_", sex, "_", age, "_", year, ".tif")
        download.file(url, file)
      }
    }
  }
}

##############################
# get world pop raster
##############################
year <- 2008
wp_pop_raw  <- raster(paste0("../Data/WorldPop-Population/mwi_ppp_", year, ".tif"))
wp_pop_raster <- raster::aggregate(wp_pop_raw, fact = 10, fun = sum, na.rm = T)
pred_dt <- as.data.frame(wp_pop_raster, xy=TRUE)
colnames(pred_dt)[3] <- "pop_all_wp"
DistrictName <- match_loc_spdfvar(loc_x = pred_dt$x,
                                        loc_y = pred_dt$y,
                                        spdf = geo,
                                        varname = c("ADM2_EN"))
pred_dt <- cbind(pred_dt, DistrictName)
pred_dt <- pred_dt[!is.na(pred_dt$DistrictName), ]
pred_dt[is.na(pred_dt$pop_all_wp), "pop_all_wp"] <- 0

get_partition <- function(pred_dt, census){
	regions <- unique(pred_dt$DistrictName)
	pred_dt$urban <- 0
	pred_dt$threshold <- NA
	for(r in regions){
		sub <- which(pred_dt$DistrictName == r)
		pop <- pred_dt$pop_all_wp[sub]
		pop_ordered <- sort(pop)
		thres <- pop_ordered[which(cumsum(pop_ordered) > sum(pop_ordered) * (1-census$urban_prop[census$Region == r]))[1]] 
		pred_dt$threshold[sub] <- thres
		pred_dt$urban[sub[pred_dt$pop_all_wp[sub] >= thres]] <- 1 
	}
	return(pred_dt)
}
pred_dt <- get_partition(pred_dt, census)
## Sanity checks
# aggregate(pop_all_wp ~ urban, data = pred_dt, FUN=sum) 
# sum(census$Urban)
# sum(census$Rural)

pred_dt$urban_factor <- "urban"
pred_dt$urban_factor[pred_dt$urban == 0] <- "rural"
# c(600, 300, 100, 10)
g1 <- ggplot(pred_dt, aes(x = x, y = y, fill = pop_all_wp)) + geom_raster() + scale_fill_viridis_c("population\ndensity", trans = 'log1p',  breaks = c(2e4, 5e3, 1000, 100, 10)) + theme_void()  
g2 <- ggplot(pred_dt, aes(x = x, y = y, fill = threshold)) + geom_raster() + theme_void()
g3 <- ggplot(pred_dt, aes(x = x, y = y, fill = urban_factor)) + geom_raster() + scale_fill_manual("", values = c("gray80", "black"))  + theme_void()
g <- g1 + g2 + g3
ggsave(g, file = "../Figures/urbanicity_mw_1km.png", width = 8, height = 5)

get_proportion <- function(year, pred_dt){

	wp_pop_raw1  <- raster(paste0("../Data/WorldPop-Population/mwi_f_0_", year, ".tif"))
	wp_pop_raw2  <- raster(paste0("../Data/WorldPop-Population/mwi_m_0_", year, ".tif"))

	wp_pop_raster1 <- raster::aggregate(wp_pop_raw1, fact = 10, fun = sum, na.rm = T)
	wp_pop_raster2 <- raster::aggregate(wp_pop_raw2, fact = 10, fun = sum, na.rm = T)

	pred_dt[, "pop_f"] <- match_loc_raster(loc_x = pred_dt$x,
                                            loc_y = pred_dt$y,
                                            raster = wp_pop_raster1)
	pred_dt$pop_f[is.na(pred_dt$pop_f)] <- 0

	pred_dt[, "pop_m"] <- match_loc_raster(loc_x = pred_dt$x,
                                            loc_y = pred_dt$y,
                                            raster = wp_pop_raster2)
	pred_dt$pop_m[is.na(pred_dt$pop_m)] <- 0
	pred_dt[, "pop_new"] <- pred_dt[, "pop_f"] + pred_dt[, "pop_m"]

	prop <- aggregate(pop_new ~ urban + DistrictName, data = pred_dt, FUN=sum)
	prop$urban <- ifelse(prop$urban == 1, "urban", "rural")
	prop <- spread(prop, urban, pop_new)
	prop$urban_prop <- prop$urban  / (prop$urban + prop$rural)
	prop$year <- year
	return(prop)
}

prop <- NULL
for(year in 2000:2019){
	tmp <- get_proportion(year, pred_dt)
	prop <- rbind(prop, tmp)
}
g4 <- ggplot(prop, aes(x = year, y = urban_prop, group = DistrictName)) + geom_line() + scale_x_continuous(breaks = seq(2000, 2019, by = 3)) + xlab("") + ylab("Proportion of population in urban clusters") + geom_vline(xintercept = 2008, color = "gray70", linetype = 2) + theme_bw()  
ggsave(g4, file = "../Figures/urbanicity_yearly_mw_under1_1km.png", width = 7, height = 5)

census$DistrictName <- census$Region
prop <- merge(prop, census[c("DistrictName", "Urban", "Rural")])
colnames(prop)[colnames(prop) == "Urban"] <- "urban_2008_census"
colnames(prop)[colnames(prop) == "Rural"] <- "rural_2008_census"

if(FALSE){
	#####################################################################
	# Table B1, FR319, Distribution of residential households DHS 2015
	#####################################################################
	household.weights <- data.frame(region = c("Chitipa","Karonga","Nkhata Bay","Rumphi","Mzimba","Likoma","Kasungu","Nkhotakota","Ntchisi","Dowa","Salima","Mchinji","Dedza","Ntcheu","Lilongwe","Mangochi","Machinga","Chiradzulu","Mwanza","Thyolo","Mulanje","Phalombe","Chikwawa","Nsanje","Balaka","Neno","Zomba","Blantyre"), 
		urban = 0.01 * c(7.74, 14.83, 5.38, 10.68, 18.29, 14.8, 7.04, 8.02, 3.28, 3.67, 7.85, 3.67, 3.08, 2.91, 35.84, 4.56, 4.61, 0.83, 15.65, 1.69, 2.55, 1.46, 2.89, 8.04, 6.66, 1.44, 11.79, 65.50))
	household.weights$rural <- 1 - household.weights$urban
	household.weights <- household.weights[weights$region != "Likoma", ]
	# colnames(household.weights) <- c("region", "1", "0")


	#####################################################################
	# Compare with the UR from census
	#####################################################################
	compare <- household.weights[, c(1, 2)]
	compare <- merge(compare, census[, c("Region", "urban_prop")], by.x = "region", by.y = "Region")
	compare$urban_to_rural_size = (1/compare$urban_prop - 1) / ((1-compare$urban)/compare$urban)
	summary(1/compare$urban_to_rural_size)
	# ggplot(compare) + aes(x = urban, y = urban_prop) + geom_point() + geom_abline(intercept = 0, slope = 1)
}
 
save(prop, census, file = "../Data/Malawi_Pop_Frac/pop_mw_under1.RData")