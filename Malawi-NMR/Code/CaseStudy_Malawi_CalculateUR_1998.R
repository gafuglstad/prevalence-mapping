########################################################### 
##
## Processing the Urban/Rural fractions in Malawi 
##   using the 1998 census
##   data from http://www.nsomalawi.mw/index.php?option=com_content&view=article&id=127%3A1998-population-and-housing-census&catid=8&Itemid=3
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
census <- read.csv("../Data/Malawi_Pop_Frac/Table1-Census1998.csv")
census$urban_prop <- census$Urban / census$Total
range(census$urban_prop)

##############################
# map district to regions
##############################
admin1_list <- data.frame(admin2 = geo$ADM2_EN, admin1 = geo$ADM1_EN)



##############################
# get world pop raster
##############################
year <- 2000
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
pred_dt$RegionName <- admin1_list$admin1[match(pred_dt$DistrictName, admin1_list$admin2)]

get_partition <- function(pred_dt, census){
	regions <- unique(pred_dt$RegionName)
	pred_dt$urban <- 0
	pred_dt$threshold <- NA
	for(r in regions){
		sub <- which(pred_dt$RegionName == r)
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
ggsave(g, file = "../Figures/urbanicity_mw_1km_1998.png", width = 8, height = 5)

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
	prop$urban[is.na(prop$urban)] <- 0
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
ggsave(g4, file = "../Figures/urbanicity_yearly_mw_under1_1km_1998.png", width = 7, height = 5)
prop.1998 <- prop
census.1998 <- census
save(prop.1998, census, file = "../Data/Malawi_Pop_Frac/pop_mw_under1_1998.RData")

