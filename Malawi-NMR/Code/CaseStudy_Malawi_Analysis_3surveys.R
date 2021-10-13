# Libraries
library(INLA)
library(rgdal)
library(SUMMER)
library(ggplot2)
library(dplyr)
library(maptools)
load("../Data-input/Malawi-Data.RData")
load("../Data-input/pop_mw_under1.RData")
load("../Data-input/pop_mw_under1_1998.RData")


#####################################################################
# Use the urban proportion from thresholding
#####################################################################
weights0 <- prop.1998[, c("DistrictName", "year", "urban_prop")]
weights0$rural_prop <- 1 - weights0$urban_prop
colnames(weights0) <- c("region", "years", "urban", "rural")
weights0$frame <- 1

weights <- prop[, c("DistrictName", "year", "urban_prop")]
weights$rural_prop <- 1 - weights$urban_prop
colnames(weights) <- c("region", "years", "urban", "rural")
weights$frame <- 2

weights <- rbind(weights0, weights)
#####################################################################
## Process DHS data
#####################################################################
DHS2004.count <- subset(DHS2004.count, age == "0")
DHS2004.count$Y <- DHS2004.count$died
DHS2004.count$years <- DHS2004.count$time
DHS2004.count$region <- DHS2004.count$admin2
DHS2004.count$strata <- DHS2004.count$v025
DHS2004.count$cluster <- DHS2004.count$v001
DHS2004.count$survey <- 1

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


# #####################################################################
# ## Fit three subnational estimates
# #####################################################################
# time.model.list <- c("ar1", "rw2", "rw1")
# fit.bb <- proj.bb <- fit.bb.strat <- proj.bb.strat <- NULL
# for(time.model in time.model.list){
# 	DHSboth <- cbind(DHS2004.count, frame = 1)
# 	DHSboth <- rbind(DHSboth, cbind(rbind(DHS2010.count, DHS2015.count), frame = 2))
# 	DHSboth$strata <- NA
# 	fit.bb[[time.model]] <- smoothCluster(data = DHSboth, Amat = malawiGraph, 
# 						family = "betabinomial", 
# 					    strata.time.effect = TRUE,
# 						age.groups = c("0"), age.n = c(1), age.rw.group = c(1),
# 						year_label = 2000:2019, 
# 						time.model = time.model, st.time.model = "ar1",
# 						pc.st.slope.u = 1, 
# 						pc.st.slope.alpha = 0.01,
# 						survey.effect = TRUE, 
# 						verbose = FALSE, 
# 						common.trend = (time.model == "ar1"))
# 	proj.bb[[time.model]] <-  getSmoothed(fit.bb[[time.model]], nsim = 1000, save.draws = TRUE) 
# 	toplot <- proj.bb[[time.model]]$stratified
# 	toplot$frame <- c("Partition defined by the 1998 Census", "Partition defined by the 2008 Census")[toplot$strata]
# 	g <- ggplot(toplot) + 
# 			  aes(x = years.num, y = median, ymin = lower, ymax = upper, color = frame, fill = frame, group = frame) + 
# 			  geom_ribbon(alpha = 0.2, color = NA) + 
# 			  geom_line() + 
# 			  geom_point() +
# 			  facet_wrap(~region, ncol = 7) +
# 			  geom_vline(xintercept = 2004, linetype=3, color = "gray10", size=.2) +
# 			  geom_vline(xintercept = 2010, linetype=3, color = "gray10", size=.2) +
# 			  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=.2) +
# 			  theme_bw()  + theme(legend.position = "bottom") +
# 			  scale_color_brewer(NULL, palette = "Set2", guide = guide_legend(ncol = 1)) +
# 			  scale_fill_brewer(NULL, palette = "Set2", guide = FALSE) + xlab("") + ylab("NMR")
# 	ggsave(g, file = paste0("../Figures/Subnational-frame-comparison-agg-no-strata-", time.model, ".pdf"), width = 25, height = 15)

# 	DHSboth <- cbind(DHS2004.count, frame = 1)
# 	DHSboth <- rbind(DHSboth, cbind(rbind(DHS2010.count, DHS2015.count), frame = 2))
# 	fit.bb.strat[[time.model]] <- smoothCluster(data = DHSboth, Amat = malawiGraph, 
# 						family = "betabinomial", 
# 					    strata.time.effect = TRUE,
# 						age.groups = c("0"), age.n = c(1), age.rw.group = c(1),
# 						year_label = 2000:2019, 
# 						time.model = time.model, st.time.model = "ar1",
# 						pc.st.slope.u = 1, 
# 						pc.st.slope.alpha = 0.01,
# 						survey.effect = TRUE, 
# 						common.trend = (time.model == "ar1"))
# 	# check this work with multiple frames
# 	proj.bb.strat[[time.model]] <-  getSmoothed(fit.bb.strat[[time.model]], weight.strata = weights, weight.frame = NULL, nsim = 1000, save.draws = TRUE) 
# 	toplot <- proj.bb.strat[[time.model]]$overall
# 	toplot$frame <- c("Partition defined by the 1998 Census", "Partition defined by the 2008 Census")[toplot$frame]
# 	g <- ggplot(toplot) + 
# 			  aes(x = years.num, y = median, ymin = lower, ymax = upper, color = frame, fill = frame, group = frame) + 
# 			  geom_ribbon(alpha = 0.2, color = NA) + 
# 			  geom_line() + 
# 			  geom_point() + 
#   			  facet_wrap(~region, ncol = 7) +
# 			  geom_vline(xintercept = 2004, linetype=3, color = "gray10", size=.2) +
# 			  geom_vline(xintercept = 2010, linetype=3, color = "gray10", size=.2) +
# 			  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=.2) +
# 			  theme_bw()  + theme(legend.position = "bottom") +
# 			  scale_color_brewer(NULL, palette = "Set2", guide = guide_legend(ncol = 1)) +
# 			  scale_fill_brewer(NULL, palette = "Set2", guide = FALSE) + xlab("") + ylab("NMR")
# 	ggsave(g, file = paste0("../Figures/Subnational-frame-comparison-agg-", time.model, ".pdf"), width = 25, height = 15)
# }
# save(fit.bb, proj.bb, fit.bb.strat, proj.bb.strat, file = "../Results/Malawi-Results.RData")


#####################################################################
## Fit three national estimates
#####################################################################
prop.nat.1998 <- prop.1998 %>% group_by(year) %>% 
					 summarize(urban = sum(urban), rural = sum(rural)) %>%
					 mutate(urban_prop = urban / (urban + rural))
prop.nat <- prop %>% group_by(year) %>% 
					 summarize(urban = sum(urban), rural = sum(rural)) %>%
					 mutate(urban_prop = urban / (urban + rural))

weights.nat.1998 <- prop.nat.1998[, c("year", "urban_prop")]
weights.nat.1998$rural_prop <- 1 - weights.nat.1998$urban_prop
colnames(weights.nat.1998) <- c("years", "urban", "rural")
weights.nat.1998$frame <- 1

weights.nat <- prop.nat[, c("year", "urban_prop")]
weights.nat$rural_prop <- 1 - weights.nat$urban_prop
colnames(weights.nat) <- c("years", "urban", "rural")
weights.nat$frame <- 2

weights.nat <- data.frame(rbind(weights.nat.1998, weights.nat))
if(FALSE){
	ggplot(weights.nat, aes(x = years, y = urban, group = frame, color = frame)) + geom_line()
}

time.model.list <- c("ar1", "rw2", "rw1")
fit.bb.nat  <- proj.bb.nat  <- fit.bb.nat.strat <- proj.bb.nat.strat <- NULL
for(time.model in time.model.list){
	# DHSboth <- cbind(DHS2004.count, frame = 1)
	# DHSboth <- rbind(DHSboth, cbind(rbind(DHS2010.count, DHS2015.count), frame = 2))
	# DHSboth$strata <- NA
	# fit.bb.nat[[time.model]] <- smoothCluster(data = DHSboth, Amat = NULL, 
	# 					family = "betabinomial", 
	# 					strata.time.effect = TRUE,
	# 					age.groups = c("0"), age.n = c(1), age.rw.group = c(1),
	# 					year_label = 2000:2019, 
	# 					time.model = time.model,
	# 					survey.effect = TRUE, 
	# 					verbose = FALSE, 
	# 					common.trend = (time.model == "ar1"))
	# proj.bb.nat[[time.model]] <-  getSmoothed(fit.bb.nat[[time.model]], nsim = 1000, save.draws = TRUE) 
	# # strata is renamed to be the same as frame. Need to make aggregation codes work in this situation in getSmoothed.
	# toplot <- proj.bb.nat[[time.model]]$stratified
	# toplot$frame <- c("Partition defined by the 1998 Census", "Partition defined by the 2008 Census")[toplot$strata]
	# g <- ggplot(toplot) + 
	# 		  aes(x = years.num, y = median, ymin = lower, ymax = upper, color = frame, fill = frame, group = frame) + 
	# 		  geom_ribbon(alpha = 0.2, color = NA) + 
	# 		  geom_line() + 
	# 		  geom_point() + 
	# 		  geom_vline(xintercept = 2004, linetype=3, color = "gray10", size=.2) +
	# 		  geom_vline(xintercept = 2010, linetype=3, color = "gray10", size=.2) +
	# 		  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=.2) +
	# 		  theme_bw()  + theme(legend.position = "bottom") +
	# 		  scale_color_brewer(NULL, palette = "Set2", guide = guide_legend(ncol = 1)) +
	# 		  scale_fill_brewer(NULL, palette = "Set2", guide = FALSE) + xlab("") + ylab("NMR")
	# ggsave(g, file = paste0("../Figures/National-frame-comparison-agg-no-strata-", time.model, ".pdf"), width = 5, height = 5)
	

	DHSboth <- cbind(DHS2004.count, frame = 1)
	DHSboth <- rbind(DHSboth, cbind(rbind(DHS2010.count, DHS2015.count), frame = 2))
	fit.bb.nat.strat[[time.model]] <- smoothCluster(data = DHSboth, Amat = NULL,
						family = "betabinomial", 
					    strata.time.effect = TRUE,
						age.groups = c("0"), age.n = c(1), age.rw.group = c(1),
						year_label = 2000:2019, 
						time.model = time.model,
						survey.effect = TRUE, 
						common.trend = (time.model == "ar1"))
	proj.bb.nat.strat[[time.model]] <-  getSmoothed(fit.bb.nat.strat[[time.model]], weight.strata = weights.nat, weight.frame = NULL, nsim = 1000, save.draws = TRUE) 
	# plot(proj.bb.nat.strat[[time.model]]$overall) + aes(color = frame)
	toplot <- proj.bb.nat.strat[[time.model]]$stratified
	toplot$frame <- "Partition defined by the 1998 Census"
	toplot$frame[toplot$strata %in% c("2-urban", "2-rural")] <- "Partition defined by the 2008 Census"
	toplot$urban <- c("urban")
	toplot$urban[toplot$strata %in% c("1-rural", "2-rural")] <- "rural"
	g <- ggplot(toplot) + 
			  aes(x = years.num, y = median, ymin = lower, ymax = upper, color = urban, fill = urban, group = urban) + 
			  geom_ribbon(alpha = 0.2, color = NA) + 
			  geom_line() + 
			  geom_point() + 
			  geom_vline(xintercept = 2004, linetype=3, color = "gray10", size=.2) +
			  geom_vline(xintercept = 2010, linetype=3, color = "gray10", size=.2) +
			  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=.2) +
			  theme_bw() + facet_wrap(~frame, ncol = 4) + 
			  scale_color_brewer(NULL, palette = "Set1") +
			  scale_fill_brewer(NULL, palette = "Set1", guide = FALSE) + 
			  theme(legend.position = "bottom") + xlab("") + ylab("NMR")
	ggsave(g, file = paste0("../Figures/National-frame-comparison-", time.model, ".pdf"), width = 8, height = 4)
	
	# toplot <- proj.bb.nat.strat[[time.model]]$overall
	# toplot$frame <- c("Partition defined by the 1998 Census", "Partition defined by the 2008 Census")[toplot$frame]
	# g <- ggplot(toplot) + 
	# 		  aes(x = years.num, y = median, ymin = lower, ymax = upper, color = frame, fill = frame, group = frame) + 
	# 		  geom_ribbon(alpha = 0.2, color = NA) + 
	# 		  geom_line() + 
	# 		  geom_point() + 
	# 		  geom_vline(xintercept = 2004, linetype=3, color = "gray10", size=.2) +
	# 		  geom_vline(xintercept = 2010, linetype=3, color = "gray10", size=.2) +
	# 		  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=.2) +
	# 		  theme_bw()  + theme(legend.position = "bottom") +
	# 		  scale_color_brewer(NULL, palette = "Set2", guide = guide_legend(ncol = 1)) +
	# 		  scale_fill_brewer(NULL, palette = "Set2", guide = FALSE) + xlab("") + ylab("NMR")
	# ggsave(g, file = paste0("../Figures/National-frame-comparison-agg-", time.model, ".pdf"), width = 5, height = 5)
 }
save(fit.bb.nat, proj.bb.nat, fit.bb.nat.strat, proj.bb.nat.strat, file = "../Results/Malawi-Results-National-2000-2019_3surveys.RData")



# #####################################################################
# ## Get Direct estimates
# #####################################################################
# direct <- getDirectList(births = list(`DHS 2010`=subset(DHS2010, age == "0"),
# 									 `DHS 2015`=subset(DHS2015, age == "0")),
# 						years = as.character(2000:2019), 
# 						regionVar = "admin2", timeVar = "time",
# 						clusterVar = "~v001", ageVar = "age", 
# 						weightsVar = "v005")
# # Combine surveys
# direct.comb <- aggregateSurvey(direct)
# save(direct, direct.comb, file = "../Results/Malawi-Results-Direct.RData")
