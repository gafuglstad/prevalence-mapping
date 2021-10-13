library(SUMMER)
library(ggplot2)
library(patchwork)
library(dplyr)

#############################
## National comparison plots
#############################

load("../Results/Malawi-Results-Direct.RData")
load("../Results/Malawi-Results-Direct-4yr.RData")
load("../Results/Malawi-Results-Smooth-Direct.RData")
load("../Results/Malawi-Results-National-2000-2019-2surveys.RData")
load("../Results/Malawi-Results-2000-2019-2surveys.RData")
load("../Data/Malawi_Pop_Frac/pop_mw_under1.RData")

igme <- read.csv("../Data/IGME/UN_IGME_Malawi.csv")
UN <- data.frame(year = as.numeric(gsub("-06", "", igme$TIME_PERIOD)), 
				 est.un = igme$OBS_VALUE / 1000, 
				 lower.un = igme$LOWER_BOUND / 1000, 
				 upper.un = igme$UPPER_BOUND / 1000 )


#####################################################################
# Use the urban proportion from thresholding
#####################################################################
weights <- prop[, c("DistrictName", "year", "urban_prop")]
weights$rural_prop <- 1 - weights$urban_prop
colnames(weights) <- c("region", "years", "urban", "rural")


dir <- "../Figures/"

ymax <- 0
for(i in 1:length(proj.bb.nat.strat)){
	ymax <- max(c(ymax, proj.bb.nat.strat[[i]]$stratified$upper), na.rm = TRUE)
}

g <- NULL
for(i in 1:length(proj.bb.nat.strat)){
	tmp <- proj.bb.nat.strat[[i]]$stratified[, c("years.num", "median", "upper", "lower", "strata")]	
	title <- NULL

	g[[i]] <- ggplot(tmp) + 
			  aes(x = years.num, y = median * 1000, ymin = lower * 1000, ymax = upper * 1000, color = strata, fill = strata, group = strata) + 
			  geom_ribbon(alpha = 0.2, color = NA) + 
			  geom_line()  + 
			  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=1) +
			  theme_bw() + 
			  theme(legend.position = c(0.2, 0.8), legend.background=element_blank()) + 
			  scale_color_brewer(NULL, palette = "Set1") +
			  scale_fill_brewer(NULL, palette = "Set1", guide = FALSE) + 
			  xlab("") + ylab("Stratum-specific NMR (per 1000 births)") + ylim(c(0.01, ymax) * 1000) + 
			  ggtitle(title) 

}
g <- g[[3]] + g[[2]] + g[[1]]
ggsave(g, file = paste0(dir, "Malawi-national-strata.pdf"), width = 9, height = 4)


g <- NULL
my_colors <- c('#e41a1c','#377eb8','#984ea3','#1b9e77','#ff7f00','#ffff33','#a65628')
for(i in 1:length(proj.bb.nat.strat)){
	tmp <- proj.bb.nat.strat[[i]]$overall[, c("years.num", "median", "upper", "lower")]
	tmp$method = paste0("Temporal Model: ", toupper(names(proj.bb.nat.strat)[i]))
	# #  Add unstratified
	# tmp1 <- proj.bb.nat[[i]]$overall[, c("years.num", "median", "upper", "lower")]
	# tmp1$method = paste0("Temporal Smoothing: ", toupper(names(proj.bb.nat.strat)[i]), " (unstratified)")
	# tmp <- rbind(tmp, tmp1)

	#  Add UN
	tmp1 <- UN
	colnames(tmp1) <- c("years.num", "median", "upper", "lower")
	tmp1$method = "UN IGME"
	tmp <- rbind(tmp, tmp1)
	# #  Add Direct
	tmp1 <- subset(direct, region == "All")
	tmp1 <- tmp1[, c("years", "mean", "lower", "upper", "surveyYears")]
	colnames(tmp1) <- c("years.num", "median",  "upper","lower", "method")
	tmp <- rbind(tmp, tmp1)
	tmp$years.num <- as.numeric(tmp$years.num)
	tmp <- subset(tmp, years.num >= 2000)
 	tmp$method <- factor(tmp$method, levels = sort(unique(tmp$method))[c(1,2,4,3)])

	# title <- paste0("Temporal smoothing: ", toupper(names(proj.bb.nat.strat)[i]))
	title <- NULL

	g[[i]] <- ggplot(tmp) + 
			  aes(x = years.num, y = median * 1000, ymin = lower * 1000, ymax = upper * 1000, color = method, fill = method, group = method) + 
			  geom_ribbon(data = subset(tmp, !method %in% c("DHS 2010", "DHS 2015")), alpha = 0.2, color = NA) + 
			  geom_line(data = subset(tmp, !method %in% c("DHS 2010", "DHS 2015"))) + 
			  geom_point(data = subset(tmp, method %in% c("DHS 2010", "DHS 2015"))) + 
			  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=1) +
			  theme_bw() + 
			  theme(legend.position = c(0.4, 0.8), legend.background=element_blank()) + 
			  scale_color_brewer(NULL, palette = "Set1") +
			  scale_fill_brewer(NULL, palette = "Set1", guide = FALSE) + 
			  xlab("")  + ylab("National NMR (per 1000 births)") + ylim(c(0.01, ymax) * 1000) + 
			  ggtitle(title) 

}
g <- g[[3]] + g[[2]] + g[[1]]
ggsave(g, file = paste0(dir, "Malawi-national-compare.pdf"), width = 9, height = 4)



# Subnational estimates
tmp <- proj.bb.strat[["ar1"]]$stratified[, c("region", "years.num", "median", "upper", "lower", "strata")]
g <- ggplot(subset(tmp, years.num >= 2000)) + 
			  aes(x = years.num, y = median * 1000, ymin = lower * 1000, ymax = upper * 1000, color = strata, fill = strata, group = strata) + 
			  geom_ribbon( alpha = 0.2, color = NA) + 
			  geom_line(size = 1) + 
			  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=1) +
			  facet_wrap(~region, ncol = 7) + 
			  theme_bw() + 
			  theme(legend.position = "bottom", legend.background=element_blank()) + 
			  scale_color_brewer(NULL, palette = "Set1") +
			  scale_fill_brewer(NULL, palette = "Set1", guide = FALSE) + 
			  scale_size_continuous(guide = FALSE) + 
			  xlab("") + ylab("Stratum-specific NMR (per 1000 births)") 
ggsave(g, file = paste0(dir, "Malawi-compare-stratified.pdf"), width = 10, height = 7) 


tmp <- proj.bb.strat[["ar1"]]$overall[, c("region", "years.num", "median", "upper", "lower")]
tmp$method <- "Beta-Binomial Model"
# add unstratified model
# tmp1 <- proj.bb[["ar1"]]$overall[, c("region", "years.num", "median", "upper", "lower")]
# tmp1$method <- "Unstratified Model"
# tmp <- rbind(tmp, tmp1)
tmp1 <- subset(direct.4yr, region != "All")
tmp1 <- tmp1[, c("region", "years", "mean", "lower", "upper", "surveyYears")]
tmp1$years <- recode(tmp1$years, "00-03" = 2001.5, "04-07" = 2005.5, "08-11" = 2009.5, "12-15" = 2013.5)
colnames(tmp1) <- c("region", "years.num", "median",  "upper","lower", "method")
# tmp1$upper <- NA
# tmp1$lower <- NA
tmp <- rbind(tmp, tmp1)
tmp$years.num <- as.numeric(tmp$years.num)
tmp <- subset(tmp, years.num >= 2000)
ylim <- range(c(tmp$median), na.rm = TRUE) * 1000
tmp$method <- factor(tmp$method, levels = c("DHS 2010", "DHS 2015", "Beta-Binomial Model"))
g <- ggplot(tmp, aes(x = years.num, y = median * 1000, ymin = lower * 1000, ymax = upper * 1000, 
					 color = method, fill = method, group = method)) + 
			  # geom_ribbon(data = subset(tmp, method %in% c("DHS 2010", "DHS 2015")), alpha = 0.1, color = NA) + 
			  geom_line(size = 0.5) + 
			  geom_line(data = subset(tmp, !method %in% c("DHS 2010", "DHS 2015")), size = 1) + 
			  geom_ribbon(data = subset(tmp, !method %in% c("DHS 2010", "DHS 2015")), alpha = 0.2, color = NA) + 
			  geom_point(data = subset(tmp, method %in% c("DHS 2010", "DHS 2015")), size = 0.8) + 
			  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=1) +
			  facet_wrap(~region, ncol = 7) + 
			  theme_bw() + 
			  theme(legend.position = "bottom", legend.background=element_blank()) + 
			  scale_color_brewer(NULL, palette = "Set1") +
			  scale_fill_brewer(NULL, palette = "Set1", guide = FALSE) + 
			  scale_size_continuous(guide = FALSE) + 
			  coord_cartesian(ylim = ylim) + xlab("")  + ylab("NMR (per 1000 births)") 
ggsave(g, file = paste0(dir, "Malawi-compare.pdf"), width = 12, height = 7)


tmp <- subset(tmp, method == "Beta-Binomial Model")
g <- ggplot(tmp, aes(x = years.num, y = median * 1000, ymin = lower * 1000, ymax = upper * 1000)) + 
			  # geom_ribbon(data = subset(tmp, method %in% c("DHS 2010", "DHS 2015")), alpha = 0.1, color = NA) + 
			  geom_line(size = 1, color = "#4daf4a") + 
			  # geom_line(data = subset(tmp, !method %in% c("DHS 2010", "DHS 2015")), size = 1) + 
			  geom_ribbon(fill =  "#4daf4a", alpha = 0.2, color = NA) + 
			  # geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=1) +
			  facet_wrap(~region, ncol = 7) + 
			  theme_bw() + 
			  theme(legend.position = "none", legend.background=element_blank()) + 
			  xlab("")  + ylab("NMR (per 1000 births)") 
ggsave(g, file = paste0(dir, "Malawi-compare-only.pdf"), width = 12, height = 7)

 

g <- mapPlot(subset(proj.bb.strat[["ar1"]]$overall, 
					years.num %in% seq(2001, 2019, by = 3)), 
			by.data = "region", geo = MalawiMap, by.geo = "ADM2_EN", 
			values = c("median"), variables = "years", is.long=TRUE,
			ncol = 7, legend.label = "NMR", direction = -1, per1000 = TRUE, 
			border = "gray70", size = 0.1)
ggsave(g, file = paste0(dir, "Malawi-NMR-map.pdf"), width = 10, height = 4)
###
### Full compare
###
tmp <- proj.bb.strat[["ar1"]]$overall[, c("region", "years.num", "median", "upper", "lower")]
tmp$method <- "Beta-Binomial Model"
tmp1 <- subset(smooth.direct, region != "All")
tmp1 <- tmp1[, c("region", "years.num", "median", "upper", "lower")]
tmp1$method <- "Smoothed Direct Model"
tmp <- rbind(tmp, tmp1)
tmp1 <- subset(direct.4yr, region != "All")
tmp1 <- tmp1[, c("region", "years", "mean", "lower", "upper", "surveyYears")]
tmp1$years <- recode(tmp1$years, "00-03" = 2001.5, "04-07" = 2005.5, "08-11" = 2009.5, "12-15" = 2013.5)
colnames(tmp1) <- c("region", "years.num", "median",  "upper","lower", "method")
tmp1$lower <- tmp1$upper <- NA
tmp <- rbind(tmp, tmp1)
tmp$years.num <- as.numeric(tmp$years.num)
tmp <- subset(tmp, years.num >= 2000)
ylim <- range(c(tmp$median), na.rm = TRUE) * 1000
tmp$method <- factor(tmp$method, levels = c("DHS 2010", "DHS 2015", "Beta-Binomial Model", "Smoothed Direct Model"))
g <- ggplot(tmp) + 
			  aes(x = years.num, y = median * 1000, ymin = lower * 1000, ymax = upper * 1000, color = method, fill = method, group = method) + 
			  # geom_ribbon(data = subset(tmp, method %in% c("DHS 2010", "DHS 2015")), alpha = 0.1, color = NA) + 
			  geom_line() + 
			  geom_line(data = subset(tmp, !method %in% c("DHS 2010", "DHS 2015")), size = 1) + 
			  geom_ribbon(data = subset(tmp, !method %in% c("DHS 2010", "DHS 2015")), alpha = 0.2, color = NA) + 
			  geom_vline(xintercept = 2015, linetype=3, color = "gray10", size=1) +
			  facet_wrap(~region, ncol = 7) + 
			  theme_bw() + 
			  theme(legend.position = "bottom", legend.background=element_blank()) + 
			  scale_color_brewer(NULL, palette = "Set1") +
			  scale_fill_brewer(NULL, palette = "Set1", guide = FALSE) + 
			  scale_size_continuous(guide = FALSE) + 
			  coord_cartesian(ylim = ylim) + xlab("")  + ylab("NMR (per 1000 births)") 
ggsave(g, file = paste0(dir, "Malawi-compare-all.pdf"), width = 12, height = 7)

###
### Odds ratio
###
T <- 20
samples <- NULL
for(i in 1:length(proj.bb.strat[["ar1"]]$draws)){
	r <- proj.bb.strat[["ar1"]]$draws[[i]]$latent[paste0("time.struct:", 1:T), ]
	u <- proj.bb.strat[["ar1"]]$draws[[i]]$latent[paste0("time.struct:", (T + 1):(2*T)), ]
	r0 <- proj.bb.strat[["ar1"]]$draws[[i]]$latent["age.intercept0:rural:1", ]
	u0 <- proj.bb.strat[["ar1"]]$draws[[i]]$latent["age.intercept0:urban:1", ]
	odds <-  exp(u + u0) / exp(r + r0)
	samples <- rbind(samples, data.frame(odds.ratio = odds, time = (1:T) + 1999, itr = i))
}
agg <-  samples %>% 
			group_by(time) %>%
			summarise(
				l = quantile(odds.ratio, .975, na.rm = T),
				u = quantile(odds.ratio, .025, na.rm = T),
				m = median(odds.ratio, na.rm = T))
g <- ggplot(agg, aes(x = time, y = m, ymin = l, ymax = u)) + 
	 geom_line(color = '#e41a1c') + 
	 geom_ribbon(color = NA, fill = '#e41a1c', alpha = 0.1) + 
	 geom_hline(yintercept = 1, linetype = 'dashed') + 
	 theme_bw() + 
	 xlab("") + ylab("Odds ratio")
ggsave(g, file = paste0(dir, "Malawi-UR-odds.pdf"), width = 5, height = 4.5)


###
### Rank plots
###
samples <- NULL
for(y in 2000:2019){
	for(re in unique(proj.bb.strat[["ar1"]]$overall$region)){
   		 cond1 <- unlist(lapply(proj.bb.strat[["ar1"]]$draws.est, function(x){x$years == y}))
   		 cond2 <- unlist(lapply(proj.bb.strat[["ar1"]]$draws.est, function(x){x$region == re}))
   		 index <- intersect(which(cond1), which(cond2))
   		 u <- ifelse(proj.bb.strat[["ar1"]]$draws.est[[index[1]]]$strata == "urban", 
   		 			index[1], index[2])
   		 r <- ifelse(proj.bb.strat[["ar1"]]$draws.est[[index[1]]]$strata == "rural", 
   		 			index[1], index[2])
   		 q <- weights[which(weights$years == y & weights$region == re), "urban"]
   		 s <- proj.bb.strat[["ar1"]]$draws.est[[u]]$draws * q + 
   		 	  proj.bb.strat[["ar1"]]$draws.est[[r]]$draws * (1 - q) 
		samples <- rbind(samples, data.frame(draws = s, years = y, region = re, itr = 1:1000))
	}
}
samples$ranking <- NA
for(i in 1:1000){
	for(y in 2000:2019){ 
		sub <- which(samples$years == y & samples$itr == i)
		ra <- order(samples[sub,]$draw, decreasing = FALSE)
		samples[sub[ra], "ranking"] <- 1:length(sub)
	}
}
pdf(paste0(dir, "Malawi-ranking.pdf"), width = 12, height = 8)
for(y in 2000:2019){
	g <- ggplot(subset(samples, years == y), aes(x = ranking)) + 
		 geom_histogram(breaks = c(1:28), fill = "#377eb8", col="white") + 
		 facet_wrap(~region, ncol = 7) + 
		 xlab("") + ylab("") + 
		 theme_bw() + scale_y_continuous(breaks=NULL) + 
		 ggtitle(paste0("Ranking of districts in ", y))
	 print(g)

}
dev.off()


###
### Save output
###
output <- proj.bb.strat[["ar1"]]$overall
output <- output[, c("region", "years", "median", "upper", "lower")]
colnames(output)[3] <- "estimate"
output$model <- "BB8"
output2 <- direct.comb
output2 <- output2[, c("region", "years", "mean", "lower", "upper")]
output2 <- subset(output2, region != "All")
colnames(output2)[3] <- "estimate"
output2$model <- "Combined Direct"
output <- rbind(output, output2)
output$estimate <- output$estimate * 1000
output$lower <- output$lower * 1000
output$upper <- output$upper * 1000
write.csv(output, file = "../Results/Malawi-NMR-admin2-2000-2019.csv", row.names = F, quote = F)