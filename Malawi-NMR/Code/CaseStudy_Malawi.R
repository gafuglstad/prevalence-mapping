##
##  Running the analysis scripts in the following order
##


# Load input data
#   The input data can also be created by source("CaseStudy_Malawi_Data.R")
#   Registration with DHS is needed to use the script to generate the data
load("../Malawi_AGG_DHS/Malawi-Agg-Data.RData")

# Organize comparison data
igme <- read.csv("../Data/IGME/UN_IGME_Malawi.csv")
UN <- data.frame(year = as.numeric(gsub("-06", "", igme$TIME_PERIOD)), 
				 est.un = igme$OBS_VALUE / 1000, 
				 lower.un = igme$LOWER_BOUND / 1000, 
				 upper.un = igme$UPPER_BOUND / 1000 )


# Fit model
source("CaseStudy_Malawi_Analysis_BetaBinomial.R")

# Make figures
source("CaseStudy_Malawi_Figures.R")
