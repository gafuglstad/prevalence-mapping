This folder contains various R scripts to reproduce the case study of Malawi NMR. Some scripts require access to the full micro-level DHS data or the large population raster files. These scripts are included only for reproducibility purpose. The aggregated/processed results of these scripts are included in this repository and can be used to reproduce the main analysis and figures directly.

### Scripts without micro-level data
``CaseStudy_Malawi.R``
    + Main script to implement the models
    
``CaseStudy_Malawi_Analysis_BetaBinomial.R``
    + Script to fit the beta-binomial models in the main analysis

``CaseStudy_Malawi_Direct_Estimation.R``
    + Script to obtain the direct and smoothed direct estimates
    + Requires registration with DHS to access the full micro-level DHS data
    + These estimates are saved in ``../Results`` folder without the raw input data

``CaseStudy_Malawi_Figures.R``
    + Script to make the figures in the paper and supplementary materials

### Scripts that require micro-level data

``CaseStudy_Malawi_Data.R``
    + Script to download and process the micro-level DHS data
    + Processed data are saved in ``../Data/``

``CaseStudy_Malawi_CalculateUR.R``
    + Script to download and process the population raster
    + Processed data are saved in ``../Data/``

``CaseStudy_Malawi_CalculateUR_1998.R``
    + Script to download and process the population raster for the urbanicity defined in the 1998 census. 
    + Additional analysis in the supplementary materials.
    + Processed data are saved in ``../Data/``

``CaseStudy_Malawi_Analysis_3surveys.R.R``
    + Script to fit the beta-binomial models in the supplementary analysis with three surveys.




