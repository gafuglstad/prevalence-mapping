# prevalence-mapping
Case study Nigeria

Code files:
    Nigeria-Vac/Code/CaseStudy_Nigeria.R:
        Main script. Runs all models and validation and creates figures and data for tables.
    Nigeria-Vac/Code/functions.R:
        Support code. Provides functions for fitting models and for aggregation.

Data sources:
    WorldPop:
        Licence: https://creativecommons.org/licenses/by/4.0/
        Poverty map for 2010 ($2 a day limit)
            Link: https://www.worldpop.org/geodata/summary?id=1267
        Population 2017:
            Link: https://data.worldpop.org/GIS/Population/Global_2000_2020/2017/NGA/nga_ppp_2017_UNadj.tif
        Children (male/female) aged 1--4 years.
            Link (female): https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2018/NGA/nga_f_1_2018.tif
            Link (male): https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2018/NGA/nga_m_1_2018.tif
    The Malaria Atlas Project
        Licence: http://creativecommons.org/licenses/by/4.0/
        Accessibility to cities 2015
        Link: https://malariaatlas.org/research-project/accessibility-to-cities/
    GADM:
        Licence: Academic use only. No commercial use.
        Shape files of Nigeria (version 3.6): admin0, admin1 and admin2
        Link: https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_NGA_shp.zip
    DHS:
        Projected urban proportions at admin2 -- Personal correspondence
        NDHS2018 survey + displaced GPS coordinates
        Licence: Can share proportions and aggregated data from survey
        

