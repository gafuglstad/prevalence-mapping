<h1>Code for "The Two Cultures for Prevalence Mapping: Small Area Estimation and Spatial Statistics"</h1>
<a href="https://arxiv.org/abs/2110.09576">arXiv preprint</a>

<h2> Description </h2>
The paper concerns spatial and spatio-temporal estimation of prevalences in administrative areas. This repository contains
the code for two case studies: spatial estimation of vaccination coverage in Nigeria in 2018, and spatio-temporal estimation of neonatal mortality rate (NMR) in Malawi in 2000&ndash;2019.

<h2> Case studies </h2>
<h3> Spatial estimation of vaccination coverage </h3>
Start by running "Nigeria-Vac/Code/CaseStudy_Nigeria.R".

<h4> Goal </h4>
We consider vaccination coverage for the first dose of measles-containing-vaccine (MCV1)
among children aged 12&ndash;23 months in Nigeria. The goal of the analysis is estimation of MCV1 coverage for
admin1, which consists of 36 states and the federal capital area, 
and admin2, which consists of 774 local government areas (LGAs), and to this end we
analyze the 2018 Nigeria Demographic and Health Survey.

<h4> Code files </h4>
<ul>
    <li> <b>Nigeria-Vac/Code/CaseStudy_Nigeria.R:</b> Main script. Runs all models and validation and creates figures and data for tables. </li>
    <li> <b>Nigeria-Vac/Code/functions.R:</b> Support code. Provides functions for fitting models and for aggregation. </li>
</ul>

<h4> Data files </h4>
<ul>
    <li> <b>Nigeria-Vac/Data/gadm36_NGA_shp:</b> Shape file for admin0, admin1 and admin2. </li>
    <li> <b>Nigeria-Vac/Data/Nigeria_2017_proj:</b> Urban proportions in updated sampling frame from 2017. </li>
    <li> <b>Nigeria-Vac/Data/Nigeria_AGG_DHS:</b> Data extracted from Nigeria 2018 DHS survey.</li>
    <li> <b>Nigeria-Vac/Data/Nigeria_pop:</b> Population rasters.</li>
    <li> <b>Nigeria-Vac/Data/preparedCovRasters:</b> Prepared covariate rasters.</li>
</ul>

<h4> Data sources </h4>
<ul>
    <li>WorldPop</li>
    <ul>
        <li> Poverty map for 2010 ($2 a day limit): https://www.worldpop.org/geodata/summary?id=1267</li>
        <li> Population (all age, all sex) 2017: https://data.worldpop.org/GIS/Population/Global_2000_2020/2017/NGA/nga_ppp_2017_UNadj.tif</li>
        <li> Population (all age, all sex) 2018: https://data.worldpop.org/GIS/Population/Global_2000_2020/2018/NGA/nga_ppp_2018_UNadj.tif</li>
        <li> Population (1&ndash;4 years, female) 2018: https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2018/NGA/nga_f_1_2018.tif</li>
        <li> Population (1&ndash;4 years, male) 2018: https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2018/NGA/nga_m_1_2018.tif</li>
        <li> Population (1&ndash;4 years, female) 2006: https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2006/NGA/nga_f_1_2006.tif</li>
        <li> Population (1&ndash;4 years, male) 2006: https://data.worldpop.org/GIS/AgeSex_structures/Global_2000_2020/2006/NGA/nga_m_1_2006.tif</li>
        <li> License: https://creativecommons.org/licenses/by/4.0/</li>
        <li> Population rasters are automatically downloaded when running "Nigeria-Vac/Code/CaseStudy_Nigeria.R"</li>
    </ul>
    <li>The Malaria Atlas Project</li>
    <ul>
        <li> Accessibility to cities 2015: https://malariaatlas.org/research-project/accessibility-to-cities/</li>
        <li> License: http://creativecommons.org/licenses/by/4.0/ </li>
    </ul>
    <li>GADM</li>
    <ul>
        <li> Shape files of Nigeria (version 3.6; admin0, admin1 and admin2): https://biogeo.ucdavis.edu/data/gadm3.6/shp/gadm36_NGA_shp.zip</li>
        <li> License: Academic use only. No commercial use or redistribution. See "Nigeria-Vac/Data/gadm36_NGA_shp/license.txt" for details.</li>
	<li> Files for Nigeria included with permission from GADM</li>
    </ul>
    <li>DHS</li>
    <ul>
        <li> Urban proportions in updated sampling frame for 2018 DHS survey:  Personal correspondence with DHS</li>
        <li> Nigeria 2018 DHS survey (NDHS2018) + displaced GPS coordinates: https://dhsprogram.com/ </li>
        <li> Acquired permission to share aggregated cluster data and displaced GPS coordinates. </li>
    </ul>
</ul>

<h3> Spatiotemporal estimation of Neonatal Mortality Rates (NMR) </h3>
Start by running "Nigeria-Vac/Code/CaseStudy_Malawi.R".

<h4> Goal </h4>
The goal of the analysis is estimation of NMR for admin2 areas, which consists of 28 regions from  2000 to 2015, and short term projection to 2019. To this end we analyze the 2010 and 2015-16 Malawi Demographic and Health Survey.

<h4> Code files </h4>
<ul>
    <li> <b>Malawi-NMR/Code/CaseStudy_Malawi.R:</b> Main script. Runs all models and creates figures and data for tables. </li>
    <li> Additional files for processing various raw input data and additional analysis, explained in <b>Malawi-NMR/Code/Readme.md</b> </li>
</ul>

<h4> Data files </h4>
<ul>
    <li> <b>Malawi-NMR/Data/shapefiles:</b> Shape file for admin2 regions.</li>
    <li> <b>Malawi-NMR/Data/IGME:</b> UN-IGME estimates of national NMR.</li>
    <li> <b>Malawi-NMR/Data/Malawi_AGG_DHS:</b> Data extracted from Malawi DHS surveys.</li>
    <li> <b>Malawi-NMR/Data/Malawi_Pop_Frac:</b> Urban/rural proportions computed from the WorldPop raster files. </li>
    <li> <b>Malawi-NMR/Data/Malawi_FULL_DHS:</b> Micro-level data extracted from Malawi DHS surveys. Not included in the repository due to data privacy. See the folder for details on how to download and process the data. </li>
    <li> <b>Malawi-NMR/Data/WorldPop-Population:</b> Population raster files from WorldPop. Not included in the repository due to size of the files. See the folder for details on how to download and process the data. </li>
</ul>

<h4> Data sources </h4>

See individual folder for details.
        

