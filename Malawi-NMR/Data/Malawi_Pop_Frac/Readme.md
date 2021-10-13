This folder contains the processed urban/rural proportion for the NMR analysis. The dataset contains only number of births and neonatal deaths aggregated at the cluster-level. The RData object can be generated with the script ``../Codes/CaseStudy_Malawi_CalculateUR.R`` and ``../Codes/CaseStudy_Malawi_CalculateUR_1998.R``.

``pop_mw_under1.RData``
    + prop: Data frame of time-varying urban/rural proportion by adm2 areas.
    + census: Data frame of census reported urban/rural proportion in Table 2 of 2008 Census.

``Table2-Census2008.csv``
    + Used in the main analysis to threshold urban/rural  

``pop_mw_under1_1998.RData``
    + prop: Data frame of time-varying urban/rural proportion by adm2 areas, according to the definition of 1998 census.
    + census: Data frame of census reported urban/rural proportion in Table 1 of 2008 Census.

``Table1-Census1998.csv``
    + Used in the additional analysis in the supplement.
