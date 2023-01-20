# Multiscale spatially varying coeﬀicient modelling using a Geographical Gaussian Process GAM

Alexis Comber<sup>1*</sup> Paul Harris<sup>2</sup> and Chris Brunsdon<sup>3</sup> 

<sup>1</sup> School of Geography, University of Leeds, Leeds, UK.\
<sup>2</sup> Sustainable Agriculture Sciences North Wyke, Rothamsted Research, Okehampton, UK.
<sup>2</sup> National Centre for Geocomputation, Maynooth University, Ireland
<sup>*</sup> contact author: a.comber@leeds.ac.uk

## Abstract
This paper proposes a novel spatially varying coefficient (SVC) regression through a Geographical Gaussian Process GAM (GGP-GAM): a Generalized Additive Model (GAM) with Gaussian Process (GP) splines parameterised at observation locations. As with Multiscale Geographically Weighted Regression (MGWR), the proposed GGP-GAM specification is multiscale but has fewer theoretical and technical limitations. Currently, MGWR is the leading SVC regression and can be similarly viewed as GAM-based but with a moving window kernel parameterisation. A GGP-GAM was applied to simulated coefficient datasets exhibiting varying degrees of spatial heterogeneity and was shown to out-perform MGWR under a range of fit metrics for both coefficient and prediction accuracy. The proposed SVC model was then applied to a case study of the Brexit vote over UK parliamentary constituencies with covariates describing a range of socio-economic, life-stage and voting factors. Resultant coefficient estimates were mapped, showing regional scales of relationship non-stationarity in the UK’s Brexit voting process. For the proposed GGP-GAM, a number of areas of further work are identified. This includes the creation of a user-friendly wrapper R package to support model creation and coefficient mapping, and to facilitate ease of comparison with alternate SVC models. Calibration issues are discussed such as GAM tuning (particularly of knots and spline smoothing parameters) and the difficulty in linking the GGP-GAM spline smoothing parameters to intuitive user understandings of process spatial heterogeneity, where corresponding parameters for MGWR (i.e., its kernel bandwidths) hold an advantage.

This paper has been submitted to IJGIS.

## Code
To run the analysis in this paper you should download the the R script `gam_gpp.R`the 2 data files (`sim1_fullpaper.RData` and `simulation_result_fullpaper.RData`) and install the packages. Package and other info is below. The data files and supporting scripts will need will need to be locally available . The code recreates the results as the same sequence in the paper. 

If you have any problems with data / code / versions etc please contact Lex Comber at the email above.

```{r}
sessionInfo()
R version 4.2.0 (2022-04-22)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Monterey 12.6.1

Matrix products: default
BLAS:   /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRblas.0.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib

locale:
[1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] spdep_1.2-5       broom_1.0.1       parlitools_0.4.1  GWmodel_2.2-9    
 [5] spatialreg_1.2-5  sf_1.0-8          Matrix_1.5-1      spData_2.2.0     
 [9] Rcpp_1.0.9        robustbase_0.95-0 maptools_1.1-4    sp_1.5-0         
[13] mgcv_1.8-40       nlme_3.1-159      cols4all_0.3      cowplot_1.1.1    
[17] forcats_0.5.2     stringr_1.4.1     dplyr_1.0.10      purrr_0.3.5      
[21] readr_2.1.2       tidyr_1.2.1       tibble_3.1.8      ggplot2_3.3.6    
[25] tidyverse_1.3.2  

loaded via a namespace (and not attached):
 [1] googledrive_2.0.0   colorspace_2.0-3    deldir_1.0-6        ellipsis_0.3.2     
 [5] class_7.3-20        snakecase_0.11.0    mnis_0.3.1          fs_1.5.2           
 [9] proxy_0.4-27        farver_2.1.1        fansi_1.0.3         lubridate_1.8.0    
[13] xml2_1.3.3          codetools_0.2-18    splines_4.2.0       knitr_1.40         
[17] jsonlite_1.8.0      dbplyr_2.2.1        hansard_0.8.0       compiler_4.2.0     
[21] httr_1.4.4          backports_1.4.1     assertthat_0.2.1    gargle_1.2.1       
[25] cli_3.4.1           s2_1.1.0            tools_4.2.0         coda_0.19-4        
[29] gtable_0.3.1        glue_1.6.2          wk_0.7.0            gmodels_2.18.1.1   
[33] cellranger_1.1.0    raster_3.6-3        vctrs_0.4.2         gdata_2.18.0.1     
[37] xfun_0.33           rvest_1.0.3         lifecycle_1.0.3     gtools_3.9.3       
[41] googlesheets4_1.0.1 terra_1.6-17        DEoptimR_1.0-11     LearnBayes_2.15.1  
[45] MASS_7.3-58.1       zoo_1.8-11          scales_1.2.1        hms_1.1.2          
[49] parallel_4.2.0      expm_0.999-6        stringi_1.7.8       highr_0.9          
[53] e1071_1.7-12        boot_1.3-28         intervals_0.15.2    rlang_1.0.6        
[57] pkgconfig_2.0.3     lattice_0.20-45     labeling_0.4.2      tidyselect_1.2.0   
[61] magrittr_2.0.3      R6_2.5.1            generics_0.1.3      DBI_1.1.3          
[65] pillar_1.8.1        haven_2.5.1         foreign_0.8-82      withr_2.5.0        
[69] units_0.8-0         xts_0.12.1          abind_1.4-5         spacetime_1.2-8    
[73] janitor_2.1.0       modelr_0.1.9        crayon_1.5.1        KernSmooth_2.23-20 
[77] utf8_1.2.2          tzdb_0.3.0          grid_4.2.0          readxl_1.4.1       
[81] FNN_1.1.3.1         reprex_2.0.2        digest_0.6.29       classInt_0.4-8     
[85] munsell_0.5.0    
```
