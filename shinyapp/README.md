# Shiny App

[![Build Status](https://travis-ci.com/saezlab/HF_meta-analysis.svg?token=PagY1pyvMyyL3AJHRy5V&branch=master)](https://travis-ci.com/saezlab/HF_meta-analysis)

```r
sessioninfo::session_info()
```

```r
─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value                       
 version  R version 4.0.0 (2020-04-24)
 os       macOS Mojave 10.14.5        
 system   x86_64, darwin17.0          
 ui       RStudio                     
 language (EN)                        
 collate  en_US.UTF-8                 
 ctype    en_US.UTF-8                 
 tz       Europe/Berlin               
 date     2020-05-01                  

─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package            * version date       lib source                                              
 AachenColorPalette * 1.1.1   2020-05-01 [1] Github (christianholland/AachenColorPalette@d4c547b)
 assertthat           0.2.1   2019-03-21 [1] CRAN (R 4.0.0)                                      
 BiocParallel         1.21.2  2019-12-21 [1] Bioconductor                                        
 cli                  2.0.2   2020-02-28 [1] CRAN (R 4.0.0)                                      
 colorspace           1.4-1   2019-03-18 [1] CRAN (R 4.0.0)                                      
 cowplot            * 1.0.0   2019-07-11 [1] CRAN (R 4.0.0)                                      
 crayon               1.3.4   2017-09-16 [1] CRAN (R 4.0.0)                                      
 crosstalk            1.1.0.1 2020-03-13 [1] CRAN (R 4.0.0)                                      
 data.table           1.12.8  2019-12-09 [1] CRAN (R 4.0.0)                                      
 digest               0.6.25  2020-02-23 [1] CRAN (R 4.0.0)                                      
 dplyr              * 0.8.5   2020-03-07 [1] CRAN (R 4.0.0)                                      
 DT                 * 0.13    2020-03-23 [1] CRAN (R 4.0.0)                                      
 ellipsis             0.3.0   2019-09-20 [1] CRAN (R 4.0.0)                                      
 fansi                0.4.1   2020-01-08 [1] CRAN (R 4.0.0)                                      
 fastmap              1.0.1   2019-10-08 [1] CRAN (R 4.0.0)                                      
 fastmatch            1.1-0   2017-01-28 [1] CRAN (R 4.0.0)                                      
 fgsea              * 1.13.4  2020-02-20 [1] Bioconductor                                        
 forcats            * 0.5.0   2020-03-01 [1] CRAN (R 4.0.0)                                      
 ggplot2            * 3.3.0   2020-03-05 [1] CRAN (R 4.0.0)                                      
 glue                 1.4.0   2020-04-03 [1] CRAN (R 4.0.0)                                      
 gridExtra            2.3     2017-09-09 [1] CRAN (R 4.0.0)                                      
 gtable               0.3.0   2019-03-25 [1] CRAN (R 4.0.0)                                      
 hms                  0.5.3   2020-01-08 [1] CRAN (R 4.0.0)                                      
 htmltools            0.4.0   2019-10-04 [1] CRAN (R 4.0.0)                                      
 htmlwidgets          1.5.1   2019-10-08 [1] CRAN (R 4.0.0)                                      
 httpuv               1.5.2   2019-09-11 [1] CRAN (R 4.0.0)                                      
 httr                 1.4.1   2019-08-05 [1] CRAN (R 4.0.0)                                      
 jsonlite             1.6.1   2020-02-02 [1] CRAN (R 4.0.0)                                      
 later                1.0.0   2019-10-04 [1] CRAN (R 4.0.0)                                      
 lattice              0.20-41 2020-04-02 [1] CRAN (R 4.0.0)                                      
 lazyeval             0.2.2   2019-03-15 [1] CRAN (R 4.0.0)                                      
 lifecycle            0.2.0   2020-03-06 [1] CRAN (R 4.0.0)                                      
 magrittr             1.5     2014-11-22 [1] CRAN (R 4.0.0)                                      
 markdown             1.1     2019-08-07 [1] CRAN (R 4.0.0)                                      
 Matrix               1.2-18  2019-11-27 [1] CRAN (R 4.0.0)                                      
 mime                 0.9     2020-02-04 [1] CRAN (R 4.0.0)                                      
 munsell              0.5.0   2018-06-12 [1] CRAN (R 4.0.0)                                      
 packrat              0.5.0   2018-11-14 [1] CRAN (R 4.0.0)                                      
 pillar               1.4.3   2019-12-20 [1] CRAN (R 4.0.0)                                      
 pkgconfig            2.0.3   2019-09-22 [1] CRAN (R 4.0.0)                                      
 plotly             * 4.9.2.1 2020-04-04 [1] CRAN (R 4.0.0)                                      
 promises             1.1.0   2019-10-04 [1] CRAN (R 4.0.0)                                      
 purrr              * 0.3.4   2020-04-17 [1] CRAN (R 4.0.0)                                      
 R6                   2.4.1   2019-11-12 [1] CRAN (R 4.0.0)                                      
 Rcpp               * 1.0.4.6 2020-04-09 [1] CRAN (R 4.0.0)                                      
 readr              * 1.3.1   2018-12-21 [1] CRAN (R 4.0.0)                                      
 rlang                0.4.5   2020-03-01 [1] CRAN (R 4.0.0)                                      
 rsconnect            0.8.16  2019-12-13 [1] CRAN (R 4.0.0)                                      
 rstudioapi           0.11    2020-02-07 [1] CRAN (R 4.0.0)                                      
 scales             * 1.1.0   2019-11-18 [1] CRAN (R 4.0.0)                                      
 sessioninfo          1.1.1   2018-11-05 [1] CRAN (R 4.0.0)                                      
 shiny              * 1.4.0.2 2020-03-13 [1] CRAN (R 4.0.0)                                      
 shinycssloaders    * 0.3     2020-01-16 [1] CRAN (R 4.0.0)                                      
 shinyhelper        * 0.3.2   2019-11-09 [1] CRAN (R 4.0.0)                                      
 shinyjs            * 1.1     2020-01-13 [1] CRAN (R 4.0.0)                                      
 shinyWidgets       * 0.5.1   2020-03-04 [1] CRAN (R 4.0.0)                                      
 stringi              1.4.6   2020-02-17 [1] CRAN (R 4.0.0)                                      
 stringr            * 1.4.0   2019-02-10 [1] CRAN (R 4.0.0)                                      
 tibble             * 3.0.1   2020-04-20 [1] CRAN (R 4.0.0)                                      
 tidyr              * 1.0.2   2020-01-24 [1] CRAN (R 4.0.0)                                      
 tidyselect           1.0.0   2020-01-27 [1] CRAN (R 4.0.0)                                      
 vctrs                0.2.4   2020-03-10 [1] CRAN (R 4.0.0)                                      
 viridisLite          0.3.0   2018-02-01 [1] CRAN (R 4.0.0)                                      
 withr                2.2.0   2020-04-20 [1] CRAN (R 4.0.0)                                      
 xfun                 0.13    2020-04-13 [1] CRAN (R 4.0.0)                                      
 xtable               1.8-4   2019-04-21 [1] CRAN (R 4.0.0)                                      
 yaml                 2.2.1   2020-02-01 [1] CRAN (R 4.0.0)   
 ```
 