---
output: github_document
---

# Purpose

This is a package for MR-rr by Yuexiang Peng (author, maintainer) and Zhilong Zhang (author).

- The URL to the GitHub (i.e., the source code) is: https://github.com/MatthewYuexiangPeng/MR.rr
- The URL to the Pkgdown webpage is: https://matthewyuexiangpeng.github.io/MR.rr/

# How to install
This package is called `MR.rr`. To install, run the following code (in R):

```R
library(devtools)
devtools::install_github("MatthewYuexiangPeng/MR.rr")
```

Upon completion, you can run the following code (in R):
```R
library(MR.rr")
```

# Dependencies

The package depends on the following packages: `MASS`, `ggplot2`, `patchwork`, and `pheatmap`.

# Session info

This package was developed in the following environment
```R
> devtools::session_info()
─ Session info ───────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.3.1 (2023-06-16 ucrt)
 os       Windows 11 x64 (build 22631)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  Chinese (Simplified)_China.utf8
 ctype    Chinese (Simplified)_China.utf8
 tz       America/Los_Angeles
 date     2024-11-18
 rstudio  2023.09.0+463 Desert Sunflower (desktop)
 pandoc   3.1.1 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────
 ! package      * version    date (UTC) lib source
   askpass        1.2.0      2023-09-03 [1] CRAN (R 4.3.2)
   brio           1.1.4      2023-12-10 [1] CRAN (R 4.3.2)
   cachem         1.0.8      2023-05-01 [1] CRAN (R 4.3.2)
   callr          3.7.5      2024-02-19 [1] CRAN (R 4.3.2)
   cli            3.6.2      2023-12-11 [1] CRAN (R 4.3.2)
   colorspace     2.1-0      2023-01-23 [1] CRAN (R 4.3.2)
   commonmark     1.9.0      2023-03-17 [1] CRAN (R 4.3.2)
   credentials    2.0.1      2023-09-06 [1] CRAN (R 4.3.2)
   curl           6.0.0      2024-11-05 [1] CRAN (R 4.3.3)
   desc           1.4.3      2023-12-10 [1] CRAN (R 4.3.2)
   devtools       2.4.5      2022-10-11 [1] CRAN (R 4.3.3)
   digest         0.6.34     2024-01-11 [1] CRAN (R 4.3.2)
   dplyr          1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
   ellipsis       0.3.2      2021-04-29 [1] CRAN (R 4.3.2)
   evaluate       1.0.1      2024-10-10 [1] CRAN (R 4.3.3)
   fansi          1.0.6      2023-12-08 [1] CRAN (R 4.3.2)
   fastmap        1.1.1      2023-02-24 [1] CRAN (R 4.3.2)
   fs             1.6.3      2023-07-20 [1] CRAN (R 4.3.2)
   generics       0.1.3      2022-07-05 [1] CRAN (R 4.3.2)
   gert           2.0.1      2023-12-04 [1] CRAN (R 4.3.2)
 V ggplot2        3.4.4      2024-04-23 [1] CRAN (R 4.3.3) (on disk 3.5.1)
   gh             1.4.0      2023-02-22 [1] CRAN (R 4.3.2)
   gitcreds       0.1.2      2022-09-08 [1] CRAN (R 4.3.2)
   glue           1.7.0      2024-01-09 [1] CRAN (R 4.3.2)
 V gtable         0.3.4      2024-10-25 [1] CRAN (R 4.3.3) (on disk 0.3.6)
   htmltools      0.5.8.1    2024-04-04 [1] CRAN (R 4.3.3)
   htmlwidgets    1.6.4      2023-12-06 [1] CRAN (R 4.3.2)
   httpuv         1.6.14     2024-01-26 [1] CRAN (R 4.3.2)
   httr2          1.0.6      2024-11-04 [1] CRAN (R 4.3.3)
   jsonlite       1.8.9      2024-09-20 [1] CRAN (R 4.3.3)
   knitr          1.49       2024-11-08 [1] CRAN (R 4.3.3)
   later          1.3.2      2023-12-06 [1] CRAN (R 4.3.2)
   lifecycle      1.0.4      2023-11-07 [1] CRAN (R 4.3.2)
   magrittr       2.0.3      2022-03-30 [1] CRAN (R 4.3.2)
   MASS           7.3-60.0.1 2024-01-13 [1] CRAN (R 4.3.3)
   memoise        2.0.1      2021-11-26 [1] CRAN (R 4.3.2)
   mime           0.12       2021-09-28 [1] CRAN (R 4.3.1)
   miniUI         0.1.1.1    2018-05-18 [1] CRAN (R 4.3.2)
   MR.rr        * 0.0.0.9000 2024-11-18 [1] local
 V munsell        0.5.0      2024-04-01 [1] CRAN (R 4.3.3) (on disk 0.5.1)
   openssl        2.1.1      2023-09-25 [1] CRAN (R 4.3.2)
 V patchwork      1.2.0      2024-09-16 [1] CRAN (R 4.3.3) (on disk 1.3.0)
   pheatmap       1.0.12     2019-01-04 [1] CRAN (R 4.3.3)
   pillar         1.9.0      2023-03-22 [1] CRAN (R 4.3.2)
   pkgbuild       1.4.3      2023-12-10 [1] CRAN (R 4.3.2)
   pkgconfig      2.0.3      2019-09-22 [1] CRAN (R 4.3.2)
   pkgdown        2.1.1      2024-09-17 [1] CRAN (R 4.3.3)
   pkgload        1.3.4      2024-01-16 [1] CRAN (R 4.3.2)
   prettyunits    1.2.0      2023-09-24 [1] CRAN (R 4.3.2)
   processx       3.8.3      2023-12-10 [1] CRAN (R 4.3.2)
   profvis        0.3.8      2023-05-02 [1] CRAN (R 4.3.2)
   promises       1.2.1      2023-08-10 [1] CRAN (R 4.3.2)
   ps             1.7.6      2024-01-18 [1] CRAN (R 4.3.2)
   purrr          1.0.2      2023-08-10 [1] CRAN (R 4.3.3)
   R6             2.5.1      2021-08-19 [1] CRAN (R 4.3.2)
   rappdirs       0.3.3      2021-01-31 [1] CRAN (R 4.3.2)
   rcmdcheck      1.4.0      2021-09-27 [1] CRAN (R 4.3.2)
   RColorBrewer   1.1-3      2022-04-03 [1] CRAN (R 4.3.1)
   Rcpp           1.0.12     2024-01-09 [1] CRAN (R 4.3.2)
   remotes        2.4.2.1    2023-07-18 [1] CRAN (R 4.3.2)
   rlang          1.1.3      2024-01-10 [1] CRAN (R 4.3.2)
   rmarkdown      2.29       2024-11-04 [1] CRAN (R 4.3.3)
   roxygen2       7.3.1      2024-01-22 [1] CRAN (R 4.3.2)
   rprojroot      2.0.4      2023-11-05 [1] CRAN (R 4.3.2)
   rstudioapi     0.15.0     2023-07-07 [1] CRAN (R 4.3.2)
   scales         1.3.0      2023-11-28 [1] CRAN (R 4.3.2)
   sessioninfo    1.2.2      2021-12-06 [1] CRAN (R 4.3.2)
   shiny          1.8.0      2023-11-17 [1] CRAN (R 4.3.2)
   stringi        1.8.3      2023-12-11 [1] CRAN (R 4.3.2)
   stringr        1.5.1      2023-11-14 [1] CRAN (R 4.3.2)
   sys            3.4.2      2023-05-23 [1] CRAN (R 4.3.2)
   testthat       3.2.1      2023-12-02 [1] CRAN (R 4.3.2)
   tibble         3.2.1      2023-03-20 [1] CRAN (R 4.3.2)
   tidyselect     1.2.0      2022-10-10 [1] CRAN (R 4.3.2)
   urlchecker     1.0.1      2021-11-30 [1] CRAN (R 4.3.2)
   usethis        3.0.0      2024-07-29 [1] CRAN (R 4.3.3)
   utf8           1.2.4      2023-10-22 [1] CRAN (R 4.3.2)
   vctrs          0.6.5      2023-12-01 [1] CRAN (R 4.3.2)
 V withr          3.0.0      2024-10-28 [1] CRAN (R 4.3.3) (on disk 3.0.2)
   xfun           0.49       2024-10-31 [1] CRAN (R 4.3.3)
   xml2           1.3.6      2023-12-04 [1] CRAN (R 4.3.2)
   xopen          1.0.0      2018-09-17 [1] CRAN (R 4.3.2)
   xtable         1.8-4      2019-04-21 [1] CRAN (R 4.3.2)
   yaml           2.3.8      2023-12-11 [1] CRAN (R 4.3.2)

 [1] C:/Users/Yuexiang Peng/AppData/Local/R/win-library/4.3
 [2] C:/Program Files/R/R-4.3.1/library

 V ── Loaded and on-disk version mismatch.

──────────────────────────────────────────────────────────────────────────────────────────────────────────
```
