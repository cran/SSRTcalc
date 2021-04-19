<!-- badges: start -->
[![R-CMD-check](https://github.com/agleontyev/SSRTcalc/workflows/R-CMD-check/badge.svg)](https://github.com/agleontyev/SSRTcalc/actions)
<!-- badges: end -->

# SSRTcalc

Tools to estimate stop-signal reaction time in R

## Introduction
This package contains functions that allow for easy estimation of stop-signal reaction time (SSRT), obtained from stop-signal task experiments.
For this package to work, the data has to be in the long format (i.e., each row represents one trial for one individual). It is possible to apply functions by an individual, if the results are in one big dataframe. For example:

```{r}
sapply(split(df, df$SubjID), integration_adaptiveSSD, stop_col = 'vol',rt_col = 'RT_exp', acc_col = 'correct', ssd_col = 'soa')
```

NB: the package is in beta version. 

## Installation
Run the following in your R/RStudio console:

```{r}
install.packages("devtools")
devtools::install_github("agleontyev/SSRTcalc")
```

## Contact me  

Any questions/concerns/suggestions are welcome. Contact me at a.g.leontiev@tamu.edu or a.g.leontiev@gmail.com
