<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{An R Markdown Vignette made with knitr}
-->

NMGS - Neutral models for community assembly
===========

### Installing and loading the development version


```r
install.packages("devtools")
library(devtools)
install_github("NMGS", "microbiome", subdir = "R/pkg/")
library(NMGS)
```


### Learning the NMGS model

Carry out the NMGS simulations with the [C
scripts](https://github.com/microbiome/NMGS/tree/master/C). These will
produce the NMGS [output
files](https://github.com/microbiome/NMGS/tree/master/C/output).


### Reading the NMGS outputs

Read the output files in R for further analysis. We assume you have
the C output files in the working directory.


```r
samples <- read_nmgs("Simulation_out.csv")
metacommunity <- read_nmgs_metacommunity("Simulation_out_m.csv")
stats <- read_nmgs_stats("Simulation_out_s.csv")
```


### Analyze the NMGS outputs

Calculate medians for all samples:


```r
medians <- nmgs_median(samples, burnin = 25000)
print(medians)
```

```
##        theta     i1    i2     i3     i4     i5    i6    i7    i8    i9
## lower  17.58 17.581 17.58 17.581 17.581 17.581 17.58 17.58 17.58 17.58
## median 21.40  2.199 13.27  8.022  8.836  8.414 12.13 23.05 13.47 19.09
## upper  25.83 25.832 25.83 25.832 25.832 25.832 25.83 25.83 25.83 25.83
##          i10   i11   i12   i13   i14   i15   i16   i17   i18   i19   i20
## lower  17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58
## median 34.32 20.49 20.23 32.91 21.81 33.42 32.47 39.18 31.13 47.72 39.79
## upper  25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83
##          i21   i22   i23   i24   i25   i26   i27   i28   i29   i30   i31
## lower  17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58
## median 34.34 39.94 47.84 43.70 44.45 40.50 52.13 65.44 37.65 60.79 49.44
## upper  25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83
##          i32   i33   i34   i35   i36   i37   i38   i39   i40   i41   i42
## lower  17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58 17.58
## median 74.26 47.84 62.69 65.06 96.50 81.19 87.16 68.11 85.92 73.24 74.72
## upper  25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83 25.83
##          i43    i44   i45    i46   i47    i48   i49   i50
## lower  17.58  17.58 17.58  17.58 17.58  17.58 17.58 17.58
## median 74.99 113.32 67.16 106.30 86.27 116.84 85.67 93.84
## upper  25.83  25.83 25.83  25.83 25.83  25.83 25.83 25.83
```


Test neutrality:


```r
# Compare the log-likelihood under the full neutral model with the
# observed data.
print(nmgs_neutrality(stats, "LN", "LO"))
```

```
##       median1       median2            n1             n pseudo.pvalue 
##    -9473.7323    -7367.6373       24.0000     2500.0000        0.0096
```

```r

# The same but for local community assembly.
print(nmgs_neutrality(stats, "LL", "LO"))
```

```
##       median1       median2            n1             n pseudo.pvalue 
##    -7343.4197    -7367.6373     1446.0000     2500.0000        0.5784
```


Calculate metapopulation average (currently only for the local communities):


```r
print(nmgs_metapopulation_average(metacommunity)$local)
```

```
##   [1] 0.0347458 0.0239980 0.0113415 0.0558154 0.0330752 0.0217902 0.0487197
##   [8] 0.0272092 0.0533554 0.0293637 0.0371854 0.0218708 0.0342013 0.0433238
##  [15] 0.0421186 0.0272166 0.0441229 0.0415109 0.0259772 0.0414386 0.0520751
##  [22] 0.0505479 0.0335348 0.0510489 0.0466658 0.0010942 0.0022312 0.0031043
##  [29] 0.0044216 0.0008509 0.0043333 0.0024925 0.0096614 0.0039874 0.0007124
##  [36] 0.0004227 0.0010531 0.0008568 0.0011608 0.0005665 0.0002745 0.0010882
##  [43] 0.0011955 0.0004347 0.0002802 0.0002851 0.0015172 0.0011871 0.0007507
##  [50] 0.0008831 0.0002929 0.0007395 0.0004119 0.0009948 0.0018817 0.0002762
##  [57] 0.0007216 0.0008666 0.0005673 0.0001367 0.0010618 0.0002776 0.0002745
##  [64] 0.0002735 0.0001350 0.0014793 0.0001346 0.0001394 0.0002763 0.0001349
##  [71] 0.0007240 0.0001352 0.0001378 0.0001405 0.0001378 0.0001371 0.0001355
##  [78] 0.0001340 0.0001444 0.0001342 0.0002860 0.0001367 0.0001403 0.0002826
##  [85] 0.0002737 0.0001335 0.0001374 0.0002759 0.0002730 0.0001367 0.0001419
##  [92] 0.0001334 0.0001330 0.0001414 0.0001292 0.0001386 0.0001384 0.0001342
##  [99] 0.0001366 0.0001365 0.0001313 0.0001351 0.0001364 0.0001383 0.0001343
## [106] 0.0001365 0.0001332 0.0001347 0.0001350 0.0001348 0.0001340 0.0001323
## [113] 0.0001320 0.0001345 0.0001342 0.0001370 0.0001318 0.0001334 0.0001363
## [120] 0.0001349 0.0001337 0.0001344 0.0001364 0.0001348 0.0001314 0.0001383
## [127] 0.0029269
```



### Licensing and Citations

This work can be freely used, modified and distributed under the 
[GNU General Public GPL>=2 license](https://en.wikipedia.org/wiki/GNU_General_Public_License).

Kindly cite the work, if appropriate, as 'Leo Lahti, Christoper Quince et al. (2013). NMGS R package. URL: https://github.com/microbiome/NMGS)'. 


### Session info

This vignette was created with


```r
sessionInfo()
```

```
## R version 2.15.2 (2012-10-26)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=C                 LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] knitr_1.2   NMGS_0.1.02
## 
## loaded via a namespace (and not attached):
## [1] digest_0.6.3   evaluate_0.4.4 formatR_0.9    stringr_0.6.2 
## [5] tools_2.15.2
```





