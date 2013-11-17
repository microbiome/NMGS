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
install_github("microbiome", "NMGS")
library(NMGS)
```


### Learning the NMGS model

Carry out the NMGS simulations with the [C
scripts](https://github.com/microbiome/NMGS/tree/master/C).


### Reading the NMGS outputs

Read the output files in R for further analysis. We assume you have
the C output files in the working directory.


```r
samples <- read_nmgs("Simulation_out.csv")
```

```
## Error: could not find function "read_nmgs"
```

```r
metacommunity <- read_nmgs_metacommunity("Simulation_out_m.csv")
```

```
## Error: could not find function "read_nmgs_metacommunity"
```

```r
stats <- read_nmgs_stats("Simulation_out_s.csv")
```

```
## Error: could not find function "read_nmgs_stats"
```


### Analyze the NMGS outputs

Calculate medians for all samples:


```r
medians <- nmgs_median(samples, burnin = 25000)
```

```
## Error: could not find function "nmgs_median"
```

```r
print(medians)
```

```
## Error: object 'medians' not found
```


Test neutrality:


```r
# Compare the log-likelihood under the full neutral model with the
# observed data.
print(nmgs_neutrality(stats, "LN", "LO"))
```

```
## Error: could not find function "nmgs_neutrality"
```

```r

# The same but for local community assembly.
print(nmgs_neutrality(stats, "LL", "LO"))
```

```
## Error: could not find function "nmgs_neutrality"
```


Calculate metapopulation average (currently only for the local communities):


```r
print(nmgs_metapopulation_average(metacommunity)$local)
```

```
## Error: could not find function "nmgs_metapopulation_average"
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
## R version 3.0.1 (2013-05-16)
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
## [1] knitr_1.2   NMGS_1.0    Rcpp_0.10.3
## 
## loaded via a namespace (and not attached):
## [1] digest_0.6.3   evaluate_0.4.3 formatR_0.7    stringr_0.6.2 
## [5] tools_3.0.1
```





