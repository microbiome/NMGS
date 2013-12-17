<!--
%\VignetteEngine{knitr}
%\VignetteIndexEntry{An R Markdown Vignette made with knitr}
-->

NMGS - Neutral models for community assembly
===========

## Installing and loading the development version


```r
install.packages("devtools")
library(devtools)
install_github("NMGS", "microbiome", subdir = "R/pkg/")
library(NMGS)
```


## Learning the NMGS model

Carry out the NMGS simulations with the [C
scripts](https://github.com/microbiome/NMGS/tree/master/C). These will
produce the NMGS [output
files](https://github.com/microbiome/NMGS/tree/master/C/output).


### Reading the NMGS outputs

Read the output files in R for further analysis. We assume you have
the C output files in the working directory.


```r
library(NMGS)
samples <- read_nmgs("Simulation_out.csv")
```

```
## Error: cannot open the connection
```

```r
metacommunity <- read_nmgs_metacommunity("Simulation_out_m.csv")  # burnin and thinning already included
```

```
## Error: cannot open the connection
```

```r
stats <- read_nmgs_stats("Simulation_out_s.csv")
```

```
## Error: cannot open the connection
```


### Description of the model outputs

_Samples:_ Posterior samples of the model parameters: MCMC samples x
 parameters matrix with the following columns: MCMC sample ID, theta
 (biodiversity parameter), i1...iN immigration rates

_stats:_ Statistics on the sampled communities (after removing burn-in
and parsing the remaining samples; by default, every 10th sample
used); Matrix of selected samples x statistics: 
  * MCMC sampleID: nth sample generated under the neutral model with fitted parameters taken from the corresponding MCMC sample
  * LN, LL, LO: log-likelihoods of the full neutral sample, the local
community sample and the observed sample, respectively;
  * HN, HL, HO: species entropies (ie. Shannon diversity; neutral, local, observed)
  * SN, SL, SO: species richnesses (neutral, local, observed)

_metacommunity:_ A list with N, p, q. Gives the generated
 metacommunities: SampleN, SL, SN; p1,...,pSL; q1,...,qSN. The N is a
 matrix with MCMC sample ID, SL (species richness in the local
 assembly sample; always constant), SN (species richness in the sample
 generated under the full neutral model); p is the metacommunity
 distribution for the local assembly sample (matrix of MCMC samples x
 SL); q is the metacommunity distribution under the full neutral model
 (a list with elements with varying length SN, each corresponding to
 one MCMC sample).


## Analyzing the NMGS outputs

### Summarizing posterior samples

Posterior summaries:


```r
summaries <- nmgs_posterior_summaries(samples, burnin = 25000, thinning = 10)
```

```
## Error: could not find function "nmgs_posterior_summaries"
```

```r
head(summaries)
```

```
## Error: object 'summaries' not found
```


Visualize posterior of a given parameter (excluding burn-in and
thinning the samples):


```r
# Define burn-in, thinning and parameter to follow
burnin <- 25000  # Leave out this many samples from the start
thinning <- 10  # Take every kth samples
param <- "theta"  # or one of the migration rates i1...iN

p <- ggplot(data.frame(x = samples[seq(burnin + 1, nrow(samples), thinning), 
    param]), aes(x = x))
```

```
## Error: could not find function "ggplot"
```

```r
p <- p + geom_density() + xlab(param)
```

```
## Error: object 'p' not found
```

```r
p <- p + ggtitle(paste(param, "histogram"))
```

```
## Error: object 'p' not found
```

```r
print(p)
```

```
## Error: object 'p' not found
```


Follow the convergence of a given parameter. The vertical line indicates the burn-in point:


```r
thinning <- 200  # Take every kth samples
s <- seq(1, nrow(samples), thinning)
```

```
## Error: object 'samples' not found
```

```r
plot(s, samples[s, param], type = "l", main = paste(param, "convergence"))
```

```
## Error: object 's' not found
```

```r
abline(v = burnin)
```

```
## Error: plot.new has not been called yet
```



### Metacommunity

Metacommunity distribution with the local assembly model


```r
boxplot(metacommunity$p, las = 1)
```

```
## Error: object 'metacommunity' not found
```


Metacommunity averages:


```r
print(nmgs_metapopulation_average(metacommunity)$local)
```

```
## Error: object 'metacommunity' not found
```


_Metacommunity distribution for the full neutral model (q) - are the species
indices directly comparable, or how to combine across MCMC samples?_


### Local communities

Generated local communities are not currently available from the
output files but we can compare diversity distributions across MCMC
models w.r.t. observed data. The vertical line indicates diversity in
the observed data.


```r
library(ggplot2)
theme_set(theme_bw(15))

full <- qplot(stats$HN, binwidth = 0.1, geom = "histogram") + geom_vline(x = unique(stats$HO), 
    linetype = 2)
```

```
## Error: object 'stats' not found
```

```r
full <- full + ggtitle("Full neutral model; model diversity vs. observed diversity")
```

```
## Error: object 'full' not found
```

```r
full <- full + xlab("Diversity")
```

```
## Error: object 'full' not found
```

```r

local <- qplot(stats$HL, binwidth = 0.1, geom = "histogram") + geom_vline(x = unique(stats$HO), 
    linetype = 2)
```

```
## Error: object 'stats' not found
```

```r
local <- local + ggtitle("Local model; model diversity vs. observed diversity")
```

```
## Error: non-numeric argument to binary operator
```

```r
local <- local + xlab("Diversity")
```

```
## Error: non-numeric argument to binary operator
```

```r

library(gridExtra)
grid.arrange(local, full, nrow = 2)
```

```
## Error: object 'full' not found
```



Richnesses in MCMC simulations (histogram) and observed data (vertical
dashed line):


```r
full <- qplot(stats$SN, binwidth = 1, geom = "histogram") + geom_vline(x = unique(stats$SO), 
    linetype = 2)
```

```
## Error: object 'stats' not found
```

```r
full <- full + ggtitle("Full neutral model; model richness vs. observed richness")
```

```
## Error: object 'full' not found
```

```r
full <- full + xlab("Richness")
```

```
## Error: object 'full' not found
```

```r

local <- qplot(stats$SL, binwidth = 1, geom = "histogram") + geom_vline(x = unique(stats$SO), 
    linetype = 2)
```

```
## Error: object 'stats' not found
```

```r
local <- local + ggtitle("Local model; model richness vs. observed richness")
```

```
## Error: non-numeric argument to binary operator
```

```r
local <- local + xlab("Richness")
```

```
## Error: non-numeric argument to binary operator
```

```r

grid.arrange(local, full, nrow = 2)
```

```
## Error: object 'full' not found
```


Add here the similar histogram for the likelihoods

### Testing neutrality

Neutrality under the full neutral model and the local assembly
model. Compares the log-likelihood under the full neutral model to the
observed data:


```r
print(nmgs_neutrality(stats, "full"))
```

```
## Error: object 'stats' not found
```

```r
print(nmgs_neutrality(stats, "local"))
```

```
## Error: object 'stats' not found
```



### Licensing and Citations

This work can be freely used, modified and distributed under the [GNU
General Public GPL>=2
license](https://en.wikipedia.org/wiki/GNU_General_Public_License).

Kindly cite the work, if appropriate, as 'Leo Lahti, Christoper Quince
et al. (2013). NMGS R package. URL:
https://github.com/microbiome/NMGS)'.


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
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
## [1] gridExtra_0.9.1 ggplot2_0.9.3.1 NMGS_0.1.02     knitr_1.2      
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-2   dichromat_2.0-0    digest_0.6.3      
##  [4] evaluate_0.4.3     formatR_0.7        gtable_0.1.2      
##  [7] labeling_0.1       MASS_7.3-26        munsell_0.4       
## [10] plyr_1.8           proto_0.3-10       RColorBrewer_1.0-5
## [13] reshape2_1.2.2     scales_0.2.3       stringr_0.6.2     
## [16] tools_3.0.1
```





