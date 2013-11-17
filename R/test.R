library(NMGS)

source("~/Rpackages/NMGS/R/pkg/R/io.R")
source("~/Rpackages/NMGS/R/pkg/R/analysis.R")

samples <- read_nmgs("~/Rpackages/NMGS/C/Simulation_out.csv")
metacommunity <- read_nmgs_metacommunity("~/Rpackages/NMGS/C/Simulation_out_m.csv")
stats <- read_nmgs_stats("~/Rpackages/NMGS/C/Simulation_out_s.csv")

# Calculate medians for all samples
medians <- nmgs_median(samples, burnin = 25000)
print(medians)

# Test neutrality
#' Compare the log-likelihood under the full neutral model with the observed data.
print(nmgs_neutrality(stats, "LN", "LO"))
#' the same but for local community assembly. 
print(nmgs_neutrality(stats, "LL", "LO"))

# Calculate metapopulation average (currently only for the local communities)
print(nmgs_metapopulation_average(metacommunity)$local)


# --------------------------------------------------

#dyn.load("/usr/lib/libgsl.so", local = FALSE, now = FALSE)
#dyn.load("/usr/lib/libgslcblas.so", local = FALSE, now = FALSE)
#dyn.load("pkg/src/NMGS.so")
#testf(1,2)
#.Call("compare_doubles", 1, 2, PACKAGE = "NMGS")

