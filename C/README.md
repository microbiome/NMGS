## Initializing the code

To make NMGS executable simply type:

make

In Ubuntu it might necessary to install the libraries gsl-bin and
libgsl0-dev first.


## Model fitting

As an example generate Gibbs samples from fitting to the simulated data set

./NMGS -in Simulation.csv -out Simulation_out -v -s

-in specifies the input abundances across samples
-out the output file stub
-v prints progress to screen
-s generates samples to determine whether community appears neutral and output sampled metacommunity relative abundances


### Output files

Produces three output files

Simulation_out.csv
Simulation_out_m.csv
Simulation_out_s.csv

Simulation_out.csv: Contains the MCMC sampled parameters in format:

iter1,theta,i1,..,iN
..
iterM,theta,i1,..,iN

where theta is the biodiversity parameter and i1,...,iN are the immigration rates to each sample

Simulation_out_m.csv: Gives the generated metacommunities with format:

SampleN,SL,SN
p1,...,pSL
q1,...,qSN

where SampleN is nth sample generated from the neutral model with parameters from the corresponding MCMC sample. SL and SN are the number of species in the local assembly sample (always S) and the sample generated under the full neutral model respectively, p and q are then the metacommunity distributions

Simulation_out_s.csv: Gives statistics on the sampled communities in format:

SampleN,LN,LL,LO,HN,HL,HO,SN,SL,SO

where SampleN is nth sample generated under the neutral model with fitted parameters taken from the corresponding MCMC sample. LN,LL,LO are the log-likelihoods of the full neutral sample, the local community sample and the observed sample respectively. Then HN,HL,HO and SN,SL,SO are the species entropies and richness's in the same order.


## Analyzing the output

Scripts are provided to process this output:

Median.pl generates medians from the sampled MCMC parameters. It takes the variable no as a first argument and the burn time as the second, the third is the sample file itself so:

./Scripts/Median.pl 1 25000 Simulation_out.csv 

Calculates the median theta values with upper and lower confidence intervals for a burn-in of 25,000 iterations. Format:

lc med uc

where lc is the lower 95% confidence interval, med is the median and uc the upper confidence interval. To calculate median for the first immigration rate:

./Scripts/Median.pl 2 25000 Simulation_out.csv

etc.

To determine if community appears neutral use the script Sig.pl this calculates proportion of one sample that exceeds the median of a second so:

./Scripts/Sig.pl 1 3 Simulation_out_s.csv

Compares the log-likelihood under the full neutral model with the observed data.

./Scripts/Sig.pl 2 3 Simulation_out_s.csv

Is the same but for local community assembly. The output format is:

med1 med2 n1 nT p

where med1 is the median for the first value, med2 second, n1 is the number of samples of the first variable that exceed the second, nT is the total number of samples and p the proportion that exceed i.e. what we take as the pseudo p-value.

Finally to average over the metapopulations we use AverageM.pl:

./Scripts/AverageM.pl < Simulation_out_m.csv

