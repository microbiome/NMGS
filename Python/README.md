### Neutral model simulation

To simulate the neutral model, run:

python ./HDPSample.py -t 10.0 -n 10 -o Test.csv -s 1 -x Test_x.csv -m 10.0

 * t theta the biodiversity parameter  
 * n is the number of sites
 * o the output frequency file of species across site
 * s random seed
 * x output for a debug file, this is the ancestors don’t worry about it
 * m is the migration rate - fixed for all sites

