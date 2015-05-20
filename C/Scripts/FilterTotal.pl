#!/usr/bin/perl

use strict;
my $minRare = $ARGV[0];
my $line = <STDIN>;
chomp($line);
my @tokens = split(/,/,$line);
shift(@tokens);

my @samples  = @tokens;
my $nSamples = scalar(@samples);
my @values   = ();
my @OTUs     = ();
my @totals   = ();

my $i = 0;
while($line = <STDIN>){
  chomp($line);

  @tokens = split(/,/,$line);

  push(@OTUs, shift(@tokens));

  my $j = 0;
  foreach my $tok(@tokens){
    $values[$j][$i] = $tok;
    $totals[$j] += $tok;
    $j++;
  }

  $i++;
}
my $nOTUs = scalar(@OTUs);


printf("OTU,");
my @tsamples = ();
for($i = 0; $i < $nSamples; $i++){
	if($totals[$i] > $minRare){
		push(@tsamples,$samples[$i]);
	}
}

my $stringName = join(",",@tsamples);
print "$stringName\n";

for(my $j = 0; $j < $nOTUs; $j++){
  #print "$totals[$j] $j\n";

    my $i = 0;
    printf("%s,",$OTUs[$j]);
    my @tvalues = ();
    for($i = 0; $i < $nSamples; $i++){
	if($totals[$i] > $minRare){
		push(@tvalues,$values[$i][$j]);
    	}
    }	 
    my $stringName = join(",",@tvalues);
   print "$stringName\n";	
}
