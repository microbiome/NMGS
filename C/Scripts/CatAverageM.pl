#!/opt/local/bin/perl

my $csvFile = $ARGV[0];
my $metFile = $ARGV[1];

my @vals1 = ();
my @vals2 = ();
my @counts2 = ();
my $count = 0;
my $ur1 = 0.0;
my $ur2 = 0.0;

open(FILE, $metFile) or die "Can't open $metFile\n";

while($line = <FILE>){
  $count++;

  my @tokens = split(/,/,$line);
  my $nS1 = $tokens[1];
  my $nS2 = $tokens[2];

  $line = <FILE>;
  chomp($line);
  my @tokens = split(/,/,$line);
  $ur1 += pop(@tokens);
  
  for($i = 0; $i < $nS1; $i++){
    $vals1[$i] += $tokens[$i];
  }

  $line = <FILE>;
  chomp($line);
  my @tokens = split(/,/,$line);
  $ur2 += pop(@tokens);
  
  for($i = 0; $i < $nS2; $i++){
    $vals2[$i] += $tokens[$i];
    $counts2[$i]++;
  }
}

close(FILE);

my $nS1 = scalar(@vals1);

open(FILE, $csvFile) or die "Can't open $csvFile\n";

my $line = <FILE>;
chomp($line);
my @tokens = split(/,/,$line);
shift(@tokens);

my @samples  = @tokens;
my $nSamples = scalar(@samples);
my @values   = ();
my @OTUs     = ();
my @totals   = ();

my $i = 0;
while($line = <FILE>){
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
push(@tsamples,"M");
for($i = 0; $i < $nSamples; $i++){
  push(@tsamples,$samples[$i]);
}

my $stringName = join(",",@tsamples);
print "$stringName\n";

for(my $j = 0; $j < $nOTUs; $j++){

    my $i = 0;
    printf("%s,",$OTUs[$j]);
    my @tvalues = ();
    
    push(@tvalues,$vals1[$j]/$count);
    for($i = 0; $i < $nSamples; $i++){
	push(@tvalues,$values[$i][$j]);
    }	 

    my $stringName = join(",",@tvalues);
   print "$stringName\n";	
}


