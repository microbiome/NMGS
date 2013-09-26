#!/usr/bin/perl

my @values1 = ();
my @values2 = ();

my $val1 = shift(@ARGV);
my $val2 = shift(@ARGV);

for $file(@ARGV){

    open(FILE, $file) or die "Can't open $file\n";

    while($line = <FILE>){
	chomp($line);

	@tokens = split(/,/,$line);

	#$time = shift(@tokens);
	if($tokens[$val1] ne "nan"){
		push(@values1, $tokens[$val1]);
	}
	if($tokens[$val2] ne "nan"){
		push(@values2, $tokens[$val2]);
	}
	
    }

    close(FILE);
}

my $median1 = &median(0.5, @values1);
my $median2 = &median(0.5, @values2);

my $upper = 0;
my $total = 0;

foreach $val(@values1){
  if($val > $median2){
    $upper++;
  }
  $total++;
} 

printf("%f %f %d %d %f\n", $median1, $median2, $upper, $total, $upper/$total);

sub median()
{
    my $quantile = shift(@_);
    my @array = @_;
    
    @array = sort {$a <=> $b } @array;

    #print "@array\n";

    my $size = scalar(@array);

    return $array[$size*$quantile];
}
