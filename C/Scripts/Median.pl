#!/usr/bin/perl

@values = ();

$var   = shift(@ARGV);

$nBurn = shift(@ARGV);

for $file(@ARGV){

    open(FILE, $file) or die "Can't open $file\n";

    $line = <FILE>;

    while($line = <FILE>){
	chomp($line);

	@tokens = split(/,/,$line);

	$time = $tokens[0];

	$value = $tokens[$var];
    
	if($time > $nBurn){
	    push(@values, $value);
	}
	
    }

    close(FILE);
}

printf("%f %f %f\n", &median(0.025, @values), &median(0.5, @values), 
       &median(0.975, @values));

sub median()
{
    my $quantile = shift(@_);
    my @array = @_;
    
    @array = sort {$a <=> $b } @array;

    #print "@array\n";

    my $size = scalar(@array);

    return $array[$size*$quantile];
}
