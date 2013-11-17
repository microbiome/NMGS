#!/usr/bin/perl

@values = ();
@theta  = ();

$nBurn = shift(@ARGV);
$size = 0;
my $count = 0;
for $file(@ARGV){

    open(FILE, $file) or die "Can't open $file\n";


    while($line = <FILE>){
	chomp($line);

	@tokens = split(/,/,$line);

	$time = shift(@tokens);

	my $tt = shift(@tokens);
	if($count == 0){
	    $size = scalar(@tokens);
	    for($i = 0; $i < $size; $i++){
	       my @temp = ();
               $values[$i] = \@temp; 
            }
	
	}

	if($time > $nBurn){
	    push(@theta,$tt);
	    for($i = 0; $i < $size; $i++){
	    	push(@{$values[$i]}, $tokens[$i]);
	    }
	    $count++;
	}
	
    }

    close(FILE);
}

my @meds = ();
$x = &median(0.5, @theta);
$xl = &median(0.025, @theta);
$xu = &median(0.975, @theta);
printf("%f,%f,%f\n",$xl,$x,$xu);
for($i = 0; $i < $size; $i++){
#	push(@meds,&median(0.5, @{$values[$i]}));
	printf("%d,%f,%f,%f\n", $i,&median(0.025, @{$values[$i]}), &median(0.5, @{$values[$i]}), &median(0.975, @{$values[$i]}));
}

for($i = 0; $i < $size; $i++){
	#print "$x $meds[$i]\n";
}

$yl = &median(0.025, @meds);
$y = &median(0.5, @meds);
$yu = &median(0.975, @meds);

#printf("%f %f %f %f %f %f\n",$x,$y,$x - $xl, $xu - $x, $y - $yl, $yu - $y);

sub median()
{
    my $quantile = shift(@_);
    my @array = @_;
    
    @array = sort {$a <=> $b } @array;

    #print "@array\n";

    my $size = scalar(@array);

    return $array[$size*$quantile];
}
