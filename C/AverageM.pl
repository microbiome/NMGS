#!/opt/local/bin/perl

my @vals1 = ();
my @vals2 = ();
my @counts2 = ();

my $ur1 = 0.0;
my $ur2 = 0.0;

while($line = <STDIN>){
  $count++;

  my @tokens = split(/,/,$line);
  my $nS1 = $tokens[1];
  my $nS2 = $tokens[2];

  $line = <STDIN>;
  chomp($line);
  my @tokens = split(/,/,$line);
  $ur1 += pop(@tokens);
  @tokens = sort {$b <=>$a} @tokens;
  for($i = 0; $i < $nS1; $i++){
    $vals1[$i] += $tokens[$i];
  }

  $line = <STDIN>;
  chomp($line);
  my @tokens = split(/,/,$line);
  $ur2 += pop(@tokens);
  @tokens = sort {$b <=>$a} @tokens;
  for($i = 0; $i < $nS2; $i++){
    $vals2[$i] += $tokens[$i];
    $counts2[$i]++;
  }
}

my $nS1 = scalar(@vals1);

for($i = 0; $i < $nS1; $i++){
  printf("%d %f\n",$i, $vals1[$i]/$count);
}
printf("%d %f\n\n",$nS1 + 1, $ur1/$count);

my $nS2 = scalar(@vals2);
for($i = 0; $i < $nS2; $i++){
  printf("%d %f\n",$i,$vals2[$i]/$counts2[$i]);
}
printf("%d %f\n\n",$nS2 + 1, $ur2/$count);
