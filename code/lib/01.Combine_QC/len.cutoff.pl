use strict;
use warnings;

die "usage:perl $0 <in:fa> <output> <min> <max>\n" unless @ARGV == 4;
my($fa,$output,$min,$max)=@ARGV;

$/="\n>";
open(OR,$fa);
open OUT,">$output";
while(my$seq=<OR>){
	chomp$seq;
	$seq=~s/^>//;
	my$head=$1 if $seq=~/^(.*)/;
	$seq=~s/^.*//;
	$seq=~s/\s+//g;
	my $len=length($seq);
	print OUT ">$head\n$seq\n" if $len > $min && $len  < $max;
}
close OR;
close OUT;
