use strict;
use warnings;
die"usage:perl $0 <in> <out>\n" unless @ARGV == 2;
my($in,$out)=@ARGV;

open(OR,$in);
open(OUT,">$out");
while(<OR>){
	chomp;
	my@or=split/\t/;
	print OUT "$or[0]\t$or[3]\n";
}
close OUT;
close OR;
