use strict;
use warnings;
die"usage:perl $0 <in> <out>\n" unless @ARGV == 2;
my($in,$out)=@ARGV;

open(OR,$in);
open(OUT,">$out");
while(<OR>){
	chomp;
	next if /^#/;
	my@or=split/\s+/;
	print OUT "$or[0]\t$or[1]\n";
}
close OUT;
close OR;
