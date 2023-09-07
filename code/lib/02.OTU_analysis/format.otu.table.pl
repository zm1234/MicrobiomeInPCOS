use strict;
use warnings;
die"usage:perl $0 <in> <out>\n" unless @ARGV == 2;
my($in,$out)=@ARGV;

open(OR,$in);
open(OUT,">$out");
while(<OR>){
	chomp;
	next if /^# Constructed from biom file/;
	$_=~s/^#//;
	my@or=split/\t/;pop@or;
	print OUT join("\t",@or)."\n";
}
close OUT;
close OR;
