use strict;
use warnings;

die"usage:perl $0 <in:txt> <in:col> <out>\n" unless @ARGV == 3;
my($in,$col,$out)=@ARGV;

my@col=("sample","OBS","ACE","Chao1","Simpson","Shannon","Coverage");
my%drug=(
"d1","d1",
"d2","d2",
"d3","d3",
"dc","Control",
"df","Model",
);
my%time=(
"L0","0",
"L3","3",
"L7","7",
"L14","14",
"L28","28",
"0","0",
"3","3",
"7","7",
"14","14",
"28","28",
);


open OUT,">$out";
print OUT "Drug\tTime\tSample\t".$col[$col-1]."\n";
open IN,$in;
<IN>;
while (<IN>) {
	chomp;
	my($sample,$l)=(split/\t/)[0,$col-1];
	my($drug,$time)=($1,$2) if $sample=~/(d[123cf])(L?\d+)\d/;
#	next if $drug eq "d1";
	$drug=$drug{$drug};
	$time=$time{$time};
	print OUT "$drug\t$time\t$sample\t$l\n";
}
close IN;
close OUT;
