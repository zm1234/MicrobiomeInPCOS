#Junru Chen, chenjr@geneworks.cn
##2016-10-21
use strict;
use warnings;
use FindBin qw($Bin);


die"usage:perl $0 <in:biom> <in:summarize> <output> <mark:min|max|median|mean> <sys:single|multiple>\n" unless @ARGV == 5;

my($in,$summarize,$out,$mark,$sys)=@ARGV;
$mark && $mark=~/^(min|max|median|mean)$/ || die "mark can just be min|max|median|mean\n";
$sys && $sys=~/^(single|multiple)$/ || die "mark can just be single|multiple\n";


#get software pathway
use lib "$Bin/../../lib/00.Commbin";
use PATHWAY;
(-s "$Bin/../../bin/Pathway_cfg.txt") || die"error: can't find config at $Bin/../../bin/, $!\n";
my($qiime,)=get_pathway("$Bin/../../bin/Pathway_cfg.txt",[qw(QIIME)]);

my%mark;
my@or=split/\s+/,`grep Min: $summarize`;
$mark{min}=sprintf("%d",(split/\s+/,`grep Min: $summarize`)[2]);
$mark{max}=sprintf("%d",(split/\s+/,`grep Max: $summarize`)[2]);
$mark{median}=sprintf("%d",(split/\s+/,`grep Median: $summarize`)[2]);
$mark{mean}=sprintf("%d",(split/\s+/,`grep Mean: $summarize`)[2]);

my$depth=$mark{$mark};
my$step=sprintf("%d",$depth/10);

if ($sys eq 'single') {
	print "#single_rarefaction\n$qiime/single_rarefaction.py -i $in -o $out -d 20000 \n";
	system("$qiime/single_rarefaction.py -i $in -o $out -d 20000 ") ;
}elsif($sys eq 'multiple' ){
	print("#multiple_rarefactions\n$qiime/multiple_rarefactions.py -i $in -m 10 -x $depth -s $step -o $out\n") ;
	system("$qiime/multiple_rarefactions.py -i $in -m 10 -x $depth -s $step -o $out ") ;
}
