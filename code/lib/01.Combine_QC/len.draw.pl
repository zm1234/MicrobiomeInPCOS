use strict;
use warnings;

die "usage:perl $0 <in:fa> <outdir> <out:prefix>\n" unless @ARGV == 3;
my($fa,$outdir,$prefix)=@ARGV;
(-s $outdir) || `mkdir -p $outdir`;

$/="\n>";
open(OR,$fa);
my%len2num;
while(my$seq=<OR>){
	chomp$seq;
	$seq=~s/^>//;
	my$head=$1 if $seq=~/^(.*)/;
	$seq=~s/^.*//;
	$seq=~s/\s+//g;
	my $len=length($seq);
	$len2num{$len}++;
}
close OR;

open(LEN,">$outdir/$prefix.txt");
print LEN "Length\tNumber\n";
foreach (sort{$a <=> $b} keys %len2num){
	print LEN "$_\t$len2num{$_}\n";
}
close LEN;

open(R,">$outdir/$prefix.R");
print R 
"setwd(\"$outdir\")
data <- read.table(\"$prefix.txt\",sep=\"\\t\",header=T)
pdf(\"$prefix.pdf\")
plot(data,type=\"h\",lwd = 2,col=\"firebrick\",xlab=\"Length(nt)\",ylab=\"Number(#)\")
dev.off()";
close R;
system("Rscript $outdir/$prefix.R");
system("cd $outdir;convert -density 300 $prefix.pdf $prefix.png");
