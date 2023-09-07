use strict;
use warnings;

die "usage:perl $0 <indir> <outdir>\n" unless @ARGV == 2;
my($indir,$outdir)=@ARGV;
-s $outdir || `mkdir -p $outdir`;

opendir DIR,$indir;
while(my $in=readdir(DIR)){
next if $in eq '.' || $in eq '..';
open(OR,"$indir/$in");
my $head=<OR>;
my%max2line;my@others;
while(<OR>){
	chomp;
	my@or=split/\t/;
	my$tax=shift@or;
	my$detail=pop@or;
	my $max=(sort{$a <=> $b} @or)[-1];
	$tax eq "Other" ? 
	push @others,"$_\n":
	push @{$max2line{$max}},"$_\n";
}
close OR;

open(OUT,">$outdir/$in");
print OUT "$head";
foreach my $max (sort{$b <=> $a} keys %max2line){
	print OUT join("",@{$max2line{$max}});
}
print OUT join("",@others) if @others;
close OUT;
}
