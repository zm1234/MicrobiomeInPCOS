use strict;
use warnings;

die"usage:perl $0 <indir> <in:mf> <outdir>\n" unless @ARGV ==3;
my($indir,$mf,$outdir)=@ARGV;
-s $outdir || `mkdir -p $outdir`;

my (%s2g,@groups);
open OR,$mf;
while(<OR>){
	chomp;
	next if /^#/;
	my@or=split/\s+/;
	$s2g{$or[0]}=$or[1];
	push @groups,$or[1] if ! grep{$or[1] eq $_ } @groups;
}
close OR;

opendir INDIR,$indir;
while(my $file=readdir(INDIR)){
	next if $file eq '.' || $file eq '..';
	next unless $file=~/tax_table.[kpcofgs].txt/;
	my $input = "$indir/$file";
	open IN,$input;
	open OUT,">$outdir/$file";
	print OUT "Tax\t".join("\t",@groups)."\tTaxonomy\n";
	my $head=<IN>;chomp$head;my@head=split/\t/,$head;
	my %g2other;
	while(<IN>){
		chomp;
		my@or=split/\t/;
		next if $or[0] eq 'Other';
		my %g2abun;
		for(my $i=1;$i < @or-1;$i++){
			$g2abun{$s2g{$head[$i]}}{'abun'} += $or[$i];
			$g2abun{$s2g{$head[$i]}}{'num'} ++;
		}
		print OUT "$or[0]";
		foreach my $g (@groups){
			if ($g2abun{$g}{'abun'}) {
				print OUT "\t".$g2abun{$g}{'abun'}/$g2abun{$g}{num};
				$g2other{$g} += $g2abun{$g}{'abun'}/$g2abun{$g}{num};
			}else{
				print OUT "\t0";
			}			
		}
		print OUT "\t$or[-1]\n";
	}
	print OUT "Other";
	foreach my $g (@groups){
		my $others= 1- $g2other{$g};
		$others < 0 ?
		print OUT "\t0" :
		print OUT "\t$others";
	}
	print OUT "\tOther\n";
	close OUT;
	close IN;
}
closedir INDIR;
