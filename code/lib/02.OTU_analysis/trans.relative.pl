#Junru Chen, chenjr@geneworks.cn
##2016-10-21
use strict;
use warnings;

die"usage:perl $0 <indir> <outdir> <outprefix>\n" unless @ARGV ==3;
my($indir,$outdir,$prefix)=@ARGV;
opendir(INDIR,"$indir");
(-s $outdir) || `mkdir -p $outdir`;
while (readdir(INDIR)) {
	my$file=$_;
	next if $file eq '.' || $file eq '..';
	$file = "$indir/$file";
	my@ranks=qw(k p c o f g s);
	my$rank=$1 if $file=~/([23456]).txt/;
	$rank=$ranks[$rank-1];
	my$output = "$outdir/$prefix.$rank.txt";
	&trans_file($file,$output);
}
closedir INDIR;


sub trans_file{
	my($infile,$outfile)=@_;
	my@ranks=qw(k p c o f g s);
	my$rank=$1 if $infile=~/([23456]).txt/;
	$rank=$ranks[$rank-1];
	open(OR,$infile);
	<OR>;
	my $head=<OR>;chomp$head;
	$head=~s/^#OTU ID/Tax/;
	my @head=split/\t/,$head;
	my@others;$others[0]='Other';
	my%max2abun;
	my %check_double_name; #add songdq 
	while (<OR>) {
		chomp;
		my@or=split/\t/;
		if ($or[0]=~/${rank}__$|;Other$/) {
			for (my $i = 1; $i < @or; $i++) {
				$others[$i] += $or[$i];
			}
		}else{
			my $detail_tax=shift@or;
			my @taxes=split/;/,$detail_tax;
			$taxes[-1]=~s/^[kpcofgs]__//;
			my $max=(sort{$a <=> $b} @or)[-1];
			#add songdq in 20190509 
			if(exists $check_double_name{$taxes[-1]}){
				#$check_double_name{$taxes[-1]} ++;
				$taxes[-1]="$taxes[-1]_$check_double_name{$taxes[-1]}";
				$check_double_name{$taxes[-1]} ++;
			}else{
				$check_double_name{$taxes[-1]}++;
			}
			my $new_line="$taxes[-1]\t".join("\t",@or)."\t$detail_tax\n";
			push @{$max2abun{$max}},$new_line;
		}
	}
	close OR;
	open(OUT,">$outfile");
	print OUT "$head\tTaxonomy\n";
	foreach my $max (sort{$b <=> $a} keys %max2abun){
		print OUT join("",@{$max2abun{$max}});
	}
	my $ar_nm=@others;
	print OUT join("\t",@others)."\tOther\n" if $ar_nm != 1;
	close OUT;
}
