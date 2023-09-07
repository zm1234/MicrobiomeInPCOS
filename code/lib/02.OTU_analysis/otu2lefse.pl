#! /usr/bin/perl
# Function: transfer otu_table.even.txt to the file lefse can accepted
# Data: 2014-06-12, Version 1.0, input is otu_table.relative.mat
# Date: 2015-01-12, fix the bug for lefse input,the script has been modified
# Date: 2015-01-12, Version 2.0, input is Relative mats folder
# Contact: chenjunru@novogene.cn

use strict;
use warnings;
use Getopt::Long;
my %opt = (category=>1);
GetOptions(\%opt,"category:i");
if (@ARGV < 3) {
	die "usage: perl $0 <input:Relative mats folder> <input:all.mf> <output: for lefse input> [vs_group:group1 group2] [-c 1]
--category  set Description col num based on mf file,default is 1\n\n";
}
my($otu_table,$mf,$output,@vs_group)=@ARGV;
my (%mf,@vs_samples);
my $mf_col = $opt{category};

open(MF,"$mf");
#<MF>;
my $sample_num;
while (my $or=<MF>) {
	chomp $or;
	my @or=split(/\s+/,$or);
	$mf{$or[0]}=$or[$mf_col];
    $sample_num=`grep -w  $or[$mf_col] $mf | wc -l`;
	if(@vs_group && $sample_num > 2) { #添加组内样本个数的判断，9月10日
		for (@vs_group){
			if ($or[$mf_col] eq $_) {
				push @vs_samples,$or[0];
				#print "@vs_samples\n";
				# body...
			}

		}

		#push my @vs_samples,$or[0];
		# body...
	}
}
close MF;

my (%tax2mat,@lastsamples);
foreach(`ls $otu_table/tax_table.[cfgop].txt`){
	chomp;
	open(OR,"$_");
	my $head=<OR>;
	chomp $head;
	my @samples=split/\t/,$head;
	@lastsamples=@samples;
	while (my $or=<OR>) {
		chomp $or;
		my @or=split/\t/,$or;
		next if($or[0] eq 'Others' or $or[0] eq 'Other');
		$or[-1]=~s/^k__//;
		$or[-1]=~s/;[cfgops]__/|/g;
		$or[-1]=~s/;$//;
		$or[-1]=~s/\s+/_/g;
		for (my $i = 1; $i < $#or; $i++) {
			$tax2mat{$or[-1]}{$samples[$i]}=$or[$i];
		}
	}
	close OR;
}



if (@vs_group) {
	print "lefse_vs_group:@vs_group\ngroup_samples:@vs_samples\n";
	open(OUT,">$output");
    print OUT "Class";
    for (my $i = 0; $i <=$#vs_samples; $i++) {
		print OUT "\t$mf{$vs_samples[$i]}";
    }
    print OUT "\n";
    foreach(sort{$a cmp $b} keys %tax2mat){
	    my $key=$_;
        my $sum_abundance=0;
        for (my $i=0; $i <= $#vs_samples;$i++){
            $sum_abundance+=$tax2mat{$key}{$vs_samples[$i]};
        }
        next if($sum_abundance==0);
	    print  OUT "$key";
	    for (my $i = 0; $i <= $#vs_samples; $i++) {
		    print OUT "\t$tax2mat{$key}{$vs_samples[$i]}";
	    }
	    print OUT "\n";
    }
    close OUT;
	# body...
} else {
	open(OUT,">$output");
    print OUT "Class";
    for (my $i = 1; $i < $#lastsamples; $i++) {
		print OUT "\t$mf{$lastsamples[$i]}";
    }
    print OUT "\n";
    foreach(sort{$a cmp $b} keys %tax2mat){
	    my $key=$_;
        my $sum_abundance=0;
        for (my $i=1; $i < $#lastsamples;$i++){
            $sum_abundance+=$tax2mat{$key}{$lastsamples[$i]};
        }
        next if($sum_abundance==0);
        print  OUT "$key";
	    for (my $i = 1; $i < $#lastsamples; $i++) {
		    print OUT "\t$tax2mat{$key}{$lastsamples[$i]}";
	    }
	    print OUT "\n";
    }
    close OUT;
	# else...
}

