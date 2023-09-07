#http://blog.csdn.net/anshiquanshu/article/details/68954536
use strict;
use warnings;

die"usage:perl $0 <in:table> <in:group> <output:prefix> <opts:cutoff> \n" unless @ARGV >=3;
my($table,$group_in,$output)=@ARGV;
my$cutoff=$ARGV[3] if $ARGV[3]; 

my %s2group;
my (@groups,@samples);
open OR,$group_in;
while (<OR>) {
	chomp;
	next if /^#/;
	my@or=split/\s+/;
	$s2group{$or[0]}=$or[1];
	push @groups,$or[1] if !grep{$or[1] eq $_ } @groups;
	push @samples,$or[0];
}
close OR;

open OUT2,">$output.median(avg).txt";
my @head_groups2=map{$_.".(median)(avg)";} @groups;
print OUT2 "ID\t".join("\t",@head_groups2),"\n";

my (%s2abun,%out2out,%out2avg,%out2table,%out2box,@tax_p,@pvalue_p);
open OR,"$table";
my $head=<OR>;chomp$head;
my @head=split/\t/,$head;pop@head;
while (my $l=<OR>) {
	chomp$l;
	my @or=split/\t/,$l;
	next if $or[0] eq 'Other';
	my (%g2abun);
	my $nzeronum=0;
	for (my $i = 1; $i < @or -1; $i++){
		next unless $s2group{$head[$i]};
		if ($or[$i]) {
			push @{$g2abun{$s2group{$head[$i]}}},$or[$i];
			$s2abun{$or[0]}{$head[$i]}=$or[$i];
			$nzeronum++ if $or[$i] ne '0.0';
		}else{
			push @{$g2abun{$s2group{$head[$i]}}},"0";
		}
	}
	next if $cutoff && $nzeronum < $cutoff;
	my (@group,@data,$flag);
	$out2out{$or[0]}="$or[0]";
	$out2avg{$or[0]}="$or[0]";
	print OUT2 "$or[0]";
	foreach my $g (@groups){
		push @group,($g) x @{$g2abun{$g}};
		push @data,@{$g2abun{$g}};
		my $median=&median(@{$g2abun{$g}});
		my $avg=&avg(@{$g2abun{$g}});
		$median=sprintf("%.2e", $median);
		$avg=sprintf("%.2e", $avg);
		$median=0 if $median eq "0.00e+00";
		$avg=0 if $avg eq "0.00e+00";
		print OUT2 "\t$median($avg)";
	}
	print OUT2 "\n";
}
close OR; 


sub avg{
	my($num,$total);
	foreach my $n (@_){
		$num++;
		$total += $n;
	}
	return $total/$num;
}

sub median{
    my @list = sort{$a<=>$b} @_;
    my $count = @list;
    if( $count == 0 )
    {
        return undef;
    }   
    if(($count%2)==1){
        return $list[int(($count-1)/2)];
    }
    elsif(($count%2)==0){
        return ($list[int(($count-1)/2)]+$list[int(($count)/2)])/2;
    }
}
